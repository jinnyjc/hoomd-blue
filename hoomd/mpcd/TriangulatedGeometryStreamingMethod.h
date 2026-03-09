// Copyright (c) 2009-2025 The Regents of the University of Michigan.
// Part of HOOMD-blue, released under the BSD 3-Clause License.

/*!
 * \file mpcd/TriangulatedGeometryStreamingMethod.h
 * \brief Declaration of mpcd::TriangulatedGeometryStreamingMethod
 */

#ifndef MPCD_TRIANGULATED_GEOMETRY_STREAMING_METHOD_H_
#define MPCD_TRIANGULATED_GEOMETRY_STREAMING_METHOD_H_

#ifdef __HIPCC__
#error This header cannot be compiled by nvcc
#endif

#include "StreamingMethod.h"
#include "TriangulatedGeometry.h"
#include <pybind11/pybind11.h>

namespace hoomd
    {
namespace mpcd
    {
//! MPCD confined streaming method
/*!
 * This method implements the base version of ballistic propagation of MPCD
 * particles in confined triangulated geometries.
 *
 * \tparam Geometry The confining triangulated geometry.
 * \tparam Force The solvent force.
 *
 * The integration scheme is essentially Verlet with specular reflections. The particle is streamed
 * forward over the time interval. If it moves outside the Geometry, it is placed back on the
 * boundary and its velocity is updated according to the boundary conditions. Streaming then
 * continues until the timestep is completed.
 *
 */

inline bool intersectTriangle(const vec3<Scalar>& pos_v,
                              const vec3<Scalar>& vel_v,
                              const vec3<Scalar>& a,
                              const vec3<Scalar>& b,
                              const vec3<Scalar>& c,
                              const Scalar dt_remain,
                              const Scalar eps,
                              Scalar& t_hit)

    {
    // calculate dimension where the ray direction is maximal
    const Scalar ax = fabs(vel_v.x);
    const Scalar ay = fabs(vel_v.y);
    const Scalar az = fabs(vel_v.z);

    int kz = 0;
    if (ay > ax)
        kz = 1;
    if ((kz == 0 && az > ax) || (kz == 1 && az > ay))
        kz = 2;
    int kx = kz + 1;
    if (kx == 3)
        kx = 0;
    int ky = kx + 1;
    if (ky == 3)
        ky = 0;

    const Scalar det = vel_v[kz];
    if (fabs(det) <= eps)
        return false;

    if (det < Scalar(0.0))
        std::swap(kx, ky);

    // calculate shear constants
    const Scalar Sx = vel_v[kx] / det;
    const Scalar Sy = vel_v[ky] / det;
    const Scalar Sz = Scalar(1.0) / det;

    // calculate vertices relative to ray origin
    const vec3<Scalar> A = a - pos_v;
    const vec3<Scalar> B = b - pos_v;
    const vec3<Scalar> C = c - pos_v;

    // apply shear and scale
    const Scalar Ax = A[kx]- Sx * A[kz];
    const Scalar Ay = A[ky]- Sy * A[kz];
    const Scalar Bx = B[kx]- Sx * B[kz];
    const Scalar By = B[ky]- Sy * B[kz];
    const Scalar Cx = C[kx]- Sx * C[kz];
    const Scalar Cy = C[ky]- Sy * C[kz];

    // calculate scaled barycentric coordinates
    Scalar u = Cx * By - Cy * Bx;
    Scalar v = Ax * Cy - Ay * Cx;
    Scalar w = Bx * Ay - By * Ax;

    if ((u < Scalar(0.0) || v < Scalar(0.0) || w < Scalar(0.0)) 
            && (u > Scalar (0.0) || v > Scalar(0.0) || w > Scalar(0.0)))
        return false;

    const Scalar inv_det = Scalar(1.0) / (u + v + w);
    if (!std::isfinite((double)inv_det))
        return false;

    // scaled z
    const Scalar Az = Sz * A[kz];
    const Scalar Bz = Sz * B[kz];
    const Scalar Cz = Sz * C[kz];

    const Scalar t = (u * Az + v * Bz + w * Cz) * inv_det;

    if (t <= eps || t > dt_remain)
        return false;

    t_hit = t;
    return true;
    }


template<class Force>
class PYBIND11_EXPORT TriangulatedGeometryStreamingMethod : public mpcd::StreamingMethod
    {
    public:
    //! Constructor
    /*!
     * \param sysdef System definition
     * \param cur_timestep Current system timestep
     * \param period Number of timesteps between collisions
     * \param phase Phase shift for periodic updates
     * \param geom Streaming geometry
     * \param force Solvent force
     */
    TriangulatedGeometryStreamingMethod(std::shared_ptr<SystemDefinition> sysdef,
                              unsigned int cur_timestep,
                              unsigned int period,
                              int phase,
                              std::shared_ptr<TriangulatedGeometry> geom,
                              std::shared_ptr<Force> force)
        : mpcd::StreamingMethod(sysdef, cur_timestep, period, phase), m_geom(geom), m_force(force)
        {
        }

    //! Implementation of the streaming rule
    void stream(uint64_t timestep) override;

    //! Get triangulated geometry
    std::shared_ptr<TriangulatedGeometry> getGeometry() const
        {
        return m_geom;
        }

    //! Set triangulated geometry
    void setGeometry(std::shared_ptr<TriangulatedGeometry> geom)
        {
        m_geom = geom;
        }

    //! Set the solvent force
    std::shared_ptr<Force> getForce() const
        {
        return m_force;
        }

    //! Get the solvent force
    void setForce(std::shared_ptr<Force> force)
        {
        m_force = force;
        }

    protected:
    std::shared_ptr<TriangulatedGeometry> m_geom; //!< Triangulated geometry
    std::shared_ptr<Force> m_force;   //!< Solvent force
    };

/*!
 * \param timestep Current time to stream
 */
template<class Force>
void TriangulatedGeometryStreamingMethod<Force>::stream(uint64_t timestep)
    {
    if (!shouldStream(timestep))
        return;

    if (!m_cl)
        {
        throw std::runtime_error("Cell list has not been set");
        }

    if (!m_geom)
        {
        throw std::runtime_error("Triangulated geometry has not been set");
        }

    const BoxDim box = m_cl->getCoverageBox();

    ArrayHandle<Scalar4> h_pos(m_mpcd_pdata->getPositions(),
                               access_location::host,
                               access_mode::readwrite);
    ArrayHandle<Scalar4> h_vel(m_mpcd_pdata->getVelocities(),
                               access_location::host,
                               access_mode::readwrite);
    ArrayHandle<Scalar3> h_vertices(m_geom->getUnwrappedVertices(),
                                    access_location::host,
                                    access_mode::read);
    ArrayHandle<uint3> h_triangles(m_geom->getUnwrappedTriangles(),
                                   access_location::host,
                                   access_mode::read);
    const Scalar mass = m_mpcd_pdata->getMass();

    // default construct a force if one is not set
    const Force force = (m_force) ? *m_force : Force();

    for (unsigned int cur_p = 0; cur_p < m_mpcd_pdata->getN(); ++cur_p)
        {
        const Scalar4 postype = h_pos.data[cur_p];
        Scalar3 pos = make_scalar3(postype.x, postype.y, postype.z);
        const unsigned int type = __scalar_as_int(postype.w);

        const Scalar4 vel_cell = h_vel.data[cur_p];
        Scalar3 vel = make_scalar3(vel_cell.x, vel_cell.y, vel_cell.z);
        // estimate next velocity based on current acceleration
        vel += Scalar(0.5) * m_mpcd_dt * force.evaluate(pos) / mass;

        // propagate the particle to its new position ballistically
        Scalar dt_remain = m_mpcd_dt;
        bool collide = true;

        vec3<Scalar> pos_v(pos);
        vec3<Scalar> vel_v(vel);
        do
            {
            bool found = false;

            // to keep the earliest hit
            Scalar best_t = dt_remain;
            vec3<Scalar> best_n(0, 1, 0);
            vec3<Scalar> best_pos = pos_v;

            const unsigned int num_triangles = m_geom->getNumUnwrappedTriangles();
            const Scalar eps = Scalar(1e-8);

            for (unsigned int cur_tri = 0; cur_tri < num_triangles; ++cur_tri)
                {
                const uint3 triangles = h_triangles.data[cur_tri];

                const vec3<Scalar> a(h_vertices.data[triangles.x]);
                const vec3<Scalar> b(h_vertices.data[triangles.y]);
                const vec3<Scalar> c(h_vertices.data[triangles.z]);

                // find intersection
                Scalar t_hit;
                if (!intersectTriangle(pos_v, vel_v, a, b, c, dt_remain, eps, t_hit))
                    continue;

                if (t_hit >= best_t)
                    continue;

                // compute triangle normal
                const vec3<Scalar> e1 = b - a;
                const vec3<Scalar> e2 = c - a;
                const vec3<Scalar> n = cross(e1, e2);
                const Scalar nn = dot(n, n);

                // degenerate triangles
                if (nn <= Scalar(0))
                    continue;

                vec3<Scalar> n_unit = n * (Scalar(1.0) / fast::sqrt(nn));

                // adjust the normal against velocity
                if (dot(vel_v, n_unit) > Scalar(0))
                    n_unit = -n_unit;

                found = true;
                best_t = t_hit;
                best_n = n_unit;
                best_pos = pos_v + vel_v * t_hit;
                }

            if (!found)
                // no collision
                {
                pos_v += dt_remain * vel_v;
                dt_remain = Scalar(0);
                collide = false;
                }
            else
                {
                // backtrack the particle for dt to get to point of contact
                pos_v = best_pos;

                // apply boundary condition
                if (m_geom->getNoSlip())
                    {
                    vel_v = -vel_v;
                    } 
                else
                    {
                    vel_v = vel_v - Scalar(2.0) * dot(vel_v, best_n) * best_n;
                    }
                
                dt_remain -= best_t;
                collide = (dt_remain > Scalar(0));
                }
            } while (dt_remain > 0 && collide);

        // finalize velocity update
        pos = vec_to_scalar3(pos_v);
        vel = vec_to_scalar3(vel_v);
        vel += Scalar(0.5) * m_mpcd_dt * force.evaluate(pos) / mass;

        // wrap and update the position
        int3 image = make_int3(0, 0, 0);
        box.wrap(pos, image);

        h_pos.data[cur_p] = make_scalar4(pos.x, pos.y, pos.z, __int_as_scalar(type));
        h_vel.data[cur_p]
            = make_scalar4(vel.x, vel.y, vel.z, __int_as_scalar(mpcd::detail::NO_CELL));
        }

    // particles have moved, so the cell cache is no longer valid
    m_mpcd_pdata->invalidateCellCache();
    }
        
namespace detail
    {
//! Export mpcd::StreamingMethod to python
/*!
 * \param m Python module to export to
 */
template<class Force> void export_TriangulatedGeometryStreamingMethod(pybind11::module& m)
    {
    const std::string name = "TriangulatedGeometryStreamingMethod" + Force::getName();
    pybind11::class_<mpcd::TriangulatedGeometryStreamingMethod<Force>,
                     mpcd::StreamingMethod,
                     std::shared_ptr<mpcd::TriangulatedGeometryStreamingMethod<Force>>>(
        m,
        name.c_str())
        .def(pybind11::init<std::shared_ptr<SystemDefinition>,
                            unsigned int,
                            unsigned int,
                            int,
                            std::shared_ptr<TriangulatedGeometry>,
                            std::shared_ptr<Force>>())
        .def_property_readonly("geometry",
                               &mpcd::TriangulatedGeometryStreamingMethod<Force>::getGeometry)
        .def_property_readonly("mpcd_particle_force",
                               &mpcd::TriangulatedGeometryStreamingMethod<Force>::getForce);
    }
    } // end namespace detail
    } // end namespace mpcd
    } // end namespace hoomd
#endif // MPCD_TRIANGULATED_GEOMETRY_STREAMING_METHOD_H_