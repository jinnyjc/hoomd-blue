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

                // find two edges
                const vec3<Scalar> e1 = b - a;
                const vec3<Scalar> e2 = c - a;

                // find the determinant
                vec3<Scalar> pvec = cross(vel_v, e2);
                const Scalar det = dot(e1, pvec);

                // ray and triangle are parallel if det is close to 0
                if (fabs(det) < eps)
                    continue;
                
                // inside-outside test
                const Scalar inv_det = Scalar(1.0) / det;
                const vec3<Scalar> tvec = pos_v - a;
                const Scalar u = dot(tvec, pvec) * inv_det;

                if (u < Scalar(0) || u > Scalar(1))
                    continue;

                const vec3<Scalar> qvec = cross(tvec, e1);
                const Scalar v = dot(vel_v, qvec) * inv_det;

                if (v < Scalar(0) || u + v > Scalar(1))
                    continue;

                // find the time the ray hits the triangle
                const Scalar t_hit = dot(e2, qvec) * inv_det;

                // continue if time  < 0 or time > dt_remain
                if (t_hit <= eps || t_hit > dt_remain)
                    continue;

                // choose the earliest hit among all triangles
                if (t_hit >= best_t)
                    continue;

                // compute triangle normal
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