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
//! Ray-triangle intersection
/*!
 * This method computes the intersection between a particle trajectory and a triangle. The particle
 * trajectory is treated as a ray starting at \a pos with direction \a vel over the remaining
 * timestep \a dt_remain. If an intersection occurs within this interval, the hit time is returned
 * in
 * \a t_hit. This implementation follows the watertight ray-triangle intersection algorithm
 * described in: Woop, S., Benthin, C., Wald, I. (2013). Watertight Ray/Triangle Intersection.
 *
 * \param pos_v Particle position
 * \param vel_v Particle direction
 * \param a, b, c Triangle vertices
 * \param dt_remain Remaining timestep
 * \param eps Numerical tolerance
 * \param t_hit Time of intersection
 *
 * This returns True if intersection occurs within dt_remain, false otherwise.
 *
 */
inline bool intersectTriangle(const Scalar3& pos_v,
                              const Scalar3& vel_v,
                              const Scalar3& a,
                              const Scalar3& b,
                              const Scalar3& c,
                              const Scalar dt_remain,
                              const Scalar eps,
                              Scalar& t_hit)

    {
    // calculate dimension where the ray direction is maximal
    const Scalar ax = std::fabs(vel_v.x);
    const Scalar ay = std::fabs(vel_v.y);
    const Scalar az = std::fabs(vel_v.z);

    int kz = 0;
    if (ay >= ax && ay >= az)
        kz = 1;
    else if (az >= ax && az >= ay)
        kz = 2;
    int kx = kz + 1;
    if (kx == 3)
        kx = 0;
    int ky = kx + 1;
    if (ky == 3)
        ky = 0;

    const Scalar vel_comp[3] = {vel_v.x, vel_v.y, vel_v.z};
    const Scalar det = vel_comp[kz];
    if (std::fabs(det) <= eps)
        return false;

    if (det < Scalar(0.0))
        std::swap(kx, ky);

    // calculate shear constants
    const Scalar Sx = vel_comp[kx] / det;
    const Scalar Sy = vel_comp[ky] / det;
    const Scalar Sz = Scalar(1.0) / det;

    // calculate vertices relative to ray origin
    const Scalar3 A = a - pos_v;
    const Scalar3 B = b - pos_v;
    const Scalar3 C = c - pos_v;

    const Scalar A_comp[3] = {A.x, A.y, A.z};
    const Scalar B_comp[3] = {B.x, B.y, B.z};
    const Scalar C_comp[3] = {C.x, C.y, C.z};

    // apply shear and scale
    const Scalar Ax = A_comp[kx] - Sx * A_comp[kz];
    const Scalar Ay = A_comp[ky] - Sy * A_comp[kz];
    const Scalar Bx = B_comp[kx] - Sx * B_comp[kz];
    const Scalar By = B_comp[ky] - Sy * B_comp[kz];
    const Scalar Cx = C_comp[kx] - Sx * C_comp[kz];
    const Scalar Cy = C_comp[ky] - Sy * C_comp[kz];

    // calculate scaled barycentric coordinates
    Scalar u = Cx * By - Cy * Bx;
    Scalar v = Ax * Cy - Ay * Cx;
    Scalar w = Bx * Ay - By * Ax;

    if (u == Scalar(0.0) || v == Scalar(0.0) || w == Scalar(0.0))
        {
        const double CxBy = (double)Cx * (double)By;
        const double CyBx = (double)Cy * (double)Bx;
        u = static_cast<Scalar>(CxBy - CyBx);

        const double AxCy = (double)Ax * (double)Cy;
        const double AyCx = (double)Ay * (double)Cx;
        v = static_cast<Scalar>(AxCy - AyCx);

        const double BxAy = (double)Bx * (double)Ay;
        const double ByAx = (double)By * (double)Ax;
        w = static_cast<Scalar>(BxAy - ByAx);
        }

    if ((u < Scalar(0.0) || v < Scalar(0.0) || w < Scalar(0.0))
        && (u > Scalar(0.0) || v > Scalar(0.0) || w > Scalar(0.0)))
        return false;

    const Scalar inv_det = Scalar(1.0) / (u + v + w);
    if (!std::isfinite((double)inv_det))
        return false;

    // scaled z
    const Scalar Az = Sz * A_comp[kz];
    const Scalar Bz = Sz * B_comp[kz];
    const Scalar Cz = Sz * C_comp[kz];

    const Scalar t = (u * Az + v * Bz + w * Cz) * inv_det;

    if (t <= Scalar(0.0) || t > dt_remain)
        return false;

    t_hit = t;
    return true;
    }

//! Build an AABB enclosing a particle path segment
/*!
 * This builds the axis-aligned bounding box enclosing the particle path from
 * \a pos to \a pos + vel * dt.
 *
 * \param pos Particle position
 * \param vel Particle direction
 * \param dt Remaining timestep
 */
inline hoomd::detail::AABB makePathAABB(const Scalar3& pos, const Scalar3& vel, const Scalar dt)
    {
    const Scalar3 pos_end = pos + dt * vel;
    const Scalar3 lower = make_scalar3(std::min(pos.x, pos_end.x),
                                       std::min(pos.y, pos_end.y),
                                       std::min(pos.z, pos_end.z));
    const Scalar3 upper = make_scalar3(std::max(pos.x, pos_end.x),
                                       std::max(pos.y, pos_end.y),
                                       std::max(pos.z, pos_end.z));

    return hoomd::detail::AABB(vec3<Scalar>(lower.x, lower.y, lower.z),
                               vec3<Scalar>(upper.x, upper.y, upper.z));
    }

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
    std::shared_ptr<Force> m_force;               //!< Solvent force
    };

/*!
 * \param timestep Current time to stream
 */
template<class Force> void TriangulatedGeometryStreamingMethod<Force>::stream(uint64_t timestep)
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
    ArrayHandle<Scalar3> h_vertices(m_geom->getVertices(),
                                    access_location::host,
                                    access_mode::read);
    ArrayHandle<uint3> h_triangles(m_geom->getTriangles(),
                                   access_location::host,
                                   access_mode::read);
    const Scalar mass = m_mpcd_pdata->getMass();

    // default construct a force if one is not set
    const Force force = (m_force) ? *m_force : Force();

    // triangle candidates through candidates
    std::vector<unsigned int> candidate_triangles;

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

        Scalar3 pos_v(pos);
        Scalar3 vel_v(vel);
        do
            {
            bool found = false;

            // to keep the earliest hit
            Scalar best_t = dt_remain;
            unsigned int best_tri = 0;
            const Scalar eps = Scalar(1e-12);

            // broad search: query candidate triangles overlapping the particle path AABB
            candidate_triangles.clear();
            const hoomd::detail::AABB path_aabb = makePathAABB(pos_v, vel_v, dt_remain);
            m_geom->getTriangleTree().query(candidate_triangles, path_aabb);

            for (const auto cur_tri : candidate_triangles)
                {
                const uint3 triangles = h_triangles.data[cur_tri];

                const Scalar3 a(h_vertices.data[triangles.x]);
                const Scalar3 b(h_vertices.data[triangles.y]);
                const Scalar3 c(h_vertices.data[triangles.z]);

                // calculate normal vector of triangle
                const Scalar3 e1 = b - a;
                const Scalar3 e2 = c - a;
                const Scalar3 n = cross(e1, e2);

                // exclude particles moving away from the triangle
                if (dot(vel_v, n) <= Scalar(0))
                    continue;

                // narrow search: exact ray-triangle intersection
                Scalar t_hit;
                if (!intersectTriangle(pos_v, vel_v, a, b, c, dt_remain, eps, t_hit))
                    continue;

                if (t_hit >= best_t)
                    continue;

                found = true;
                best_t = t_hit;
                best_tri = cur_tri;
                }

            if (found)
                {
                // retrieve the triangle that produces earliest hit
                const uint3 triangle = h_triangles.data[best_tri];

                const Scalar3 a(h_vertices.data[triangle.x]);
                const Scalar3 b(h_vertices.data[triangle.y]);
                const Scalar3 c(h_vertices.data[triangle.z]);

                // compute triangle normal
                const Scalar3 e1 = b - a;
                const Scalar3 e2 = c - a;
                const Scalar3 n = cross(e1, e2);
                const Scalar nn = dot(n, n);

                Scalar3 n_unit = n * (Scalar(1.0) / fast::sqrt(nn));

                // backtrack the particle for dt to get to point of contact
                pos_v += vel_v * best_t;

                // apply boundary condition
                if (m_geom->getNoSlip())
                    {
                    vel_v = -vel_v;
                    }
                else
                    {
                    vel_v = vel_v - Scalar(2.0) * dot(vel_v, n_unit) * n_unit;
                    }

                dt_remain -= best_t;
                collide = (dt_remain > Scalar(0));
                }
            else
                {
                pos_v += dt_remain * vel_v;
                dt_remain = Scalar(0);
                collide = false;
                }
            } while (dt_remain > 0 && collide);

        // finalize velocity update
        pos = pos_v;
        vel = vel_v;
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
