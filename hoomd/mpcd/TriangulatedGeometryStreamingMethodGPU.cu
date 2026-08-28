// Copyright (c) 2009-2026 The Regents of the University of Michigan.
// Part of HOOMD-blue, released under the BSD 3-Clause License.

/*!
 * \file mpcd/TriangulatedGeometryStreamingMethodGPU.cu
 * \brief CUDA kernel implementations for mpcd::TriangulatedGeometryStreamingMethodGPU
 */

#include "ParticleDataUtilities.h"
#include "TriangulatedGeometryStreamingMethodGPU.cuh"

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wconversion"
#include <hipcub/hipcub.hpp>
#pragma GCC diagnostic pop

#include <neighbor/neighbor.h>

namespace hoomd
    {
namespace mpcd
    {
namespace gpu
    {
#define DEVICE __device__ __forceinline__

//! Insert operation for triangulated geometry
/*!
 * Provides an axis-aligned bounding box per triangle to neighbor::LBVH during
 * construction. The bounding box of each triangle is the AABB of its three vertices.
 * One primitive in the LBVH corresponds to one triangle, so the LBVH primitive index
 * is the triangle index.
 */
struct TriangleInsertOp
    {
    //! Constructor
    /*!
     * \param vertices_ List of vertices.
     * \param triangles_ List of triangles.
     * \param N_ Number of primitives (triangles) to insert.
     */
    TriangleInsertOp(const ShortReal3* vertices_, const uint3* triangles_, unsigned int N_)
        : vertices(vertices_), triangles(triangles_), N(N_)
        {
        }

    //! Construct bounding box for a triangle
    /*!
     * \param idx Nominal index of the primitive [0,N).
     * \returns A neighbor::BoundingBox covering the three vertices.
     */
    DEVICE neighbor::BoundingBox get(unsigned int idx) const
        {
        const uint3 tri = triangles[idx];
        const ShortReal3 a = vertices[tri.x];
        const ShortReal3 b = vertices[tri.y];
        const ShortReal3 c = vertices[tri.z];

        const float3 lo = make_float3(min(a.x, min(b.x, c.x)),
                                      min(a.y, min(b.y, c.y)),
                                      min(a.z, min(b.z, c.z)));
        const float3 hi = make_float3(max(a.x, max(b.x, c.x)),
                                      max(a.y, max(b.y, c.y)),
                                      max(a.z, max(b.z, c.z)));

        return neighbor::BoundingBox(lo, hi);
        }

    //! Get the number of primitives
    __host__ DEVICE unsigned int size() const
        {
        return N;
        }

    const ShortReal3* vertices; //!< Triangle vertices
    const uint3* triangles;     //!< Triangle vertex indices
    const unsigned int N;       //!< Number of primitives
    };

//! Path query operation for active particles
/*!
 * For each currently active particle (one whose \a dt_remain > 0), the query
 * bounding volume is the AABB swept by the segment from \a pos to
 * \a pos + \a dt_remain * \a vel. Any LBVH leaf whose box overlaps this AABB
 * is forwarded to CollisionOutputOp::process as a candidate triangle for the
 * ray-triangle intersection test.
 */
struct PathQueryOp
    {
    //! Constructor
    /*!
     * \param d_pos_ Particle positions.
     * \param d_vel_ Particle velocities.
     * \param d_dt_remain_ Remaining timestep.
     * \param d_active_idx_ Active particle indices.
     * \param N_ Number of active particles
     */
    PathQueryOp(const Scalar4* d_pos_,
                const Scalar4* d_vel_,
                const Scalar* d_dt_remain_,
                const unsigned int* d_active_idx_,
                unsigned int N_)
        : d_pos(d_pos_), d_vel(d_vel_), d_dt_remain(d_dt_remain_), d_active_idx(d_active_idx_),
          N(N_)
        {
        }

    //! Data stored per thread for traversal
    /*!
     * Holds the particle's position, velocity, remaining timestep, and true index.
     */
    struct ThreadData
        {
        DEVICE ThreadData(Scalar3 pos_, Scalar3 vel_, Scalar dt_remain_, unsigned int idx_)
            : pos(pos_), vel(vel_), dt_remain(dt_remain_), idx(idx_)
            {
            }

        Scalar3 pos;      //!< Particle position.
        Scalar3 vel;      //!< Particle velocity.
        Scalar dt_remain; //!< Remaining timestep.
        unsigned int idx; //!< True particle index.
        };

    // Traversal volume type used by neighbor::LBVHTraverser
    typedef neighbor::BoundingBox Volume;

    //! Loads the per-thread data
    /*!
     * \param i Traversal thread index in [0, num_active).
     * \returns The ThreadData required for traversal.
     *
     * The threadData is only loaded for active particles. The true particle
     * index is mapped through \a d_active_idx together with its position, velocity,
     * and remaining timestep.
     */
    DEVICE ThreadData setup(unsigned int i) const
        {
        const unsigned int pidx = d_active_idx[i];
        const Scalar4 postype = d_pos[pidx];
        const Scalar4 vel_cell = d_vel[pidx];
        return ThreadData(make_scalar3(postype.x, postype.y, postype.z),
                          make_scalar3(vel_cell.x, vel_cell.y, vel_cell.z),
                          d_dt_remain[pidx],
                          pidx);
        }

    //! Return the traversal volume subject to a translation
    /*!
     * \param q The current thread data.
     * \param image The image vector for traversal.
     * \returns The traversal bounding volume.
     *
     * The ThreadData is converted to a search volume. The search volume is the
     * AABB of the path from \a q.pos to \a q.pos + \a q.dt_remain * \a q.vel.
     * The \a image argumet is ignored as triangles are already unwrapped.
     */
    DEVICE Volume get(const ThreadData& q, const float3&) const
        {
        const Scalar3 end = q.pos + q.dt_remain * q.vel;
        const Scalar3 lo
            = make_scalar3(min(q.pos.x, end.x), min(q.pos.y, end.y), min(q.pos.z, end.z));
        const Scalar3 hi
            = make_scalar3(max(q.pos.x, end.x), max(q.pos.y, end.y), max(q.pos.z, end.z));
        return neighbor::BoundingBox(lo, hi);
        }

    //! Perform the overlap test with the LBVH
    /*!
     * \param v Traversal volume.
     * \param box Box in LBVH to intersect with.
     * \returns True if the volume and box overlap.
     *
     * This performs broad-phase search through the overlap test between particle's
     * path AABB and triangle AABB.
     */
    DEVICE bool overlap(const Volume& v, const neighbor::BoundingBox& box) const
        {
        return v.overlap(box);
        }

    //! Refine the rough overlap test with a primitive
    /*!
     * \param q The current thread data.
     * \param primitive Index of the intersected primitive.
     * \returns True if the volumes still overlap after refinement.
     *
     * Every candidate is forwarded to CollisionOutputOp::process to perform the
     * exact ray-triangle intersection test, so this always returns \a True.
     */
    DEVICE bool refine(const ThreadData&, const int) const
        {
        return true;
        }

    //! Get the number of query threads
    __host__ DEVICE unsigned int size() const
        {
        return N;
        }

    const Scalar4* d_pos;             //!< Particle positions
    const Scalar4* d_vel;             //!< Particle velocities
    const Scalar* d_dt_remain;        //!< Remaining timesteps
    const unsigned int* d_active_idx; //!< Active particle indices
    const unsigned int N;             //!< Number of active particles
    };

//! Collision output operation for active particles
/*!
 * For each candidate triangle visited during traversal, process() performs
 * a Woop-style watertight ray-triangle intersection between the particle's
 * remaining path and the triangle, keeping track of the earliest valid hit.
 * Once traversal is done, finalize() advances the particles, reflects its
 * velocity based on the boundary condition.
 */
struct CollisionOutputOp
    {
    //! Constructor
    /*!
     * \param vertices_ Triangle vertices
     * \param triangles_ Triangle vertex indices
     * \param d_pos_ Particle positions
     * \param d_vel_ Particle velocities
     * \param d_dt_remain_ Remaining timesteps
     * \param d_flags_ flag set to 1 if the particle's dt_remain > 0
     * \param no_slip Boundary condition at the wall
     * \param eps_ Numerical tolerance
     */
    CollisionOutputOp(const ShortReal3* vertices_,
                      const uint3* triangles_,
                      Scalar4* d_pos_,
                      Scalar4* d_vel_,
                      Scalar* d_dt_remain_,
                      unsigned int* d_flags_,
                      bool no_slip_,
                      Scalar eps_)
        : vertices(vertices_), triangles(triangles_), d_pos(d_pos_), d_vel(d_vel_),
          d_dt_remain(d_dt_remain_), d_flags(d_flags_), no_slip(no_slip_), eps(eps_)
        {
        }

    //! Thread-local data
    /*!
     * Stores the particle state during traversal and accumulates the earliest valid
     * triangle intersection found for the current particle. Two indices are kept:
     *  1) \a i is the traversal thread index, used to write \a d_flags, which runs in
     *      parallel to the active list.
     *  2) \a idx is the true particle index, used to write \a d_pos, \a d_vel, and
     *      \a d_dt_remain.
     */
    struct ThreadData
        {
        //! Constructor
        /*!
         * \param i_ Thread index
         * \param idx_ True particle index
         * \param pos_ Particle positions
         * \param vel_ Particle velocities
         * \param dt_remain_ Remaining timesteps
         */
        DEVICE ThreadData(unsigned int i_,
                          unsigned int idx_,
                          Scalar3 pos_,
                          Scalar3 vel_,
                          Scalar dt_remain_)
            : i(i_), idx(idx_), pos(pos_), vel(vel_), dt_remain(dt_remain_), best_t(dt_remain_),
              best_tri(0xffffffff), found(false)
            {
            }

        unsigned int i;        //!< Thread index
        unsigned int idx;      //!< True particle index
        Scalar3 pos;           //!< Particle position
        Scalar3 vel;           //!< Particle velocity
        Scalar dt_remain;      //!< Remaining timestep
        Scalar best_t;         //!< Time of the earliest hit
        unsigned int best_tri; //!< Triangle index of the earliest hit
        bool found;            //!< Whether any valid hit was found
        };

    //! Setup the thread data
    /*!
     * \param i Thread index
     * \param q Thread-local query data.
     * \returns The ThreadData for output.
     *
     * \tparam Type of QueryData.
     *
     * This setup function can poach data from the query data in order to save loads.
     * The thread index \a i is kept so finalize() can write d_flags[i]
     */
    template<class QueryDataT> DEVICE ThreadData setup(unsigned int i, const QueryDataT& q) const
        {
        return ThreadData(i, q.idx, q.pos, q.vel, q.dt_remain);
        }

    //! Processes a newly intersected primitive.
    /*!
     * \param t The thread's output data.
     * \param primitive The index of the candidate triangle to test.
     *
     * Performs the watertight ray-triangle intersection algorithm described in:
     * Woop, S., Benthin, C., Wald, I. (2013). Watertight Ray/Triangle Intersection.
     */
    DEVICE void process(ThreadData& t, const int primitive) const
        {
        const uint3 tri = triangles[primitive];
        const ShortReal3 a = vertices[tri.x];
        const ShortReal3 b = vertices[tri.y];
        const ShortReal3 c = vertices[tri.z];

        const Scalar3 aa = make_scalar3(a.x, a.y, a.z);
        const Scalar3 bb = make_scalar3(b.x, b.y, b.z);
        const Scalar3 cc = make_scalar3(c.x, c.y, c.z);

        const Scalar3 e1 = bb - aa;
        const Scalar3 e2 = cc - aa;
        const Scalar3 n = cross(e1, e2);

        // exclude particles moving away from the triangle
        if (dot(t.vel, n) <= Scalar(0))
            return;

        // calculate dimension where the ray direction is maximal
        const Scalar ax = fabs(t.vel.x);
        const Scalar ay = fabs(t.vel.y);
        const Scalar az = fabs(t.vel.z);

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

        // vel components by axis
        const Scalar vel_comp[3] = {t.vel.x, t.vel.y, t.vel.z};
        const Scalar det = vel_comp[kz];
        if (fabs(det) <= eps)
            return;

        if (det < Scalar(0.0))
            {
            int tmp = kx;
            kx = ky;
            ky = tmp;
            }

        // calculate shear constants
        const Scalar Sx = vel_comp[kx] / det;
        const Scalar Sy = vel_comp[ky] / det;
        const Scalar Sz = Scalar(1.0) / det;

        // calculate vertices relative to ray origin
        const Scalar3 A = aa - t.pos;
        const Scalar3 B = bb - t.pos;
        const Scalar3 C = cc - t.pos;

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
            return;

        const Scalar inv_det = Scalar(1.0) / (u + v + w);
        if (!isfinite((double)inv_det))
            return;

        // scaled z
        const Scalar Az = Sz * A_comp[kz];
        const Scalar Bz = Sz * B_comp[kz];
        const Scalar Cz = Sz * C_comp[kz];

        const Scalar t_hit = (u * Az + v * Bz + w * Cz) * inv_det;

        if (t_hit <= Scalar(0.0) || t_hit > t.dt_remain)
            return;

        if (t_hit < t.best_t)
            {
            t.best_t = t_hit;
            t.best_tri = static_cast<unsigned int>(primitive);
            t.found = true;
            }
        }

    //! Finish the output job once the thread is ready to terminate.
    /*!
     * \param t The thread's output data
     *
     * If a hit was found, the particle is advanced to the contact point, its
     * remaining timestep is reduced by the time consumed, and its velocity is
     * reflected based on the boundary condition (slip or no-slip). Otherwise, the
     * particle advances by its full remaining timestep and dt_remain becomes zero.
     */
    DEVICE void finalize(ThreadData& t) const
        {
        if (t.found)
            {
            // retrieve the triangle that produces earliest hit
            const uint3 tri = triangles[t.best_tri];
            const ShortReal3 a = vertices[tri.x];
            const ShortReal3 b = vertices[tri.y];
            const ShortReal3 c = vertices[tri.z];

            const Scalar3 aa = make_scalar3(a.x, a.y, a.z);
            const Scalar3 bb = make_scalar3(b.x, b.y, b.z);
            const Scalar3 cc = make_scalar3(c.x, c.y, c.z);

            // compute triangle normal
            const Scalar3 e1 = bb - aa;
            const Scalar3 e2 = cc - aa;
            const Scalar3 n = cross(e1, e2);
            const Scalar3 n_unit = n * (Scalar(1) / fast::sqrt(dot(n, n)));

            // backtrack the particle for dt to get to point of contact
            t.pos += t.vel * t.best_t;
            t.dt_remain -= t.best_t;

            if (no_slip)
                t.vel = -t.vel;
            else
                t.vel = t.vel - Scalar(2) * dot(t.vel, n_unit) * n_unit;
            }
        else
            {
            t.pos += t.dt_remain * t.vel;
            t.dt_remain = Scalar(0);
            }

        const Scalar4 postype = d_pos[t.idx];
        const Scalar4 vel_cell = d_vel[t.idx];
        d_pos[t.idx] = make_scalar4(t.pos.x, t.pos.y, t.pos.z, postype.w);
        d_vel[t.idx] = make_scalar4(t.vel.x, t.vel.y, t.vel.z, vel_cell.w);
        d_dt_remain[t.idx] = t.dt_remain;
        d_flags[t.i] = (t.dt_remain > Scalar(0)) ? 1u : 0u;
        }

    const ShortReal3* vertices; //!< Triangle vertices
    const uint3* triangles;     //!< Triangle vertex indices
    Scalar4* d_pos;             //!< Particle positions
    Scalar4* d_vel;             //!< Particle velocities
    Scalar* d_dt_remain;        //!< Remaining timestep
    unsigned int* d_flags;      //!< Flag set to 1 if the particle's dt_remain > 0
    bool no_slip;               //!< Boundary condition at the wall
    Scalar eps;                 //!< Numerical tolerance
    };

//! Driver function implementations
hipError_t triangulated_compact_active(void* cub_tmp,
                                       size_t& cub_tmp_bytes,
                                       unsigned int* d_active_idx,
                                       unsigned int* d_next_active_idx,
                                       unsigned int* d_flags,
                                       unsigned int* d_num_active,
                                       unsigned int num_active)
    {
    hipcub::DevicePartition::Flagged(cub_tmp,
                                     cub_tmp_bytes,
                                     d_active_idx,
                                     d_flags,
                                     d_next_active_idx,
                                     d_num_active,
                                     num_active);
    return hipSuccess;
    }

//! Host function to convert a double to a float in round-down mode
float double2float_rd(double x)
    {
    float xf = static_cast<float>(x);
    if (static_cast<double>(xf) > x)
        {
        xf = std::nextafterf(xf, -std::numeric_limits<float>::infinity());
        }
    return xf;
    }

//! Host function to convert a double to a float in round-up mode
float double2float_ru(double x)
    {
    float xf = static_cast<float>(x);
    if (static_cast<double>(xf) < x)
        {
        xf = std::nextafterf(xf, std::numeric_limits<float>::infinity());
        }
    return xf;
    }

/*!
 * Initializes the shared pointer for the underlying LBVH.
 */
TriangleLBVHWrapper::TriangleLBVHWrapper()
    {
    lbvh_ = new neighbor::LBVH();
    }

TriangleLBVHWrapper::~TriangleLBVHWrapper()
    {
    delete lbvh_;
    }

/*!
 * \param vertices Triangle vertices
 * \param triangles Triangle vertex indices
 * \param num_triangles Number of triangles
 * \param stream CUDA stream for execution
 */
void TriangleLBVHWrapper::setup(const ShortReal3* vertices,
                                const uint3* triangles,
                                unsigned int num_triangles,
                                hipStream_t stream)
    {
    TriangleInsertOp insert(vertices, triangles, num_triangles);
    lbvh_->setup(stream, insert);
    }

/*!
 * \param vertices Triangle vertices
 * \param triangles Triangle vertex indices
 * \param num_triangles Number of triangles
 * \param lo Lower bound of box
 * \param hi Upper bound of box
 * \param stream CUDA stream for execution
 * \param block_size CUDA block size for execution
 *
 * If HOOMD is using double-precision Scalars, then the lo and hi bounds of the
 * box are internally converted to floats using round-down and round-up modes,
 * respectively, which conserves the original box.
 */
void TriangleLBVHWrapper::build(const ShortReal3* vertices,
                                const uint3* triangles,
                                unsigned int num_triangles,
                                const Scalar3& lo,
                                const Scalar3& hi,
                                hipStream_t stream,
                                unsigned int block_size)
    {
#if HOOMD_LONGREAL_SIZE == 64
    float3 lof = make_float3(double2float_rd(lo.x), double2float_rd(lo.y), double2float_rd(lo.z));
    float3 hif = make_float3(double2float_ru(hi.x), double2float_ru(hi.y), double2float_ru(hi.z));
#else
    float3 lof = lo;
    float3 hif = hi;
#endif

    TriangleInsertOp insert(vertices, triangles, num_triangles);
    lbvh_->build(neighbor::LBVH::LaunchParameters(block_size, stream), insert, lof, hif);
    }

unsigned int TriangleLBVHWrapper::getN() const
    {
    return lbvh_->getN();
    }

std::vector<unsigned int> TriangleLBVHWrapper::getTunableParameters() const
    {
    return lbvh_->getTunableParameters();
    }

/*!
 * Initializes the shared pointer for the underlying LBVHTraverser.
 */
TriangleLBVHTraverserWrapper::TriangleLBVHTraverserWrapper()
    {
    trav_ = new neighbor::LBVHTraverser();
    }

TriangleLBVHTraverserWrapper::~TriangleLBVHTraverserWrapper()
    {
    delete trav_;
    }

/*!
 * \param lbvh LBVH to traverse
 * \param stream CUDA stream for execution
 *
 * A NullTransformOp is used as the LBVH primitive order is the triangle order.
 */
void TriangleLBVHTraverserWrapper::setup(neighbor::LBVH& lbvh, hipStream_t stream)
    {
    neighbor::NullTransformOp transform;
    trav_->setup(stream, lbvh, transform);
    }

/*!
 * \param args Common arguments for a streaming kernel
 * \param d_vertices Triangle vertices
 * \param d_triangles Triangle vertex indices
 * \param d_dt_remain Remaining timestep
 * \param d_active_idx Active particle indices
 * \param d_flags Flags set to 1 if the particle's dt_remain > 0
 * \param num_active Number of active particles
 * \param no_slip Boundary condition at the wall
 * \param lbvh LBVH to traverse
 * \param stream CUDA stream for execution
 * \param block_size CUDA block size for execution
 *
 * A single zero image (neighbor::SelfOp) is used as the triangle is already unwrapped
 * given the unwrap distance. A NullTransformOp is used as no primitive mapping is needed.
 */
void TriangleLBVHTraverserWrapper::traverse(const triangulated_stream_args_t& args,
                                            const ShortReal3* d_vertices,
                                            const uint3* d_triangles,
                                            Scalar* d_dt_remain,
                                            const unsigned int* d_active_idx,
                                            unsigned int* d_flags,
                                            unsigned int num_active,
                                            bool no_slip,
                                            neighbor::LBVH& lbvh,
                                            hipStream_t stream,
                                            unsigned int block_size)
    {
    PathQueryOp query(args.d_pos, args.d_vel, d_dt_remain, d_active_idx, num_active);

    CollisionOutputOp output(d_vertices,
                             d_triangles,
                             args.d_pos,
                             args.d_vel,
                             d_dt_remain,
                             d_flags,
                             no_slip,
                             Scalar(1e-12));

    neighbor::NullTransformOp transform;
    neighbor::SelfOp translate;

    trav_->traverse(neighbor::LBVHTraverser::LaunchParameters(block_size, stream),
                    lbvh,
                    query,
                    output,
                    translate,
                    transform);
    }

std::vector<unsigned int> TriangleLBVHTraverserWrapper::getTunableParameters() const
    {
    return trav_->getTunableParameters();
    }

    } // end namespace gpu
    } // end namespace mpcd
    } // end namespace hoomd
