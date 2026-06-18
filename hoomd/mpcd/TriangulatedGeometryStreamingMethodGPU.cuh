// Copyright (c) 2009-2025 The Regents of the University of Michigan.
// Part of HOOMD-blue, released under the BSD 3-Clause License.

#ifndef MPCD_TRIANGULATED_GEOMETRY_STREAMING_METHOD_GPU_CUH_
#define MPCD_TRIANGULATED_GEOMETRY_STREAMING_METHOD_GPU_CUH_

/*!
 * \file mpcd/TriangulatedGeometryStreamingMethodGPU.cuh
 * \brief Declaration of CUDA kernels for mpcd::TriangulatedGeometryStreamingMethodGPU
 */

#include <vector>

#include <hip/hip_runtime.h>

#include "ParticleDataUtilities.h"
#include "hoomd/BoxDim.h"
#include "hoomd/HOOMDMath.h"

// forward declaration
namespace neighbor
    {
class LBVH;
class LBVHTraverser;
    } // namespace neighbor

namespace hoomd
    {
namespace mpcd
    {
namespace gpu
    {
//! Common arguments passed to all streaming kernels
struct triangulated_stream_args_t
    {
    //! Constructor
    triangulated_stream_args_t(Scalar4* _d_pos,
                               Scalar4* _d_vel,
                               const Scalar _mass,
                               const BoxDim& _box,
                               const Scalar _dt,
                               const unsigned int _N)
        : d_pos(_d_pos), d_vel(_d_vel), mass(_mass), box(_box), dt(_dt), N(_N)
        {
        }

    Scalar4* d_pos;       //!< Particle positions
    Scalar4* d_vel;       //!< Particle velocities
    const Scalar mass;    //!< Particle mass
    const BoxDim box;     //!< Simulation box
    const Scalar dt;      //!< Timestep
    const unsigned int N; //!< Number of particles
    };

//! Wrapper around neighbor::LBVH for triangulated geometry
/*!
 * This wrapper only exposes data types that are natively supported in HOOMD
 * so that all neighbor-specific templates and structs can be handled only
 * in CUDA code.
 */
class TriangleLBVHWrapper
    {
    public:
    //! Constructor
    TriangleLBVHWrapper();

    //! Destructor
    ~TriangleLBVHWrapper();

    //! Setup the LBVH
    void setup(const Scalar3* vertices,
               const uint3* triangles,
               unsigned int n_triangles,
               hipStream_t stream);

    //! Build the LBVH
    void build(const Scalar3* vertices,
               const uint3* triangles,
               unsigned int n_triangles,
               const Scalar3& lo,
               const Scalar3& hi,
               hipStream_t stream,
               unsigned int block_size);

    //! Get underlying LBVH
    neighbor::LBVH* get()
        {
        return lbvh_;
        }

    //! Get number of primitives in the LBVH
    unsigned int getN() const;

    //! Get the list of tunable parameters
    std::vector<unsigned int> getTunableParameters() const;

    private:
    neighbor::LBVH* lbvh_; //!< Underlying neighbor::LBVH
    };

//! Wrapper around neighbor::LBVHTraverser for triangulated geometry
/*!
 * This wrapper only exposes data types that are natively supported in HOOMD
 * so that all neighbor-specific templates and structs can be handled only
 * in CUDA code.
 */
class TriangleLBVHTraverserWrapper
    {
    public:
    //! Constructor
    TriangleLBVHTraverserWrapper();

    /// Destructor
    ~TriangleLBVHTraverserWrapper();

    //! Setup the LBVH traverser
    void setup(neighbor::LBVH& lbvh, hipStream_t stream);

    //! Traverse the LBVH
    void traverse(const triangulated_stream_args_t& args,
                  const Scalar3* d_vertices,
                  const uint3* d_triangles,
                  Scalar* d_dt_remain,
                  const unsigned int* d_active_idx,
                  unsigned int* d_flags,
                  unsigned int num_active,
                  bool no_slip,
                  neighbor::LBVH& lbvh,
                  hipStream_t stream,
                  unsigned int block_size);

    //! Get the list of tunable parameters
    std::vector<unsigned int> getTunableParameters() const;

    private:
    neighbor::LBVHTraverser* trav_; //!< Underlying neighbor::LBVHTraverser
    };

//! Kernel driver to apply first half-kick + initialize dt_remain and active idx
template<class Force>
hipError_t triangulated_stream_init(const triangulated_stream_args_t& args,
                                    const Force& solvent_force,
                                    Scalar* d_dt_remain,
                                    unsigned int* d_active_idx,
                                    unsigned int block_size);

//! Compact active list using CUB DevicePartition::Flagged
/*!
 * \param cub_tmp             Temporary storage
 * \param cub_tmp_bytes       Number of bytes in temporary storage
 * \param d_active_idx        Current active particle indices
 * \param d_next_active_idx   Compacted output active indices
 * \param d_flags             Flag set to 1 if the particle's dt_remain > 0
 * \param d_num_active        Output count of active particles
 * \param num_active          Current number of active particles
 */
hipError_t triangulated_compact_active(void* cub_tmp,
                                       size_t& cub_tmp_bytes,
                                       unsigned int* d_active_idx,
                                       unsigned int* d_next_active_idx,
                                       unsigned int* d_flags,
                                       unsigned int* d_num_active,
                                       unsigned int num_active);

//! Kernel driver to apply final half-kick + box wrap
template<class Force>
hipError_t triangulated_stream_finalize(const triangulated_stream_args_t& args,
                                        const Force& solvent_force,
                                        unsigned int block_size);

#ifdef __HIPCC__
namespace kernel
    {
template<class Force>
__global__ void triangulated_stream_init_kernel(Scalar4* d_pos,
                                                Scalar4* d_vel,
                                                const Scalar mass,
                                                const Scalar dt,
                                                const unsigned int N,
                                                const Force force,
                                                Scalar* d_dt_remain,
                                                unsigned int* d_active_idx)
    {
    const unsigned int idx = blockDim.x * blockIdx.x + threadIdx.x;
    if (idx >= N)
        return;

    const Scalar4 postype = d_pos[idx];
    const Scalar3 pos = make_scalar3(postype.x, postype.y, postype.z);

    const Scalar4 vel_cell = d_vel[idx];
    Scalar3 vel = make_scalar3(vel_cell.x, vel_cell.y, vel_cell.z);

    vel += Scalar(0.5) * dt * force.evaluate(pos) / mass;

    d_vel[idx] = make_scalar4(vel.x, vel.y, vel.z, vel_cell.w);
    d_dt_remain[idx] = dt;
    d_active_idx[idx] = idx;
    }

template<class Force>
__global__ void triangulated_stream_finalize_kernel(Scalar4* d_pos,
                                                    Scalar4* d_vel,
                                                    const Scalar mass,
                                                    const BoxDim box,
                                                    const Scalar dt,
                                                    const unsigned int N,
                                                    const Force force)
    {
    const unsigned int idx = blockDim.x * blockIdx.x + threadIdx.x;
    if (idx >= N)
        return;

    const Scalar4 postype = d_pos[idx];
    Scalar3 pos = make_scalar3(postype.x, postype.y, postype.z);

    const Scalar4 vel_cell = d_vel[idx];
    Scalar3 vel = make_scalar3(vel_cell.x, vel_cell.y, vel_cell.z);

    vel += Scalar(0.5) * dt * force.evaluate(pos) / mass;

    int3 image = make_int3(0, 0, 0);
    box.wrap(pos, image);

    d_pos[idx] = make_scalar4(pos.x, pos.y, pos.z, postype.w);
    d_vel[idx] = make_scalar4(vel.x, vel.y, vel.z, __int_as_scalar(mpcd::detail::NO_CELL));
    }

    } // end namespace kernel

template<class Force>
hipError_t triangulated_stream_init(const triangulated_stream_args_t& args,
                                    const Force& solvent_force,
                                    Scalar* d_dt_remain,
                                    unsigned int* d_active_idx,
                                    unsigned int block_size)
    {
    unsigned int max_block_size;
    hipFuncAttributes attr;
    hipFuncGetAttributes(
        &attr,
        reinterpret_cast<const void*>(kernel::triangulated_stream_init_kernel<Force>));
    max_block_size = attr.maxThreadsPerBlock;

    const unsigned int run_block_size = min(block_size, max_block_size);
    hipLaunchKernelGGL(kernel::triangulated_stream_init_kernel<Force>,
                       dim3(args.N / run_block_size + 1),
                       dim3(run_block_size),
                       0,
                       0,
                       args.d_pos,
                       args.d_vel,
                       args.mass,
                       args.dt,
                       args.N,
                       solvent_force,
                       d_dt_remain,
                       d_active_idx);
    return hipSuccess;
    }

template<class Force>
hipError_t triangulated_stream_finalize(const triangulated_stream_args_t& args,
                                        const Force& solvent_force,
                                        unsigned int block_size)
    {
    unsigned int max_block_size;
    hipFuncAttributes attr;
    hipFuncGetAttributes(
        &attr,
        reinterpret_cast<const void*>(kernel::triangulated_stream_finalize_kernel<Force>));
    max_block_size = attr.maxThreadsPerBlock;

    const unsigned int run_block_size = min(block_size, max_block_size);
    hipLaunchKernelGGL(kernel::triangulated_stream_finalize_kernel<Force>,
                       dim3(args.N / run_block_size + 1),
                       dim3(run_block_size),
                       0,
                       0,
                       args.d_pos,
                       args.d_vel,
                       args.mass,
                       args.box,
                       args.dt,
                       args.N,
                       solvent_force);
    return hipSuccess;
    }
#endif // __HIPCC__

    } // end namespace gpu
    } // end namespace mpcd
    } // end namespace hoomd
#endif // MPCD_TRIANGULATED_GEOMETRY_STREAMING_METHOD_GPU_CUH_
