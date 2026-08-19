// Copyright (c) 2009-2026 The Regents of the University of Michigan.
// Part of HOOMD-blue, released under the BSD 3-Clause License.

/*!
 * \file mpcd/RejectionVirtualParticleFillerGPU.cuh
 * \brief Declaration of CUDA kernels for mpcd::RejectionVirtualParticleFillerGPU
 */

#ifndef MPCD_REJECTION_VIRTUAL_PARTICLE_FILLER_GPU_CUH_
#define MPCD_REJECTION_VIRTUAL_PARTICLE_FILLER_GPU_CUH_

#include <cuda_runtime.h>

#include "ParticleDataUtilities.h"
#include "hoomd/BoxDim.h"
#include "hoomd/HOOMDMath.h"

#include "hoomd/RNGIdentifiers.h"
#include "hoomd/RandomNumbers.h"

namespace hoomd
    {
namespace mpcd
    {
namespace gpu
    {
//! Common arguments for classifying which cells need virtual particle filling
struct classify_cells_args_t
    {
    classify_cells_args_t(bool* _d_flags,
                          const BoxDim& _global_box,
                          const BoxDim& _draw_box,
                          const uint3 _local_dim,
                          const uint3 _global_dim,
                          const int3 _origin,
                          const Scalar3 _inv_dim,
                          const unsigned int _num_trials,
                          const unsigned int _seed,
                          const unsigned int _filler_id,
                          const unsigned int _block_size)
        : d_flags(_d_flags), global_box(_global_box), draw_box(_draw_box), local_dim(_local_dim),
          global_dim(_global_dim), origin(_origin), inv_dim(_inv_dim), num_trials(_num_trials),
          seed(_seed), filler_id(_filler_id), block_size(_block_size)
        {
        }

    bool* d_flags;
    const BoxDim global_box;
    const BoxDim draw_box;
    const uint3 local_dim;
    const uint3 global_dim;
    const int3 origin;
    const Scalar3 inv_dim;
    const unsigned int num_trials;
    const unsigned int seed;
    const unsigned int filler_id;
    const unsigned int block_size;
    };

//! Common arguments for drawing virtual particles cell by cell
struct draw_virtual_particles_args_t
    {
    draw_virtual_particles_args_t(Scalar4* _d_tmp_pos,
                                  Scalar4* _d_tmp_vel,
                                  bool* _d_keep_particles,
                                  const unsigned int* _d_fill_cells,
                                  const unsigned int _num_fill_cells,
                                  const BoxDim& _global_box,
                                  const BoxDim& _draw_box,
                                  const uint3 _local_dim,
                                  const uint3 _global_dim,
                                  const int3 _origin,
                                  const Scalar3 _inv_dim,
                                  const Scalar3 _grid_shift,
                                  const Scalar _mean_per_cell,
                                  const unsigned int _max_per_cell,
                                  const Scalar _vel_factor,
                                  const unsigned int _type,
                                  const uint64_t _timestep,
                                  const unsigned int _seed,
                                  const unsigned int _filler_id,
                                  const unsigned int _block_size)
        : d_tmp_pos(_d_tmp_pos), d_tmp_vel(_d_tmp_vel), d_keep_particles(_d_keep_particles),
          d_fill_cells(_d_fill_cells), num_fill_cells(_num_fill_cells), global_box(_global_box),
          draw_box(_draw_box), local_dim(_local_dim), global_dim(_global_dim), origin(_origin),
          inv_dim(_inv_dim), grid_shift(_grid_shift), mean_per_cell(_mean_per_cell),
          max_per_cell(_max_per_cell), vel_factor(_vel_factor), type(_type), timestep(_timestep),
          seed(_seed), filler_id(_filler_id), block_size(_block_size)
        {
        }

    Scalar4* d_tmp_pos;
    Scalar4* d_tmp_vel;
    bool* d_keep_particles;
    const unsigned int* d_fill_cells;
    const unsigned int num_fill_cells;
    const BoxDim global_box;
    const BoxDim draw_box;
    const uint3 local_dim;
    const uint3 global_dim;
    const int3 origin;
    const Scalar3 inv_dim;
    const Scalar3 grid_shift;
    const Scalar mean_per_cell;
    const unsigned int max_per_cell;
    const Scalar vel_factor;
    const unsigned int type;
    const uint64_t timestep;
    const unsigned int seed;
    const unsigned int filler_id;
    const unsigned int block_size;
    };

// Function declarations
template<class Geometry>
cudaError_t __attribute__((visibility("default"))) classify_cells(const classify_cells_args_t& args,
                                                                  const Geometry& geom);

template<class Geometry>
cudaError_t __attribute__((visibility("default")))
draw_virtual_particles(const draw_virtual_particles_args_t& args, const Geometry& geom);

cudaError_t __attribute__((visibility("default")))
compact_virtual_particle_indices(void* d_tmp,
                                 size_t& tmp_bytes,
                                 const bool* d_keep_particles,
                                 const unsigned int num_particles,
                                 unsigned int* d_keep_indices,
                                 unsigned int* d_num_keep);

cudaError_t __attribute__((visibility("default")))
copy_virtual_particles(unsigned int* d_keep_indices,
                       Scalar4* d_pos,
                       Scalar4* d_vel,
                       unsigned int* d_tags,
                       const Scalar4* d_tmp_pos,
                       const Scalar4* d_tmp_vel,
                       const unsigned int first_idx,
                       const unsigned int first_tag,
                       const unsigned int n_virtual,
                       const unsigned int block_size);

#ifdef __HIPCC__
namespace kernel
    {

//! Flatten global cell coordinates into a 1D index for seeding the random number generator
__device__ inline unsigned int globalCellIndex(int gi, int gj, int gk, const uint3& gdim)
    {
    if (gi < 0)
        gi += (int)gdim.x;
    else if (gi >= (int)gdim.x)
        gi -= (int)gdim.x;
    if (gj < 0)
        gj += (int)gdim.y;
    else if (gj >= (int)gdim.y)
        gj -= (int)gdim.y;
    if (gk < 0)
        gk += (int)gdim.z;
    else if (gk >= (int)gdim.z)
        gk -= (int)gdim.z;

    return (unsigned int)gi + gdim.x * ((unsigned int)gj + gdim.y * (unsigned int)gk);
    }

//! Kernel to flag the cells that can hold both fluid and solid
/*!
 * \param d_flags Flags marking the cells that need filling, one per cell
 * \param global_box Global simulation box
 * \param draw_box Box with the extent of the sampling region and the tilt of the global box
 * \param local_dim Cell dimensions owned by this rank
 * \param global_dim Global cell dimensions
 * \param origin Index of the first local cell in the global grid
 * \param inv_dim Fractional width of one cell along each lattice vector
 * \param num_trials Number of trial points drawn per cell
 * \param seed User seed for RNG
 * \param filler_id Identifier for the filler (rng argument)
 * \param geom Confining geometry
 *
 * \tparam Geometry type of the confined geometry \a geom
 *
 * \b implementation
 * We assign one thread per cell. Trial points are drawn in the region the cell can occupy under
 * grid shifting, and the cell is flagged as soon as points are found on both sides of the boundary.
 */
template<class Geometry>
__global__ void classify_cells(bool* d_flags,
                               const BoxDim global_box,
                               const BoxDim draw_box,
                               const uint3 local_dim,
                               const uint3 global_dim,
                               const int3 origin,
                               const Scalar3 inv_dim,
                               const unsigned int num_trials,
                               const unsigned int seed,
                               const unsigned int filler_id,
                               const Geometry geom)
    {
    // one thread per cell
    const unsigned int idx = blockIdx.x * blockDim.x + threadIdx.x;
    const unsigned int num_cells = local_dim.x * local_dim.y * local_dim.z;
    if (idx >= num_cells)
        return;

    const unsigned int i = idx % local_dim.x;
    const unsigned int j = (idx / local_dim.x) % local_dim.y;
    const unsigned int k = idx / (local_dim.x * local_dim.y);

    const int gi = (int)i + origin.x;
    const int gj = (int)j + origin.y;
    const int gk = (int)k + origin.z;

    const Scalar3 f_center = make_scalar3((gi + Scalar(0.5)) * inv_dim.x,
                                          (gj + Scalar(0.5)) * inv_dim.y,
                                          (gk + Scalar(0.5)) * inv_dim.z);
    const Scalar3 cell_center = global_box.makeCoordinates(f_center);

    const Scalar3 half_w = make_scalar3(Scalar(0.5) * draw_box.getL().x,
                                        Scalar(0.5) * draw_box.getL().y,
                                        Scalar(0.5) * draw_box.getL().z);

    hoomd::RandomGenerator rng(
        hoomd::Seed(hoomd::RNGIdentifier::VirtualParticleFiller, 0, seed),
        hoomd::Counter(globalCellIndex(gi, gj, gk, global_dim), filler_id, 1));

    bool found = false;

    // code = -1/+1 after the first usable point; flag the cell when a point lands
    // on the other side
    int code = 0;
    for (unsigned int n = 0; n < num_trials; ++n)
        {
        Scalar3 point = make_scalar3(hoomd::UniformDistribution<Scalar>(-half_w.x, half_w.x)(rng),
                                     hoomd::UniformDistribution<Scalar>(-half_w.y, half_w.y)(rng),
                                     hoomd::UniformDistribution<Scalar>(-half_w.z, half_w.z)(rng));

        int3 img = make_int3(0, 0, 0);
        draw_box.wrap(point, img);
        point += cell_center;

        const Scalar3 f = global_box.makeFraction(point);
        if (f.x < Scalar(0.0) || f.x >= Scalar(1.0) || f.y < Scalar(0.0) || f.y >= Scalar(1.0)
            || f.z < Scalar(0.0) || f.z >= Scalar(1.0))
            continue;

        const bool is_outside = geom.isOutside(point);

        if (code == 0)
            {
            code = is_outside ? 1 : -1;
            continue;
            }
        if ((code == -1 && is_outside) || (code == 1 && !is_outside))
            {
            found = true;
            break;
            }
        }
    d_flags[idx] = found;
    }

//! Kernel to draw a Poisson number of virtual particles in each cell that needs filling
/*!
 * \param d_tmp_pos Temporary positions
 * \param d_tmp_vel Temporary velocities
 * \param d_keep_particles Particle tracking - in/out of given geometry
 * \param d_fill_cells Local 1D indices of the cells that need filling
 * \param num_fill_cells Number of cells that need filling
 * \param global_box Global simulation box
 * \param draw_box Box with the extent of one cell and the tilt of the global box
 * \param local_dim Cell dimensions owned by this rank
 * \param global_dim Global cell dimensions
 * \param origin Index of the first local cell in the global grid
 * \param inv_dim Fractional width of one cell along each lattice vector
 * \param grid_shift Fractional grid shift at this timestep
 * \param mean_per_cell Mean of the Poisson distribution for the number of particles per cell
 * \param max_per_cell Maximum number of particles drawn per cell
 * \param vel_factor Scale factor for uniform normal velocities consistent with particle mass /
 * temperature
 * \param type Particle type for filling
 * \param timestep Current timestep
 * \param seed User seed for RNG
 * \param filler_id Identifier for the filler (rng argument)
 * \param geom Confining geometry
 *
 * \tparam Geometry type of the confined geometry \a geom
 *
 * \b implementation
 * We assign one thread per cell. Each cell owns a fixed block of max_per_cell entries so that
 * the result does not depend on how the threads are scheduled, and the unused entries are
 * flagged so that the compaction drops them.
 */
template<class Geometry>
__global__ void draw_virtual_particles(Scalar4* d_tmp_pos,
                                       Scalar4* d_tmp_vel,
                                       bool* d_keep_particles,
                                       const unsigned int* d_fill_cells,
                                       const unsigned int num_fill_cells,
                                       const BoxDim global_box,
                                       const BoxDim draw_box,
                                       const uint3 local_dim,
                                       const uint3 global_dim,
                                       const int3 origin,
                                       const Scalar3 inv_dim,
                                       const Scalar3 grid_shift,
                                       const Scalar mean_per_cell,
                                       const unsigned int max_per_cell,
                                       const Scalar vel_factor,
                                       const unsigned int type,
                                       const uint64_t timestep,
                                       const unsigned int seed,
                                       const unsigned int filler_id,
                                       const Geometry geom)
    {
    // one thread per cell
    const unsigned int n = blockIdx.x * blockDim.x + threadIdx.x;
    if (n >= num_fill_cells)
        return;

    const unsigned int cell = d_fill_cells[n];
    const unsigned int i = cell % local_dim.x;
    const unsigned int j = (cell / local_dim.x) % local_dim.y;
    const unsigned int k = cell / (local_dim.x * local_dim.y);

    const int gi = (int)i + origin.x;
    const int gj = (int)j + origin.y;
    const int gk = (int)k + origin.z;

    Scalar3 f_center = make_scalar3((gi + Scalar(0.5)) * inv_dim.x + grid_shift.x,
                                    (gj + Scalar(0.5)) * inv_dim.y + grid_shift.y,
                                    (gk + Scalar(0.5)) * inv_dim.z + grid_shift.z);
    f_center.x -= slow::floor(f_center.x);
    f_center.y -= slow::floor(f_center.y);
    f_center.z -= slow::floor(f_center.z);
    const Scalar3 cell_center = global_box.makeCoordinates(f_center);

    const Scalar3 half_w = make_scalar3(Scalar(0.5) * draw_box.getL().x,
                                        Scalar(0.5) * draw_box.getL().y,
                                        Scalar(0.5) * draw_box.getL().z);

    hoomd::RandomGenerator rng(
        hoomd::Seed(hoomd::RNGIdentifier::VirtualParticleFiller, timestep, seed),
        hoomd::Counter(globalCellIndex(gi, gj, gk, global_dim), filler_id, 0));

    unsigned int num_in_cell = hoomd::PoissonDistribution<Scalar>(mean_per_cell)(rng);
    if (num_in_cell > max_per_cell)
        num_in_cell = max_per_cell;

    const unsigned int base = n * max_per_cell;
    for (unsigned int p = 0; p < max_per_cell; ++p)
        {
        if (p >= num_in_cell)
            {
            d_keep_particles[base + p] = false;
            continue;
            }
        Scalar3 pos = make_scalar3(hoomd::UniformDistribution<Scalar>(-half_w.x, half_w.x)(rng),
                                   hoomd::UniformDistribution<Scalar>(-half_w.y, half_w.y)(rng),
                                   hoomd::UniformDistribution<Scalar>(-half_w.z, half_w.z)(rng));

        int3 img = make_int3(0, 0, 0);
        draw_box.wrap(pos, img);
        pos += cell_center;

        img = make_int3(0, 0, 0);
        global_box.wrap(pos, img);

        const bool is_outside = geom.isOutside(pos);
        d_keep_particles[base + p] = is_outside;

        if (!is_outside)
            continue;

        d_tmp_pos[base + p] = make_scalar4(pos.x, pos.y, pos.z, __int_as_scalar(type));

        hoomd::NormalDistribution<Scalar> gen(vel_factor, 0.0);
        Scalar3 vel;
        gen(vel.x, vel.y, rng);
        vel.z = gen(rng);
        geom.addToVirtualParticleVelocity(vel, pos);
        d_tmp_vel[base + p]
            = make_scalar4(vel.x, vel.y, vel.z, __int_as_scalar(mpcd::detail::NO_CELL));
        }
    }

    } // end namespace kernel

/*!
 * \param args Common arguments for all geometries
 * \param geom Confined geometry
 *
 * \tparam Geometry type of the confined geometry \a geom
 *
 * \sa mpcd::gpu::kernel::classify_cells
 */
template<class Geometry>
cudaError_t classify_cells(const classify_cells_args_t& args, const Geometry& geom)
    {
    const unsigned int num_cells = args.local_dim.x * args.local_dim.y * args.local_dim.z;
    if (num_cells == 0)
        return cudaSuccess;

    cudaFuncAttributes attr;
    cudaFuncGetAttributes(&attr, (const void*)mpcd::gpu::kernel::classify_cells<Geometry>);
    const unsigned int max_block_size = attr.maxThreadsPerBlock;

    unsigned int run_block_size = min(args.block_size, max_block_size);
    dim3 grid(num_cells / run_block_size + 1);
    mpcd::gpu::kernel::classify_cells<Geometry><<<grid, run_block_size>>>(args.d_flags,
                                                                          args.global_box,
                                                                          args.draw_box,
                                                                          args.local_dim,
                                                                          args.global_dim,
                                                                          args.origin,
                                                                          args.inv_dim,
                                                                          args.num_trials,
                                                                          args.seed,
                                                                          args.filler_id,
                                                                          geom);
    return cudaSuccess;
    }

/*!
 * \param args Common arguments for all geometries
 * \param geom Confined geometry
 *
 * \tparam Geometry type of the confined geometry \a geom
 *
 * \sa mpcd::gpu::kernel::draw_virtual_particles
 */
template<class Geometry>
cudaError_t draw_virtual_particles(const draw_virtual_particles_args_t& args, const Geometry& geom)
    {
    if (args.num_fill_cells == 0)
        return cudaSuccess;

    cudaFuncAttributes attr;
    cudaFuncGetAttributes(&attr, (const void*)mpcd::gpu::kernel::draw_virtual_particles<Geometry>);
    const unsigned int max_block_size = attr.maxThreadsPerBlock;

    unsigned int run_block_size = min(args.block_size, max_block_size);
    dim3 grid(args.num_fill_cells / run_block_size + 1);
    mpcd::gpu::kernel::draw_virtual_particles<Geometry>
        <<<grid, run_block_size>>>(args.d_tmp_pos,
                                   args.d_tmp_vel,
                                   args.d_keep_particles,
                                   args.d_fill_cells,
                                   args.num_fill_cells,
                                   args.global_box,
                                   args.draw_box,
                                   args.local_dim,
                                   args.global_dim,
                                   args.origin,
                                   args.inv_dim,
                                   args.grid_shift,
                                   args.mean_per_cell,
                                   args.max_per_cell,
                                   args.vel_factor,
                                   args.type,
                                   args.timestep,
                                   args.seed,
                                   args.filler_id,
                                   geom);
    return cudaSuccess;
    }

#endif // __HIPCC__

    } // end namespace gpu
    } // end namespace mpcd
    } // namespace hoomd

#endif // MPCD_REJECTION_VIRTUAL_PARTICLE_FILLER_GPU_CUH_
