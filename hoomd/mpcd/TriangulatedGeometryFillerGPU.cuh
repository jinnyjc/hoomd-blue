// Copyright (c) 2009-2026 The Regents of the University of Michigan.
// Part of HOOMD-blue, released under the BSD 3-Clause License.

#ifndef MPCD_TRIANGULATED_GEOMETRY_FILLER_GPU_CUH_
#define MPCD_TRIANGULATED_GEOMETRY_FILLER_GPU_CUH_

/*!
 * \file mpcd/TriangulatedGeometryFillerGPU.cuh
 * \brief Declaration of CUDA kernels for mpcd::TriangulatedGeometryFillerGPU
 */

#include <vector>

#include <hip/hip_runtime.h>

#include "RejectionVirtualParticleFillerGPU.cuh"
#include "TriangulatedGeometry.h"
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
//! Common arguments for classifying which cells need virtual particle filling
struct triangulated_classify_cells_args_t
    {
    //! Constructor
    triangulated_classify_cells_args_t(const BoxDim& _global_box,
                                       const BoxDim& _draw_box,
                                       const Index3D _ci,
                                       const uint3 _global_dim,
                                       const int3 _origin,
                                       const Scalar3 _inv_dim,
                                       const Scalar3 _max_shift,
                                       const unsigned int _num_classify_trials,
                                       const unsigned int _seed,
                                       const unsigned int _filler_id)
        : global_box(_global_box), draw_box(_draw_box), ci(_ci), global_dim(_global_dim),
          origin(_origin), inv_dim(_inv_dim), max_shift(_max_shift),
          num_classify_trials(_num_classify_trials), seed(_seed), filler_id(_filler_id)
        {
        }

    const BoxDim global_box;
    const BoxDim draw_box;
    const Index3D ci;
    const uint3 global_dim;
    const int3 origin;
    const Scalar3 inv_dim;
    const Scalar3 max_shift;
    const unsigned int num_classify_trials;
    const unsigned int seed;
    const unsigned int filler_id;
    };

//! Wrapper around neighbor::LBVHTraverser for triangulated geometry filler
/*!
 * The query volume of each cell is the AABB of the cell grown by a margin,
 * see CellQueryOp. This wrapper only exposes data types that are natively
 * supported in HOOMD so that all neighbor-specific templates and structs
 * can be handled only in CUDA code.
 */
class CellLBVHTraverserWrapper
    {
    public:
    //! Constructor
    CellLBVHTraverserWrapper();

    //! Destructor
    ~CellLBVHTraverserWrapper();

    //! Setup the LBVH traverser
    void setup(neighbor::LBVH& lbvh, hipStream_t stream);

    //! Traverse the LBVH
    void traverse(const triangulated_classify_cells_args_t& args,
                  const unsigned int* d_cells,
                  unsigned int num_cells,
                  unsigned int* d_num_tri_per_cell,
                  unsigned int* d_tri_list,
                  unsigned int max_tri_per_cell,
                  bool* d_candidate_flags,
                  unsigned int* d_max_tri_per_cell,
                  neighbor::LBVH& lbvh,
                  hipStream_t stream,
                  unsigned int block_size);

    //! Get the list of tunable parameters
    std::vector<unsigned int> getTunableParameters() const;

    private:
    neighbor::LBVHTraverser* trav_; //!< Underlying neighbor::LBVHTraverser
    };

hipError_t triangulated_classify_cells(const triangulated_classify_cells_args_t& args,
                                       bool* d_cut_flags,
                                       bool* d_center_outside,
                                       const unsigned int* d_tri_list,
                                       const unsigned int* d_num_tri_per_cell,
                                       unsigned int max_tri_per_cell,
                                       const unsigned int* d_candidate_cells,
                                       unsigned int num_candidate_cells,
                                       const ShortReal3* d_vertices,
                                       const uint3* d_triangles,
                                       unsigned int block_size);

hipError_t triangulated_flag_fill_cells(bool* d_fill_flags,
                                        const unsigned int* d_num_tri_per_cell,
                                        const unsigned int* d_cut_candidates,
                                        const bool* d_center_outside,
                                        unsigned int num_cut_cells,
                                        unsigned int block_size);

hipError_t compact_flagged_cells(void* d_tmp,
                                 size_t& tmp_bytes,
                                 const bool* d_flags,
                                 const unsigned int* d_cells,
                                 unsigned int num_cells,
                                 unsigned int* d_out_cells,
                                 unsigned int* d_out_indices,
                                 unsigned int* d_num_out);

hipError_t triangulated_draw_virtual_particles(const draw_virtual_particles_args_t& args,
                                               const unsigned int* d_fill_candidates,
                                               const unsigned int* d_tri_list,
                                               const unsigned int* d_num_tri_per_cell,
                                               unsigned int max_tri_per_cell,
                                               const ShortReal3* d_vertices,
                                               const uint3* d_triangles);
    } // end namespace gpu
    } // end namespace mpcd
    } // end namespace hoomd
#endif // MPCD_TRIANGULATED_GEOMETRY_FILLER_GPU_CUH_
