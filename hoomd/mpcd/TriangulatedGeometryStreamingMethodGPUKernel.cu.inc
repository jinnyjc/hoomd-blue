// Copyright (c) 2009-2025 The Regents of the University of Michigan.
// Part of HOOMD-blue, released under the BSD 3-Clause License.

/*!
 * \file mpcd/TriangulatedGeometryStreamingMethod@_force@GPU.cu
 * \brief Template instantation for TriangulatedGeometryStreamingMethodGPU driver (and so kernel)
 * with @_force@.
 *
 * The force class is filled in by CMake using configure_file(). See
 * mpcd/CMakeLists.txt for list of values that are used.
 */

// clang-format off
#include "hoomd/mpcd/@_force@.h"
#include "hoomd/mpcd/TriangulatedGeometryStreamingMethodGPU.cuh"

#define FORCE_CLASS @_force@
// clang-format on

namespace hoomd
    {
namespace mpcd
    {
namespace gpu
    {
//! Template instantiation
template hipError_t __attribute__((visibility("default")))
triangulated_stream_init<FORCE_CLASS>(const triangulated_stream_args_t& args,
                                      const FORCE_CLASS& solvent_force,
                                      Scalar* d_dt_remain,
                                      unsigned int* d_active_idx,
                                      unsigned int block_size);

template hipError_t __attribute__((visibility("default")))
triangulated_stream_finalize<FORCE_CLASS>(const triangulated_stream_args_t& args,
                                          const FORCE_CLASS& solvent_force,
                                          unsigned int block_size);

    } // end namespace gpu
    } // end namespace mpcd
    } // end namespace hoomd
