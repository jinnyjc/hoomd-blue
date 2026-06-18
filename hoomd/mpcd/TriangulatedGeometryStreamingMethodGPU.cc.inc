// Copyright (c) 2009-2025 The Regents of the University of Michigan.
// Part of HOOMD-blue, released under the BSD 3-Clause License.

/*!
 * \file mpcd/TriangulatedGeometryStreamingMethod@_force@GPU.cc
 * \brief Template instantation for TriangulatedGeometryStreamingMethodGPU with @_force@
 *
 * The force classes are filled in by CMake using configure_file(). See
 * mpcd/CMakeLists.txt for list of values that are used.
 */

// clang-format off
#include "hoomd/mpcd/@_force@.h"
#include "hoomd/mpcd/TriangulatedGeometryStreamingMethodGPU.h"

#define FORCE_CLASS @_force@
#define EXPORT_FUNCTION export_TriangulatedGeometryStreamingMethod@_force@GPU
// clang-format on

namespace hoomd
    {
namespace mpcd
    {

// Explicit instantiation of GPU class
template class TriangulatedGeometryStreamingMethodGPU<FORCE_CLASS>;

namespace detail
    {
void EXPORT_FUNCTION(pybind11::module& m)
    {
    export_TriangulatedGeometryStreamingMethodGPU<FORCE_CLASS>(m);
    }
    } // end namespace detail
    } // end namespace mpcd
    } // end namespace hoomd
