// Copyright (c) 2009-2025 The Regents of the University of Michigan.
// Part of HOOMD-blue, released under the BSD 3-Clause License.

/*!
 * \file mpcd/TriangulatedGeometryStreamingMethod@_force@.cc
 * \brief Template instantation for TriangulatedGeometryStreamingMethod with @_force@
 *
 * The geometry and force classes are filled in by CMake using configure_file(). See
 * mpcd/CMakeLists.txt for list of values that are used.
 */

// clang-format off
#include "hoomd/mpcd/@_force@.h"
#include "hoomd/mpcd/TriangulatedGeometryStreamingMethod.h"

#define FORCE_CLASS @_force@
#define EXPORT_FUNCTION export_TriangulatedGeometryStreamingMethod@_force@
// clang-format on

namespace hoomd
    {
namespace mpcd
    {

// Explicit instantiation of CPU class
template class TriangulatedGeometryStreamingMethod<FORCE_CLASS>;

namespace detail
    {
void EXPORT_FUNCTION(pybind11::module& m)
    {
    export_TriangulatedGeometryStreamingMethod<FORCE_CLASS>(m);
    }
    } // end namespace detail
    } // end namespace mpcd
    } // end namespace hoomd
