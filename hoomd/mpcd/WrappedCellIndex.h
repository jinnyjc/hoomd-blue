// Copyright (c) 2009-2026 The Regents of the University of Michigan.
// Part of HOOMD-blue, released under the BSD 3-Clause License.

/*!
 * \file mpcd/WrappedCellIndex.h
 * \brief Definition of mpcd::wrappedCellIndex
 */

#ifndef MPCD_WRAPPED_CELL_INDEX_H_
#define MPCD_WRAPPED_CELL_INDEX_H_

#include "hoomd/HOOMDMath.h"

#ifdef __HIPCC__
#define HOSTDEVICE __host__ __device__
#else
#define HOSTDEVICE
#endif // __HIPCC__

namespace hoomd
    {
namespace mpcd
    {
//! Flatten wrapped 3D cell coordinates into a 1D index
/*!
 * \param i Global cell coordinate along x
 * \param j Global cell coordinate along y
 * \param k Global cell coordinate along z
 * \param dim Global grid dimensions
 * \returns The flattened 1D index of the wrapped cell
 */
HOSTDEVICE inline unsigned int wrappedCellIndex(int i, int j, int k, const uint3& dim)
    {
    if (i < 0)
        i += static_cast<int>(dim.x);
    else if (i >= static_cast<int>(dim.x))
        i -= static_cast<int>(dim.x);
    if (j < 0)
        j += static_cast<int>(dim.y);
    else if (j >= static_cast<int>(dim.y))
        j -= static_cast<int>(dim.y);
    if (k < 0)
        k += static_cast<int>(dim.z);
    else if (k >= static_cast<int>(dim.z))
        k -= static_cast<int>(dim.z);

    return static_cast<unsigned int>(i)
           + dim.x * (static_cast<unsigned int>(j) + dim.y * static_cast<unsigned int>(k));
    }

    } // end namespace mpcd
    } // end namespace hoomd
#undef HOSTDEVICE

#endif // MPCD_WRAPPED_CELL_INDEX_H_
