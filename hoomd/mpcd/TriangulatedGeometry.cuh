// Copyright (c) 2009-2026 The Regents of the University of Michigan.
// Part of HOOMD-blue, released under the BSD 3-Clause License.

#ifndef MPCD_TRIANGULATED_GEOMETRY_CUH_
#define MPCD_TRIANGULATED_GEOMETRY_CUH_

/*!
 * \file mpcd/TriangulatedGeometry.cuh
 * \brief Declaration of the GPU data structures for mpcd::TriangulatedGeometry
 */

#include <vector>

#include <hip/hip_runtime.h>

#include "hoomd/HOOMDMath.h"

// forward declaration
namespace neighbor
    {
class LBVH;
    } // namespace neighbor

namespace hoomd
    {
namespace mpcd
    {
namespace gpu
    {
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
    void setup(const ShortReal3* vertices,
               const uint3* triangles,
               unsigned int num_triangles,
               hipStream_t stream);

    //! Build the LBVH
    void build(const ShortReal3* vertices,
               const uint3* triangles,
               unsigned int num_triangles,
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

    } // end namespace gpu
    } // end namespace mpcd
    } // end namespace hoomd
#endif // MPCD_TRIANGULATED_GEOMETRY_CUH_
