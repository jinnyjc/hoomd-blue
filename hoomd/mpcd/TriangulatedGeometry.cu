// Copyright (c) 2009-2026 The Regents of the University of Michigan.
// Part of HOOMD-blue, released under the BSD 3-Clause License.

/*!
 * \file mpcd/TriangulatedGeometry.cu
 * \brief Definition of the GPU data structures for mpcd::TriangulatedGeometry
 */

#include "TriangulatedGeometry.cuh"

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
     * \param vertices_ Triangle vertices
     * \param triangles_ Triangle vertex indices
     * \param N_ Number of primitives (triangles) to insert
     */
    TriangleInsertOp(const ShortReal3* vertices_, const uint3* triangles_, unsigned int N_)
        : vertices(vertices_), triangles(triangles_), N(N_)
        {
        }

    //! Construct bounding box for a triangle
    /*!
     * \param idx Nominal index of the primitive [0,N)
     * \returns A neighbor::BoundingBox covering the three vertices
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
 * Initializes the underlying LBVH.
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

    } // end namespace gpu
    } // end namespace mpcd
    } // end namespace hoomd
