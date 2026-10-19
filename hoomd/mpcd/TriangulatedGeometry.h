// Copyright (c) 2009-2026 The Regents of the University of Michigan.
// Part of HOOMD-blue, released under the BSD 3-Clause License.

/*!
 * \file mpcd/TriangulatedGeometry.h
 * \brief Declaration of the MPCD triangulated geometry
 */

#ifndef MPCD_TRIANGULATED_GEOMETRY_H_
#define MPCD_TRIANGULATED_GEOMETRY_H_

#include "hoomd/HOOMDMath.h"

#ifdef __HIPCC__
#define HOSTDEVICE __host__ __device__ inline
#else
#include "hoomd/AABBTree.h"
#include "hoomd/BoxDim.h"
#include "hoomd/ExecutionConfiguration.h"
#include "hoomd/GPUArray.h"
#include "hoomd/PythonLocalDataAccess.h"
#include "hoomd/SystemDefinition.h"
#ifdef ENABLE_HIP
#include "TriangulatedGeometry.cuh"
#endif // ENABLE_HIP
#define HOSTDEVICE inline __attribute__((always_inline))
#endif // __HIPCC__

namespace hoomd
    {
namespace mpcd
    {
//! Closest point on a triangle to a point
/*!
 * \param p Query point
 * \param a First vertex of the triangle
 * \param b Second vertex of the triangle
 * \param c Third vertex of the triangle
 * \returns The point on triangle abc (including its edges and vertices) closest to p
 *
 * The plane containing the triangle is partitioned into seven Voronoi regions (three
 * vertex regions, three edge regions, and the face). The region containing p is found
 * from the signs of six dot products, and p is projected accordingly: onto a vertex,
 * onto an edge segment, or onto the face using barycentric coordinates.
 *
 * This implementation follows the method described in:
 * Ericson, C. (2005). Real-Time Collision Detection. Section 5.1.5.
 */
HOSTDEVICE Scalar3 closestPointOnTriangle(const Scalar3& p,
                                          const Scalar3& a,
                                          const Scalar3& b,
                                          const Scalar3& c)
    {
    const Scalar3 ab = b - a;
    const Scalar3 ac = c - a;
    const Scalar3 ap = p - a;

    // vertex region a
    const Scalar d1 = dot(ab, ap);
    const Scalar d2 = dot(ac, ap);
    if (d1 <= Scalar(0.0) && d2 <= Scalar(0.0))
        return a;

    // vertex region b
    const Scalar3 bp = p - b;
    const Scalar d3 = dot(ab, bp);
    const Scalar d4 = dot(ac, bp);
    if (d3 >= Scalar(0.0) && d4 <= d3)
        return b;

    // vertex region c
    const Scalar3 cp = p - c;
    const Scalar d5 = dot(ab, cp);
    const Scalar d6 = dot(ac, cp);
    if (d6 >= Scalar(0.0) && d5 <= d6)
        return c;

    // edge region ab
    const Scalar vc = d1 * d4 - d3 * d2;
    if (vc <= Scalar(0.0) && d1 >= Scalar(0.0) && d3 <= Scalar(0.0))
        {
        const Scalar v = d1 / (d1 - d3);
        return a + v * ab;
        }

    // edge region ac
    const Scalar vb = d5 * d2 - d1 * d6;
    if (vb <= Scalar(0.0) && d2 >= Scalar(0.0) && d6 <= Scalar(0.0))
        {
        const Scalar w = d2 / (d2 - d6);
        return a + w * ac;
        }

    // edge region bc
    const Scalar va = d3 * d6 - d5 * d4;
    if (va <= Scalar(0.0) && (d4 - d3) >= Scalar(0.0) && (d5 - d6) >= Scalar(0.0))
        {
        const Scalar w = (d4 - d3) / ((d4 - d3) + (d5 - d6));
        return b + w * (c - b);
        }

    // face region: interpolate with normalized barycentric coordinates
    const Scalar denom = Scalar(1.0) / (va + vb + vc);
    const Scalar v = vb * denom;
    const Scalar w = vc * denom;
    return a + v * ab + w * ac;
    }

//! Test whether a point lies outside the fluid using a list of candidate triangles
/*!
 * \param p Query point
 * \param candidate_tris Candidate triangle indices
 * \param num_candidate_tris Number of candidates
 * \param vertices Triangle vertices
 * \param triangles Triangle vertex indices
 * \returns True if p lies on the side that the normal of the nearest candidate points to
 *
 * The sign of the dot product between (p - q) and the triangle normal decides the side,
 * where q is the closest point on the nearest candidate. Normals point out of the fluid,
 * and a point exactly on the surface is treated as fluid.
 */
HOSTDEVICE bool isOutside(const Scalar3& p,
                          const unsigned int* candidate_tris,
                          unsigned int num_candidate_tris,
                          const ShortReal3* vertices,
                          const uint3* triangles)
    {
    Scalar best_d2(0);
    Scalar best_sign(0);
    for (unsigned int t = 0; t < num_candidate_tris; ++t)
        {
        const uint3 tri = triangles[candidate_tris[t]];

        const ShortReal3 a = vertices[tri.x];
        const ShortReal3 b = vertices[tri.y];
        const ShortReal3 c = vertices[tri.z];
        const Scalar3 aa = make_scalar3(a.x, a.y, a.z);
        const Scalar3 bb = make_scalar3(b.x, b.y, b.z);
        const Scalar3 cc = make_scalar3(c.x, c.y, c.z);

        const Scalar3 q = closestPointOnTriangle(p, aa, bb, cc);
        const Scalar3 d = p - q;
        const Scalar d2 = dot(d, d);
        if (t == 0 || d2 < best_d2)
            {
            best_d2 = d2;
            const Scalar3 n = cross(bb - aa, cc - aa);
            best_sign = dot(d, n);
            }
        }
    return best_sign > Scalar(0.0);
    }

#ifndef __HIPCC__
//! Triangulated geometry
class PYBIND11_EXPORT TriangulatedGeometry
    {
    public:
    //! Constructor
    /*!
     * \param sysdef System definition
     * \param num_vertices Number of vertices
     * \param vertices Triangle vertices
     * \param num_triangles Number of triangles
     * \param triangles Triangle vertex indices
     * \param unwrap_distance Distance to unwrap the triangles
     * \param no_slip Boundary condition at the wall (slip or no-slip)
     */
    TriangulatedGeometry(std::shared_ptr<SystemDefinition> sysdef,
                         unsigned int num_vertices,
                         const Scalar3* vertices,
                         unsigned int num_triangles,
                         const uint3* triangles,
                         const Scalar unwrap_distance,
                         bool no_slip);

    //! Destructor
    ~TriangulatedGeometry();

    //! Get the number of original vertices
    unsigned int getNumVertices() const;

    //! Get the number of original triangles
    unsigned int getNumTriangles() const;

    //! Get the number of total vertices (original + unwrapped)
    unsigned int getNumTotalVertices() const;

    //! Get the number of total triangles (original + unwrapped)
    unsigned int getNumTotalTriangles() const;

    //! Get the triangle vertices
    const GPUArray<ShortReal3>& getVertices() const;

    //! Get the triangle vertex indices
    const GPUArray<uint3>& getTriangles() const;

    //! Get the unwrap distance
    const Scalar getUnwrapDistance() const;

    //! Get the wall boundary condition
    bool getNoSlip() const;

    //! Get triangle BVH Tree
    const hoomd::detail::AABBTree& getTriangleTree() const;
#ifdef ENABLE_HIP
    //! Get triangle LBVH Tree for GPU
    gpu::TriangleLBVHWrapper& getTriangleLBVH() const;
#endif

    private:
    std::shared_ptr<SystemDefinition> m_sysdef;                //!< System definition
    std::shared_ptr<const ExecutionConfiguration> m_exec_conf; //!< Execution configuration

    unsigned int m_num_vertices;  //!< Number of original vertices
    unsigned int m_num_triangles; //!< Number of original triangles

    unsigned int m_num_total_vertices;  //!< Number of total vertices
    unsigned int m_num_total_triangles; //!< Number of total triangles

    GPUArray<ShortReal3> m_vertices; //!< Vertex list
    GPUArray<uint3> m_triangles;     //!< Triangle list

    const Scalar m_unwrap_distance; //!< Distance used to unwrap triangles

    bool m_no_slip; //!< Boundary condition

    std::vector<hoomd::detail::AABB> m_triangle_aabbs; //!< AABBs for triangles
    hoomd::detail::AABBTree m_triangle_tree;           //!< BVH for triangles
#ifdef ENABLE_HIP
    std::unique_ptr<gpu::TriangleLBVHWrapper> m_lbvh; //!< LBVH for triangles
#endif                                                // ENABLE_HIP

    //! Unwrap the triangles in unwrap distance
    void unwrapTriangles(std::vector<Scalar3>& vertices, std::vector<uint3>& triangles);

    //! Build one AABB for each total triangle
    void buildTriangleAABBs();

    //! Build the BVH over all triangles
    void buildTree();
#ifdef ENABLE_HIP
    //! Build the LBVH over all triangles on the GPU
    void buildLBVH();
#endif // ENABLE_HIP
    };

//! Local data access to the triangle vertices and vertex indices
template<class Output>
class TriangulatedGeometryAccess : public LocalDataAccess<Output, TriangulatedGeometry>
    {
    public:
    TriangulatedGeometryAccess(TriangulatedGeometry& geometry, bool unwrapped = false)
        : LocalDataAccess<Output, TriangulatedGeometry>(geometry), m_vertices_handle(),
          m_triangles_handle(), m_unwrapped(unwrapped)
        {
        }

    virtual ~TriangulatedGeometryAccess() = default;

    Output getVertices()
        {
        const size_t n
            = m_unwrapped ? this->m_data.getNumTotalVertices() : this->m_data.getNumVertices();
        return this->template getBuffer<ShortReal3, ShortReal>(m_vertices_handle,
                                                               &TriangulatedGeometry::getVertices,
                                                               std::vector<size_t> {n, 3},
                                                               true);
        }

    Output getTriangles()
        {
        const size_t n
            = m_unwrapped ? this->m_data.getNumTotalTriangles() : this->m_data.getNumTriangles();
        return this->template getBuffer<uint3, uint>(m_triangles_handle,
                                                     &TriangulatedGeometry::getTriangles,
                                                     std::vector<size_t> {n, 3},
                                                     true);
        }

    protected:
    void clear()
        {
        m_vertices_handle.reset(nullptr);
        m_triangles_handle.reset(nullptr);
        }

    private:
    std::unique_ptr<ArrayHandle<ShortReal3>> m_vertices_handle;
    std::unique_ptr<ArrayHandle<uint3>> m_triangles_handle;
    bool m_unwrapped;
    };

namespace detail
    {
void export_TriangulatedGeometry(pybind11::module& m);
void export_TriangulatedGeometryAccessHost(pybind11::module& m);
#ifdef ENABLE_HIP
void export_TriangulatedGeometryAccessDevice(pybind11::module& m);
#endif // ENABLE_HIP

//! Export mpcd::TriangulatedGeometryAccess to python
template<class Output> void export_TriangulatedGeometryAccess(pybind11::module& m, std::string name)
    {
    pybind11::class_<TriangulatedGeometryAccess<Output>,
                     std::shared_ptr<TriangulatedGeometryAccess<Output>>>(m, name.c_str())
        .def(pybind11::init<TriangulatedGeometry&, bool>())
        .def("getVertices", &TriangulatedGeometryAccess<Output>::getVertices)
        .def("getTriangles", &TriangulatedGeometryAccess<Output>::getTriangles)
        .def("enter", &TriangulatedGeometryAccess<Output>::enter)
        .def("exit", &TriangulatedGeometryAccess<Output>::exit);
    }
    } // end namespace detail
#endif // __HIPCC__

    } // end namespace mpcd
    } // end namespace hoomd
#undef HOSTDEVICE

#endif // MPCD_TRIANGULATED_GEOMETRY_H_
