// Copyright (c) 2009-2026 The Regents of the University of Michigan.
// Part of HOOMD-blue, released under the BSD 3-Clause License.

/*!
 * \file mpcd/TriangulatedGeometry.h
 * \brief Definition of the MPCD triangulated geometry
 */

#ifndef MPCD_TRIANGULATED_GEOMETRY_H_
#define MPCD_TRIANGULATED_GEOMETRY_H_

#include "hoomd/AABBTree.h"
#include "hoomd/BoxDim.h"
#include "hoomd/ExecutionConfiguration.h"
#include "hoomd/GPUArray.h"
#include "hoomd/HOOMDMath.h"
#include "hoomd/PythonLocalDataAccess.h"
#include "hoomd/SystemDefinition.h"

#include <memory>

namespace hoomd
    {
namespace mpcd
    {
//! Triangulated geometry
class TriangulatedGeometry
    {
    public:
    //! Constructor
    /*!
     * \param sysdef System definition
     * \param num_vertices Number of vertices
     * \param vertices List of vertices
     * \param num_triangles Number of triangles
     * \param triangles List of triangles
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

    //! Get the number of original vertices
    unsigned int getNumVertices() const;

    //! Get the number of original triangles
    unsigned int getNumTriangles() const;

    //! Get the number of total vertices (original + unwrapped)
    unsigned int getNumTotalVertices() const;

    //! Get the number of total triangles (original + unwrapped)
    unsigned int getNumTotalTriangles() const;

    //! Get the vertex list
    const GPUArray<ShortReal3>& getVertices() const;

    //! Get the triangle list
    const GPUArray<uint3>& getTriangles() const;

    //! Get the unwrap distance
    const Scalar getUnwrapDistance() const;

    //! Get the wall boundary condition
    bool getNoSlip() const;

    //! Get triangle BVH Tree
    const hoomd::detail::AABBTree& getTriangleTree() const;

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

    //! Unwrap the triangles in unwrap distance
    void unwrapTriangles(std::vector<Scalar3>& vertices, std::vector<uint3>& triangles);

    //! build one AABB for each total triangle
    void buildTriangleAABBs();

    //! build the BVH over all triangles
    void buildTree();
    };

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
                                                               false);
        }

    Output getTriangles()
        {
        const size_t n
            = m_unwrapped ? this->m_data.getNumTotalTriangles() : this->m_data.getNumTriangles();
        return this->template getBuffer<uint3, uint>(m_triangles_handle,
                                                     &TriangulatedGeometry::getTriangles,
                                                     std::vector<size_t> {n, 3},
                                                     false);
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
void export_TriangulatedGeometryAccessDevice(pybind11::module& m);

//! Export local access
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
    } // end namespace mpcd
    } // end namespace hoomd
#endif // MPCD_TRIANGULATED_GEOMETRY_H_
