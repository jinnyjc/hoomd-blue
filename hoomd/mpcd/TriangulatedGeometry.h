// Copyright (c) 2009-2025 The Regents of the University of Michigan.
// Part of HOOMD-blue, released under the BSD 3-Clause License.

/*!
 * \file mpcd/TriangulatedGeometry.h
 * \brief Definition of the MPCD triangulated geometry
 */

#ifndef MPCD_TRIANGULATED_GEOMETRY_H_
#define MPCD_TRIANGULATED_GEOMETRY_H_

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
     * \param num_triangles Number of triangles
     * \param no_slip Boundary condition at the wall (slip or no-slip)
     */
    TriangulatedGeometry(std::shared_ptr<SystemDefinition> sysdef,
                         unsigned int num_vertices,
                         unsigned int num_triangles,
                         bool no_slip);

    //! Constructor
    /*!
     * \param sysdef System definition
     * \param num_vertices Number of vertices
     * \param vertices List of vertices
     * \param num_triangles Number of triangles
     * \param triangles List of triangles
     * \param no_slip Boundary condition at the wall (slip or no-slip)
     */
    TriangulatedGeometry(std::shared_ptr<SystemDefinition> sysdef,
                         unsigned int num_vertices,
                         const Scalar3* vertices,
                         unsigned int num_triangles,
                         const uint3* triangles,
                         bool no_slip);

    //! Number of vertices
    unsigned int getNumVertices() const;

    //! Number of triangles
    unsigned int getNumTriangles() const;

    const GPUArray<Scalar3>& getVertices() const;
    const GPUArray<uint3>& getTriangles() const;

    bool getNoSlip() const;

    private:
    std::shared_ptr<SystemDefinition> m_sysdef;
    std::shared_ptr<const ExecutionConfiguration> m_exec_conf;

    unsigned int m_num_vertices;
    unsigned int m_num_triangles;

    GPUArray<Scalar3> m_vertices;
    GPUArray<uint3> m_triangles;

    bool m_no_slip;
    };

template<class Output>
class TriangulatedGeometryAccess : public LocalDataAccess<Output, TriangulatedGeometry>
    {
    public:
    TriangulatedGeometryAccess(TriangulatedGeometry& geometry)
        : LocalDataAccess<Output, TriangulatedGeometry>(geometry), m_vertices_handle(),
          m_triangles_handle()
        {
        }

    virtual ~TriangulatedGeometryAccess() = default;

    Output getVertices()
        {
        return this->template getBuffer<Scalar3, Scalar>(
            m_vertices_handle,
            &TriangulatedGeometry::getVertices,
            std::vector<size_t> {this->m_data.getNumVertices(), 3},
            false);
        }

    Output getTriangles()
        {
        return this->template getBuffer<uint3, uint>(
            m_triangles_handle,
            &TriangulatedGeometry::getTriangles,
            std::vector<size_t> {this->m_data.getNumTriangles(), 3},
            false);
        }

    protected:
    void clear()
        {
        m_vertices_handle.reset(nullptr);
        m_triangles_handle.reset(nullptr);
        }

    private:
    std::unique_ptr<ArrayHandle<Scalar3>> m_vertices_handle;
    std::unique_ptr<ArrayHandle<uint3>> m_triangles_handle;
    };

namespace detail
    {
void export_TriangulatedGeometry(pybind11::module& m);
;
/// Export local access
template<class Output> void export_TriangulatedGeometryAccess(pybind11::module& m, std::string name)
    {
    pybind11::class_<TriangulatedGeometryAccess<Output>,
                     std::shared_ptr<TriangulatedGeometryAccess<Output>>>(m, name.c_str())
        .def(pybind11::init<TriangulatedGeometry&>())
        .def("getVertices", &TriangulatedGeometryAccess<Output>::getVertices)
        .def("getTriangles", &TriangulatedGeometryAccess<Output>::getTriangles)
        .def("enter", &TriangulatedGeometryAccess<Output>::enter)
        .def("exit", &TriangulatedGeometryAccess<Output>::exit);
    }

    } // end namespace detail
    } // end namespace mpcd
    } // end namespace hoomd
#endif // MPCD_TRIANGULATED_GEOMETRY_H_
