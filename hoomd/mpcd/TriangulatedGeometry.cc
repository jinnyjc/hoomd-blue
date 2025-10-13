// Copyright (c) 2009-2025 The Regents of the University of Michigan.
// Part of HOOMD-blue, released under the BSD 3-Clause License.

/*!
 * \file mpcd/TriangulatedGeometry.cc
 * \brief Export function MPCD triangulated geometry.
 */

#include "TriangulatedGeometry.h"

namespace hoomd
    {
namespace mpcd
    {

TriangulatedGeometry::TriangulatedGeometry(std::shared_ptr<SystemDefinition> sysdef,
                                           unsigned int num_vertices,
                                           unsigned int num_triangles)
    : m_sysdef(sysdef), m_exec_conf(m_sysdef->getParticleData()->getExecConf()),
      m_num_vertices(num_vertices), m_num_triangles(num_triangles),
      m_vertices(m_num_vertices, m_exec_conf), m_triangles(m_num_triangles, m_exec_conf)
    {
    }

unsigned int TriangulatedGeometry::getNumVertices() const
    {
    return m_num_vertices;
    }

unsigned int TriangulatedGeometry::getNumTriangles() const
    {
    return m_num_triangles;
    }

const GPUArray<Scalar3>& TriangulatedGeometry::getVertices() const
    {
    return m_vertices;
    }

const GPUArray<uint3>& TriangulatedGeometry::getTriangles() const
    {
    return m_triangles;
    }

namespace detail
    {
void export_TriangulatedGeometry(pybind11::module& m)
    {
    pybind11::class_<TriangulatedGeometry, std::shared_ptr<TriangulatedGeometry>>(
        m,
        "TriangulatedGeometry")
        .def(pybind11::init<std::shared_ptr<SystemDefinition>, unsigned int, unsigned int>())
        .def_property_readonly("num_vertices", &TriangulatedGeometry::getNumVertices)
        .def_property_readonly("num_triangles", &TriangulatedGeometry::getNumTriangles);
    }
    } // end namespace detail
    } // end namespace mpcd
    } // end namespace hoomd
