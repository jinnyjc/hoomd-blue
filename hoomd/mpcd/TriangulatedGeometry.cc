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
                                           unsigned int num_triangles,
                                           bool no_slip)
    : m_sysdef(sysdef), m_exec_conf(m_sysdef->getParticleData()->getExecConf()),
      m_num_vertices(num_vertices), m_num_triangles(num_triangles),
      m_vertices(m_num_vertices, m_exec_conf), m_triangles(m_num_triangles, m_exec_conf),
      m_no_slip(no_slip)
    {
    }

TriangulatedGeometry::TriangulatedGeometry(std::shared_ptr<SystemDefinition> sysdef,
                                           unsigned int num_vertices,
                                           const Scalar3* vertices,
                                           unsigned int num_triangles,
                                           const uint3* triangles,
                                           bool no_slip)
    : m_sysdef(sysdef), m_exec_conf(m_sysdef->getParticleData()->getExecConf()),
      m_num_vertices(num_vertices), m_num_triangles(num_triangles),
      m_vertices(m_num_vertices, m_exec_conf), m_triangles(m_num_triangles, m_exec_conf),
      m_no_slip(no_slip)
    {
    if (m_num_vertices > 0)
        {
        ArrayHandle<Scalar3> h_vertices(m_vertices, access_location::host, access_mode::overwrite);
        std::copy(vertices, vertices + m_num_vertices, h_vertices.data);
        }

    if (m_num_triangles > 0)
        {
        ArrayHandle<uint3> h_triangles(m_triangles, access_location::host, access_mode::overwrite);
        std::copy(triangles, triangles + m_num_triangles, h_triangles.data);
        }
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

bool TriangulatedGeometry::getNoSlip() const
    {
    return m_no_slip;
    }

namespace detail
    {
void export_TriangulatedGeometry(pybind11::module& m)
    {
    pybind11::class_<TriangulatedGeometry, std::shared_ptr<TriangulatedGeometry>>(
        m,
        "TriangulatedGeometry")
        .def(pybind11::init(
            [](std::shared_ptr<SystemDefinition> sysdef,
               pybind11::array_t<Scalar, pybind11::array::c_style | pybind11::array::forcecast>
                   vertices,
               pybind11::array_t<unsigned int,
                                 pybind11::array::c_style | pybind11::array::forcecast> triangles,
               bool no_slip)
            {
                if (vertices.shape(1) != 3)
                    {
                    throw std::runtime_error("Vertices must have shape (N, 3)");
                    }

                if (triangles.shape(1) != 3)
                    {
                    throw std::runtime_error("Triangles must have shape (N, 3)");
                    }

                unsigned int num_vertices = static_cast<unsigned int>(vertices.shape(0));
                unsigned int num_triangles = static_cast<unsigned int>(triangles.shape(0));

                std::vector<Scalar3> v(num_vertices);
                auto v_in = vertices.unchecked<2>();
                for (unsigned int i = 0; i < num_vertices; ++i)
                    v[i] = make_scalar3(v_in(i, 0), v_in(i, 1), v_in(i, 2));

                std::vector<uint3> t(num_triangles);
                auto t_in = triangles.unchecked<2>();
                for (unsigned int i = 0; i < num_triangles; ++i)
                    t[i] = make_uint3(t_in(i, 0), t_in(i, 1), t_in(i, 2));

                return std::make_shared<TriangulatedGeometry>(sysdef,
                                                              num_vertices,
                                                              v.data(),
                                                              num_triangles,
                                                              t.data(),
                                                              no_slip);
            }))
        .def_property_readonly("num_vertices", &TriangulatedGeometry::getNumVertices)
        .def_property_readonly("num_triangles", &TriangulatedGeometry::getNumTriangles)
        .def_property_readonly("no_slip", &TriangulatedGeometry::getNoSlip);
    }
    } // end namespace detail
    } // end namespace mpcd
    } // end namespace hoomd
