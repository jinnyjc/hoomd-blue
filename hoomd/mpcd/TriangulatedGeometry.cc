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
                                           const Scalar3* vertices,
                                           unsigned int num_triangles,
                                           const uint3* triangles,
                                           const Scalar3 unwrap_distance,
                                           bool no_slip)
    : m_sysdef(sysdef), m_exec_conf(m_sysdef->getParticleData()->getExecConf()),
      m_num_vertices(num_vertices), m_num_triangles(num_triangles),
      m_vertices(m_num_vertices, m_exec_conf), m_triangles(m_num_triangles, m_exec_conf),
      m_unwrap_distance(unwrap_distance), m_unwrapped_vertices(0, m_exec_conf),
      m_unwrapped_triangles(0, m_exec_conf), m_no_slip(no_slip)
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
    
    unwrapTriangles();
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

const Scalar3 TriangulatedGeometry::getUnwrapDistance() const
    {
    return m_unwrap_distance;
    }

unsigned int TriangulatedGeometry::getNumUnwrappedVertices() const
    {
    return m_num_unwrapped_vertices;
    }

unsigned int TriangulatedGeometry::getNumUnwrappedTriangles() const
    {
    return m_num_unwrapped_triangles;
    }

const GPUArray<Scalar3>& TriangulatedGeometry::getUnwrappedVertices() const
    {
    return m_unwrapped_vertices;
    }

const GPUArray<uint3>& TriangulatedGeometry::getUnwrappedTriangles() const
    {
    return m_unwrapped_triangles;
    }

void TriangulatedGeometry::unwrapTriangles()
    {
    const BoxDim box = m_sysdef->getParticleData()->getBox();
    const uchar3 periodic = box.getPeriodic();

    const int nx = (periodic.x && m_unwrap_distance.x > Scalar(0.0)) ? 1 : 0;
    const int ny = (periodic.y && m_unwrap_distance.y > Scalar(0.0)) ? 1 : 0;
    const int nz = (periodic.z && m_unwrap_distance.z > Scalar(0.0)) ? 1 : 0;

    const unsigned int n_images = 
        static_cast<unsigned int>((2 * nx + 1) * (2 * ny + 1) * (2 * nz + 1));

    m_num_unwrapped_vertices = m_num_vertices * n_images;
    m_num_unwrapped_triangles = m_num_triangles * n_images;

    m_unwrapped_vertices.resize(m_num_unwrapped_vertices);
    m_unwrapped_triangles.resize(m_num_unwrapped_triangles);

    ArrayHandle<Scalar3> h_base_vertices(m_vertices, access_location::host, access_mode::read);
    ArrayHandle<uint3> h_base_triangles(m_triangles, access_location::host, access_mode::read);

    ArrayHandle<Scalar3> h_unwrapped_vertices(m_unwrapped_vertices, access_location::host, access_mode::overwrite);
    ArrayHandle<uint3> h_unwrapped_triangles(m_unwrapped_triangles, access_location::host, access_mode::overwrite);

    unsigned int img = 0;
    for (int i = -nx; i <= nx; ++i)
        {
        for (int j = -ny; j <= ny; ++j)
            {
            for (int k = -nz; k <= nz; ++k)
                {
                const int3 s = make_int3(i, j ,k);
                const Scalar3 shift = box.shift(make_scalar3(0, 0, 0), s);

                const unsigned int v_off = img * m_num_vertices;
                const unsigned int t_off = img * m_num_triangles;

                for (unsigned vi = 0; vi < m_num_vertices; ++vi)
                    {
                    const Scalar3 vv = h_base_vertices.data[vi];
                    h_unwrapped_vertices.data[v_off + vi] =
                        make_scalar3(vv.x + shift.x, vv.y +shift.y, vv.z +shift.z);
                    }

                for (unsigned ti = 0; ti < m_num_triangles; ++ti)
                    {
                    const uint3 tt = h_base_triangles.data[ti];
                    h_unwrapped_triangles.data[t_off + ti] =
                        make_uint3(tt.x + v_off, tt.y + v_off, tt.z + v_off);
                    }

                ++ img;
                }
            }
        }
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
               pybind11::object unwrap_distance_obj,
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

                
                pybind11::sequence seq = unwrap_distance_obj.cast<pybind11::sequence>();
                const Scalar3 unwrap_distance = make_scalar3(
                    seq[0].cast<Scalar>(),
                    seq[1].cast<Scalar>(),
                    seq[2].cast<Scalar>());

                return std::make_shared<TriangulatedGeometry>(sysdef,
                                                              num_vertices,
                                                              v.data(),
                                                              num_triangles,
                                                              t.data(),
                                                              unwrap_distance,
                                                              no_slip);
            }))
        .def_property_readonly("num_vertices", &TriangulatedGeometry::getNumVertices)
        .def_property_readonly("num_triangles", &TriangulatedGeometry::getNumTriangles)
        .def_property_readonly("no_slip", &TriangulatedGeometry::getNoSlip);
    }
    } // end namespace detail
    } // end namespace mpcd
    } // end namespace hoomd
