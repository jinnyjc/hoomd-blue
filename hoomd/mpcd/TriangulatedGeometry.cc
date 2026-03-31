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
                                           const Scalar unwrap_distance,
                                           bool no_slip)
    : m_sysdef(sysdef), m_exec_conf(m_sysdef->getParticleData()->getExecConf()),
      m_num_vertices(num_vertices), m_num_triangles(num_triangles),
      m_num_total_vertices(num_vertices), m_num_total_triangles(num_triangles),
      m_vertices(m_num_vertices, m_exec_conf), m_triangles(m_num_triangles, m_exec_conf),
      m_unwrap_distance(unwrap_distance), m_no_slip(no_slip)
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

unsigned int TriangulatedGeometry::getNumTotalVertices() const
    {
    return m_num_total_vertices;
    }

unsigned int TriangulatedGeometry::getNumTotalTriangles() const
    {
    return m_num_total_triangles;
    }

const GPUArray<Scalar3>& TriangulatedGeometry::getVertices() const
    {
    return m_vertices;
    }

const GPUArray<uint3>& TriangulatedGeometry::getTriangles() const
    {
    return m_triangles;
    }

const Scalar TriangulatedGeometry::getUnwrapDistance() const
    {
    return m_unwrap_distance;
    }

void TriangulatedGeometry::unwrapTriangles()
    {
    const BoxDim box = m_sysdef->getParticleData()->getBox();
    const uchar3 periodic = box.getPeriodic();

    const int nx = (periodic.x && m_unwrap_distance > Scalar(0.0)) ? 1 : 0;
    const int ny = (periodic.y && m_unwrap_distance > Scalar(0.0)) ? 1 : 0;
    const int nz = (periodic.z && m_unwrap_distance > Scalar(0.0)) ? 1 : 0;

    std::vector<Scalar3> extra_vertices;
    std::vector<uint3> extra_triangles;

        {
        ArrayHandle<Scalar3> h_base_vertices(m_vertices, access_location::host, access_mode::read);
        ArrayHandle<uint3> h_base_triangles(m_triangles, access_location::host, access_mode::read);

        const Scalar3 ghost_width
            = make_scalar3(m_unwrap_distance, m_unwrap_distance, m_unwrap_distance);

        for (int i = -nx; i <= nx; ++i)
            {
            for (int j = -ny; j <= ny; ++j)
                {
                for (int k = -nz; k <= nz; ++k)
                    {
                    // skip the original image
                    if (i == 0 && j == 0 && k == 0)
                        continue;

                    const int3 image = make_int3(i, j, k);

                    // for each image, unwrap all the vertices and check if they are inside of
                    // the box expanded by unwrap_distance
                    std::vector<Scalar3> shifted_vertices(m_num_vertices);
                    std::vector<bool> in_buffer(m_num_vertices, false);

                    for (unsigned vi = 0; vi < m_num_vertices; ++vi)
                        {
                        const Scalar3 vv = h_base_vertices.data[vi];

                        const Scalar3 shifted = box.shift(vv, image);
                        shifted_vertices[vi] = shifted;

                        const Scalar3 f = box.makeFraction(shifted, ghost_width);

                        in_buffer[vi]
                            = (f.x >= Scalar(0.0) && f.x <= Scalar(1.0) && f.y >= Scalar(0.0)
                               && f.y <= Scalar(1.0) && f.z >= Scalar(0.0) && f.z <= Scalar(1.0));
                        }

                    // build a compact vertex list
                    std::map<unsigned int, unsigned int> vertex_map;
                    for (unsigned ti = 0; ti < m_num_triangles; ++ti)
                        {
                        const uint3 tt = h_base_triangles.data[ti];

                        // skip if none of the vertex is in buffer
                        if (!in_buffer[tt.x] && !in_buffer[tt.y] && !in_buffer[tt.z])
                            continue;

                        const unsigned int old_idx[3] = {tt.x, tt.y, tt.z};
                        unsigned int new_idx[3];

                        for (unsigned int m = 0; m < 3; ++m)
                            {
                            auto entry = vertex_map.find(old_idx[m]);
                            if (entry == vertex_map.end())
                                {
                                const unsigned int idx
                                    = static_cast<unsigned int>(extra_vertices.size());
                                extra_vertices.push_back(shifted_vertices[old_idx[m]]);
                                vertex_map[old_idx[m]] = idx;
                                new_idx[m] = idx;
                                }
                            else
                                {
                                new_idx[m] = entry->second;
                                }
                            }
                        extra_triangles.push_back(make_uint3(new_idx[0], new_idx[1], new_idx[2]));
                        }
                    }
                }
            }
        }
    m_num_total_vertices = m_num_vertices + static_cast<unsigned int>(extra_vertices.size());
    m_num_total_triangles = m_num_triangles + static_cast<unsigned int>(extra_triangles.size());

    m_vertices.resize(m_num_total_vertices);
    m_triangles.resize(m_num_total_triangles);

    // append extra unwrapped to original
    if (extra_vertices.size() > 0)
        {
        ArrayHandle<Scalar3> h_vertices(m_vertices, access_location::host, access_mode::overwrite);
        for (unsigned int vi = 0; vi < extra_vertices.size(); ++vi)
            {
            h_vertices.data[vi + m_num_vertices] = extra_vertices[vi];
            }
        }

    if (extra_triangles.size() > 0)
        {
        ArrayHandle<uint3> h_triangles(m_triangles, access_location::host, access_mode::overwrite);
        for (unsigned int ti = 0; ti < extra_triangles.size(); ++ti)
            {
            const uint3 tri = extra_triangles[ti];
            h_triangles.data[m_num_triangles + ti] = make_uint3(tri.x + m_num_vertices,
                                                                tri.y + m_num_vertices,
                                                                tri.z + m_num_vertices);
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
               Scalar unwrap_distance,
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
                                                              unwrap_distance,
                                                              no_slip);
            }))
        .def_property_readonly("num_vertices", &TriangulatedGeometry::getNumVertices)
        .def_property_readonly("num_triangles", &TriangulatedGeometry::getNumTriangles)
        .def_property_readonly("unwrap_distance", &TriangulatedGeometry::getUnwrapDistance)
        .def_property_readonly("no_slip", &TriangulatedGeometry::getNoSlip);
    }

void export_TriangulatedGeometryAccessHost(pybind11::module& m)
    {
    export_TriangulatedGeometryAccess<HOOMDHostBuffer>(m, "TriangulatedGeometryAccessHost");
    }

void export_TriangulatedGeometryAccessDevice(pybind11::module& m)
    {
    export_TriangulatedGeometryAccess<HOOMDDeviceBuffer>(m, "TriangulatedGeometryAccessDevice");
    }

    } // end namespace detail
    } // end namespace mpcd
    } // end namespace hoomd
