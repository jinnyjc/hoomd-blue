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
     */
    TriangulatedGeometry(std::shared_ptr<SystemDefinition> sysdef,
                         unsigned int num_vertices,
                         unsigned int num_triangles);

    //! Number of vertices
    unsigned int getNumVertices() const;

    //! Number of triangles
    unsigned int getNumTriangles() const;

    const GPUArray<Scalar3>& getVertices() const;
    const GPUArray<uint3>& getTriangles() const;

    private:
    std::shared_ptr<SystemDefinition> m_sysdef;
    std::shared_ptr<const ExecutionConfiguration> m_exec_conf;

    unsigned int m_num_vertices;
    unsigned int m_num_triangles;

    GPUArray<Scalar3> m_vertices;
    GPUArray<uint3> m_triangles;
    };

    } // end namespace mpcd
    } // end namespace hoomd
#endif // MPCD_TRIANGULATED_GEOMETRY_H_
