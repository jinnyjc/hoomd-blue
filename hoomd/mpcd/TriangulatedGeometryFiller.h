// Copyright (c) 2009-2026 The Regents of the University of Michigan.
// Part of HOOMD-blue, released under the BSD 3-Clause License.

/*!
 * \file mpcd/TriangulatedGeometryFiller.h
 * \brief Declaration of mpcd::TriangulatedGeometryFiller
 */

#ifndef MPCD_TRIANGULATED_GEOMETRY_FILLER_H_
#define MPCD_TRIANGULATED_GEOMETRY_FILLER_H_

#ifdef __HIPCC__
#error This header cannot be compiled by nvcc
#endif

#include "TriangulatedGeometry.h"
#include "VirtualParticleFiller.h"
#include "WrappedCellIndex.h"

#include <pybind11/pybind11.h>

namespace hoomd
    {
namespace mpcd
    {
//! Adds virtual particles to MPCD particle data for a triangulated geometry
/*!
 * Here we implement a virtual particle filler for a triangulated geometry using the rejection
 * sampling method. The filler first identifies the collision cells that are cut by the
 * triangulated geometry, i.e., the cells that can contain both fluid and solid, and then
 * only draws particles in those cells.
 *
 * To identify those cells, the AABB of each cell grown by the maximum grid shift (swept region)
 * is queried against the BVH of the triangles. A cell whose swept region holds no triangles
 * cannot be cut and is skipped. Otherwise, trial points are drawn in the swept region and tested
 * against those triangles. A point is classified by finding the nearest triangle and checking
 * the sign of the dot product between (p - q) and the triangle normal, where q is the
 * closest point on that triangle. Normals point out of the fluid, so a positive sign means that
 * point is outside the fluid. A cell with trial points on both sides needs filling.
 *
 * The triangles are static, so each fill cell stores its candidate triangles once during
 * classification and reuses them every fill step. Particles are drawn in the cell itself, so
 * only the triangles whose AABBs overlap the cell, not the swept region, are stored. A fill cell
 * with no candidates lies entirely on one side of the surface: it is kept if its center is outside
 * the fluid and dropped otherwise. The classification is redone only when the simulation box
 * changes.
 *
 * During filling, the number of particles drawn in each cell is Poisson distributed with mean
 * equal to the fill density times the cell volume. Each drawn particle is kept only if it lies
 * outside the fluid according to the candidate triangles of its cell, or without a test if the
 * cell has none. The per-cell capacity of the temporary arrays is grown whenever a draw exceeds
 * it.
 */
class PYBIND11_EXPORT TriangulatedGeometryFiller : public mpcd::VirtualParticleFiller
    {
    public:
    //! Constructor
    /*!
     * \param sysdef System definition
     * \param type Type of the virtual particles
     * \param density Fill density
     * \param T Temperature variant
     * \param geom Confining triangulated geometry
     * \param num_classify_trials Number of trial points drawn per cell
     */
    TriangulatedGeometryFiller(std::shared_ptr<SystemDefinition> sysdef,
                               const std::string& type,
                               Scalar density,
                               std::shared_ptr<Variant> T,
                               std::shared_ptr<const TriangulatedGeometry> geom,
                               unsigned int num_classify_trials);

    //! Destructor
    virtual ~TriangulatedGeometryFiller();

    //! Get the triangulated geometry
    std::shared_ptr<const TriangulatedGeometry> getGeometry() const
        {
        return m_geom;
        }

    //! Set the triangulated geometry
    void setGeometry(std::shared_ptr<const TriangulatedGeometry> geom)
        {
        m_geom = geom;
        m_need_classify = true;
        }

    //! Get the number of trial points used per cell during classification (0 = automatic)
    unsigned int getNumClassifyTrials() const
        {
        return m_num_classify_trials;
        }

    //! Get the number of cells that are filled
    unsigned int getNumFillCells() const
        {
        return m_num_fill_cells;
        }

    //! Get the list of cells that are filled
    const GPUArray<unsigned int>& getFillCells() const
        {
        return m_fill_cells;
        }

    //! Fill the particles outside the confinement
    void fill(uint64_t timestep) override;

    protected:
    std::shared_ptr<const TriangulatedGeometry> m_geom; //!< Confining triangulated geometry
    GPUArray<Scalar4> m_tmp_pos;                        //!< Temporary positions
    GPUArray<Scalar4> m_tmp_vel;                        //!< Temporary velocities

    unsigned int m_num_classify_trials;  //!< Trial points per cell during classification (0=auto)
    GPUArray<unsigned int> m_fill_cells; //!< Local 1D indices of the cells that need filling
    unsigned int m_num_fill_cells;       //!< Number of cells that need filling
    bool m_need_classify;                //!< True if the cell classification is out of date
    unsigned int m_alloc_per_cell;       //!< Per-cell particle capacity of the temporary arrays

    GPUArray<unsigned int> m_tri_list;         //!< Candidate triangles of each fill cell
    GPUArray<unsigned int> m_num_tri_per_cell; //!< Number of candidate triangles per fill cell
    unsigned int m_max_tri_per_cell;           //!< Uniform per-cell capacity of m_tri_list

    //! Flag the cell classification as out of date
    void setNeedClassify()
        {
        m_need_classify = true;
        }

    //! Determine which cells need virtual particles and store their triangles
    virtual void classifyCells();
    };

namespace detail
    {
//! Export TriangulatedGeometryFiller to python
void export_TriangulatedGeometryFiller(pybind11::module& m);
    } // end namespace detail
    } // end namespace mpcd
    } // end namespace hoomd
#endif // MPCD_TRIANGULATED_GEOMETRY_FILLER_H_
