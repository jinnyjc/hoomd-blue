// Copyright (c) 2009-2026 The Regents of the University of Michigan.
// Part of HOOMD-blue, released under the BSD 3-Clause License.

/*!
 * \file mpcd/TriangulatedGeometryFillerGPU.h
 * \brief Declaration of mpcd::TriangulatedGeometryFillerGPU
 */

#ifndef MPCD_TRIANGULATED_GEOMETRY_FILLER_GPU_H_
#define MPCD_TRIANGULATED_GEOMETRY_FILLER_GPU_H_

#ifdef __HIPCC__
#error This header cannot be compiled by nvcc
#endif

#include "TriangulatedGeometryFiller.h"
#include "TriangulatedGeometryFillerGPU.cuh"

namespace hoomd
    {
namespace mpcd
    {

//! Add virtual particles to MPCD particle data for a triangulated geometry using the GPU
class PYBIND11_EXPORT TriangulatedGeometryFillerGPU : public mpcd::TriangulatedGeometryFiller
    {
    public:
    //! Constructor
    TriangulatedGeometryFillerGPU(std::shared_ptr<SystemDefinition> sysdef,
                                  const std::string& type,
                                  Scalar density,
                                  std::shared_ptr<Variant> T,
                                  std::shared_ptr<const TriangulatedGeometry> geom,
                                  unsigned int num_classify_trials);

    //! Fill the volume outside the confinement
    void fill(uint64_t timestep) override;

    protected:
    //! Determine which cells need virtual particles
    void classifyCells() override;

    private:
    // Autotuners
    std::shared_ptr<Autotuner<1>> m_traverse_tuner; //!< Tuner for traversal
    std::shared_ptr<Autotuner<1>> m_classify_tuner; //!< Tuner for cell classification
    std::shared_ptr<Autotuner<2>> m_draw_tuner;     //!< Tuner for drawing particles
    std::shared_ptr<Autotuner<1>> m_tag_tuner;      //!< Tuner for particle tagging

    // LBVH
    std::unique_ptr<gpu::CellLBVHTraverserWrapper> m_traverser; //!< LBVH traverser

    GPUArray<unsigned int> m_fill_candidates; //!< Index of each fill cell in the candidate arrays
    GPUArray<bool> m_keep_particles;       //!< Track whether particles are in/out of the geometry
    GPUArray<unsigned int> m_keep_indices; //!< Indices for particles out the geometry
    GPUFlags<unsigned int> m_num_keep;     //!< Number of particles to keep
    GPUFlags<unsigned int> m_max_observed; //!< Largest per-cell draw observed by the kernel
    };

namespace detail
    {
//! Export mpcd::TriangulatedGeometryFillerGPU to python
/*!
 * \param m Python module to export to
 */
void export_TriangulatedGeometryFillerGPU(pybind11::module& m);
    } // end namespace detail
    } // end namespace mpcd
    } // end namespace hoomd
#endif // MPCD_TRIANGULATED_GEOMETRY_FILLER_GPU_H_
