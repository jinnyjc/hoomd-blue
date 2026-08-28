// Copyright (c) 2009-2026 The Regents of the University of Michigan.
// Part of HOOMD-blue, released under the BSD 3-Clause License.

/*!
 * \file mpcd/RejectionVirtualParticleFillerGPU.h
 * \brief Declaration of mpcd::RejectionVirtualParticleFillerGPU
 */

#ifndef MPCD_REJECTION_VIRTUAL_PARTICLE_FILLER_GPU_H_
#define MPCD_REJECTION_VIRTUAL_PARTICLE_FILLER_GPU_H_

#ifdef __HIPCC__
#error This header cannot be compiled by nvcc
#endif

#include "RejectionVirtualParticleFiller.h"
#include "RejectionVirtualParticleFillerGPU.cuh"

#include "hoomd/Autotuner.h"
#include "hoomd/CachedAllocator.h"
#include <pybind11/pybind11.h>

#include <iterator>

namespace hoomd
    {
namespace mpcd
    {

//! Adds virtual particles to the MPCD particle data for various confining geometries using the GPU
template<class Geometry>
class PYBIND11_EXPORT RejectionVirtualParticleFillerGPU
    : public mpcd::RejectionVirtualParticleFiller<Geometry>
    {
    public:
    //! Constructor
    RejectionVirtualParticleFillerGPU(std::shared_ptr<SystemDefinition> sysdef,
                                      const std::string& type,
                                      Scalar density,
                                      std::shared_ptr<Variant> T,
                                      std::shared_ptr<const Geometry> geom,
                                      unsigned int num_trials,
                                      unsigned int max_per_cell)
        : mpcd::RejectionVirtualParticleFiller<Geometry>(sysdef,
                                                         type,
                                                         density,
                                                         T,
                                                         geom,
                                                         num_trials,
                                                         max_per_cell),
          m_cell_flags(this->m_exec_conf), m_keep_particles(this->m_exec_conf),
          m_keep_indices(this->m_exec_conf), m_num_keep(this->m_exec_conf)
        {
        m_tuner1.reset(new Autotuner<1>({AutotunerBase::makeBlockSizeRange(this->m_exec_conf)},
                                        this->m_exec_conf,
                                        "mpcd_rejection_filler_classify_cells"));
        m_tuner2.reset(new Autotuner<2>({AutotunerBase::makeBlockSizeRange(this->m_exec_conf),
                                         AutotunerBase::getTppListPow2(this->m_exec_conf)},
                                        this->m_exec_conf,
                                        "mpcd_rejection_filler_draw_particles"));
        m_tuner3.reset(new Autotuner<1>({AutotunerBase::makeBlockSizeRange(this->m_exec_conf)},
                                        this->m_exec_conf,
                                        "mpcd_rejection_filler_tag_particles"));
        this->m_autotuners.insert(this->m_autotuners.end(), {m_tuner1, m_tuner2, m_tuner3});
        }

    protected:
    //! Determine which cells need virtual particles
    void classifyCells() override;

    //! Fill the volume outside the confinement
    void fill(uint64_t timestep) override;

    private:
    GPUArray<bool> m_cell_flags;     //!< One flag per cell, set if the cell must be filled
    GPUArray<bool> m_keep_particles; //!< Track whether particles are in/out of bounds for geometry
    GPUArray<unsigned int> m_keep_indices;  //!< Indices for particles out of bounds for geometry
    GPUFlags<unsigned int> m_num_keep;      //!< Number of particles to keep
    std::shared_ptr<Autotuner<1>> m_tuner1; //!< Autotuner for cell classification
    std::shared_ptr<Autotuner<2>> m_tuner2; //!< Autotuner for drawing particles
    std::shared_ptr<Autotuner<1>> m_tuner3; //!< Autotuner for particle tagging
    };

template<class Geometry> void RejectionVirtualParticleFillerGPU<Geometry>::classifyCells()
    {
    this->m_cl->computeDimensions();

    const BoxDim& global_box = this->m_pdata->getGlobalBox();
    const Scalar3 global_L = global_box.getL();

    const uint3 local_dim = this->m_cl->getDim();
    const uint3 global_dim = this->m_cl->getGlobalDim();
    const int3 origin = this->m_cl->getOriginIndex();
    const unsigned int num_cells = local_dim.x * local_dim.y * local_dim.z;

    const Scalar3 inv_dim = make_scalar3(Scalar(1.0) / global_dim.x,
                                         Scalar(1.0) / global_dim.y,
                                         Scalar(1.0) / global_dim.z);

    const Scalar3 max_shift
        = this->m_cl->isGridShifting() ? this->m_cl->getMaxGridShift() : make_scalar3(0, 0, 0);

    const Scalar3 width = make_scalar3(inv_dim.x + Scalar(2.0) * max_shift.x,
                                       inv_dim.y + Scalar(2.0) * max_shift.y,
                                       inv_dim.z + Scalar(2.0) * max_shift.z);
    BoxDim draw_box(make_scalar3(width.x * global_L.x, width.y * global_L.y, width.z * global_L.z));
    draw_box.setTiltFactors(global_box.getTiltFactorXY(),
                            global_box.getTiltFactorXZ(),
                            global_box.getTiltFactorYZ());

    if (num_cells > m_cell_flags.getNumElements())
        {
        GPUArray<bool> cell_flags(num_cells, this->m_exec_conf);
        m_cell_flags.swap(cell_flags);
        }

    if (num_cells > this->m_fill_cells.getNumElements())
        {
        GPUArray<unsigned int> fill_cells(num_cells, this->m_exec_conf);
        this->m_fill_cells.swap(fill_cells);
        }

        {
        ArrayHandle<bool> d_cell_flags(m_cell_flags,
                                       access_location::device,
                                       access_mode::overwrite);
        ArrayHandle<unsigned int> d_fill_cells(this->m_fill_cells,
                                               access_location::device,
                                               access_mode::overwrite);

        mpcd::gpu::classify_cells_args_t args(d_cell_flags.data,
                                              global_box,
                                              draw_box,
                                              local_dim,
                                              global_dim,
                                              origin,
                                              inv_dim,
                                              this->m_num_trials,
                                              this->m_sysdef->getSeed(),
                                              this->m_filler_id,
                                              m_tuner1->getParam()[0]);
        m_tuner1->begin();
        mpcd::gpu::classify_cells<Geometry>(args, *(this->m_geom));
        if (this->m_exec_conf->isCUDAErrorCheckingEnabled())
            CHECK_CUDA_ERROR();
        m_tuner1->end();

            {
            // compact the flagged cells into a list of indices with CUB
            void* d_tmp_storage = NULL;
            size_t tmp_storage_bytes = 0;
            mpcd::gpu::compact_virtual_particle_indices(d_tmp_storage,
                                                        tmp_storage_bytes,
                                                        d_cell_flags.data,
                                                        num_cells,
                                                        d_fill_cells.data,
                                                        m_num_keep.getDeviceFlags());
            ScopedAllocation<unsigned char> d_tmp_alloc(this->m_exec_conf->getCachedAllocator(),
                                                        (tmp_storage_bytes > 0) ? tmp_storage_bytes
                                                                                : 1);
            d_tmp_storage = (void*)d_tmp_alloc();

            // run selection
            mpcd::gpu::compact_virtual_particle_indices(d_tmp_storage,
                                                        tmp_storage_bytes,
                                                        d_cell_flags.data,
                                                        num_cells,
                                                        d_fill_cells.data,
                                                        m_num_keep.getDeviceFlags());
            }
        }
    this->m_num_fill_cells = m_num_keep.readFlags();

    this->m_exec_conf->msg->notice(6)
        << "MPCD RejectionVirtualParticleFiller: filling " << this->m_num_fill_cells << " of "
        << num_cells << " cells" << std::endl;
    }

template<class Geometry> void RejectionVirtualParticleFillerGPU<Geometry>::fill(uint64_t timestep)
    {
    this->m_cl->computeDimensions();

    if (this->m_need_classify)
        {
        classifyCells();
        this->m_need_classify = false;
        }

    const BoxDim& global_box = this->m_pdata->getGlobalBox();
    const Scalar3 global_L = global_box.getL();

    const uint3 local_dim = this->m_cl->getDim();
    const uint3 global_dim = this->m_cl->getGlobalDim();
    const int3 origin = this->m_cl->getOriginIndex();

    const Scalar3 inv_dim = make_scalar3(Scalar(1.0) / global_dim.x,
                                         Scalar(1.0) / global_dim.y,
                                         Scalar(1.0) / global_dim.z);
    const Scalar3 grid_shift = this->m_cl->getGridShift();

    BoxDim draw_box(
        make_scalar3(inv_dim.x * global_L.x, inv_dim.y * global_L.y, inv_dim.z * global_L.z));
    draw_box.setTiltFactors(global_box.getTiltFactorXY(),
                            global_box.getTiltFactorXZ(),
                            global_box.getTiltFactorYZ());

    const Scalar cell_volume
        = global_box.getVolume() / (global_dim.x * global_dim.y * global_dim.z);
    const Scalar mean_per_cell = this->m_density * cell_volume;
    const unsigned int max_per_cell
        = (this->m_max_per_cell > 0) ? this->m_max_per_cell
                                     : static_cast<unsigned int>(std::ceil(
                                           mean_per_cell + Scalar(8.0) * std::sqrt(mean_per_cell)));

    // Step 1: size the temporary arrays with a fixed block of max_per_cell entries per fill cell
    const unsigned int num_virtual_max = this->m_num_fill_cells * max_per_cell;
    if (num_virtual_max > this->m_tmp_pos.getNumElements())
        {
        GPUArray<Scalar4> tmp_pos(num_virtual_max, this->m_exec_conf);
        this->m_tmp_pos.swap(tmp_pos);
        GPUArray<Scalar4> tmp_vel(num_virtual_max, this->m_exec_conf);
        this->m_tmp_vel.swap(tmp_vel);
        GPUArray<bool> keep_particles(num_virtual_max, this->m_exec_conf);
        m_keep_particles.swap(keep_particles);
        GPUArray<unsigned int> keep_indices(num_virtual_max, this->m_exec_conf);
        m_keep_indices.swap(keep_indices);
        }

    // Step 2: draw a Poisson number of particles in each cell that needs filling, then keep only
    // the ones that are outside the geometry
    const Scalar vel_factor = fast::sqrt((*this->m_T)(timestep) / this->m_mpcd_pdata->getMass());
    unsigned int num_selected = 0;
        {
        ArrayHandle<Scalar4> d_tmp_pos(this->m_tmp_pos,
                                       access_location::device,
                                       access_mode::overwrite);
        ArrayHandle<Scalar4> d_tmp_vel(this->m_tmp_vel,
                                       access_location::device,
                                       access_mode::overwrite);
        ArrayHandle<bool> d_keep_particles(m_keep_particles,
                                           access_location::device,
                                           access_mode::overwrite);
        ArrayHandle<unsigned int> d_keep_indices(m_keep_indices,
                                                 access_location::device,
                                                 access_mode::overwrite);
        ArrayHandle<unsigned int> d_fill_cells(this->m_fill_cells,
                                               access_location::device,
                                               access_mode::read);

        mpcd::gpu::draw_virtual_particles_args_t args(d_tmp_pos.data,
                                                      d_tmp_vel.data,
                                                      d_keep_particles.data,
                                                      d_fill_cells.data,
                                                      this->m_num_fill_cells,
                                                      global_box,
                                                      draw_box,
                                                      local_dim,
                                                      global_dim,
                                                      origin,
                                                      inv_dim,
                                                      grid_shift,
                                                      mean_per_cell,
                                                      max_per_cell,
                                                      vel_factor,
                                                      this->m_type,
                                                      timestep,
                                                      this->m_sysdef->getSeed(),
                                                      this->m_filler_id,
                                                      m_tuner2->getParam()[0],
                                                      m_tuner2->getParam()[1]);
        m_tuner2->begin();
        mpcd::gpu::draw_virtual_particles<Geometry>(args, *(this->m_geom));
        if (this->m_exec_conf->isCUDAErrorCheckingEnabled())
            CHECK_CUDA_ERROR();
        m_tuner2->end();

            {
            // compact the selected particles down with CUB
            void* d_tmp_storage = NULL;
            size_t tmp_storage_bytes = 0;
            mpcd::gpu::compact_virtual_particle_indices(d_tmp_storage,
                                                        tmp_storage_bytes,
                                                        d_keep_particles.data,
                                                        num_virtual_max,
                                                        d_keep_indices.data,
                                                        m_num_keep.getDeviceFlags());
            ScopedAllocation<unsigned char> d_tmp_alloc(this->m_exec_conf->getCachedAllocator(),
                                                        (tmp_storage_bytes > 0) ? tmp_storage_bytes
                                                                                : 1);
            d_tmp_storage = (void*)d_tmp_alloc();

            // run selection
            mpcd::gpu::compact_virtual_particle_indices(d_tmp_storage,
                                                        tmp_storage_bytes,
                                                        d_keep_particles.data,
                                                        num_virtual_max,
                                                        d_keep_indices.data,
                                                        m_num_keep.getDeviceFlags());
            }
        num_selected = m_num_keep.readFlags();
        }

    // Step 3: Allocate memory for the new virtual particles and copy them in. The tags can only be
    // assigned now because the number that survived rejection is not known in advance.
    const unsigned int first_tag = this->computeFirstTag(num_selected);
    const unsigned int first_idx = this->m_mpcd_pdata->addVirtualParticles(num_selected);
        {
        ArrayHandle<Scalar4> d_tmp_pos(this->m_tmp_pos, access_location::device, access_mode::read);
        ArrayHandle<Scalar4> d_tmp_vel(this->m_tmp_vel, access_location::device, access_mode::read);
        ArrayHandle<unsigned int> d_keep_indices(m_keep_indices,
                                                 access_location::device,
                                                 access_mode::read);
        ArrayHandle<Scalar4> d_pos(this->m_mpcd_pdata->getPositions(),
                                   access_location::device,
                                   access_mode::readwrite);
        ArrayHandle<Scalar4> d_vel(this->m_mpcd_pdata->getVelocities(),
                                   access_location::device,
                                   access_mode::readwrite);
        ArrayHandle<unsigned int> d_tag(this->m_mpcd_pdata->getTags(),
                                        access_location::device,
                                        access_mode::readwrite);
        m_tuner3->begin();
        mpcd::gpu::copy_virtual_particles(d_keep_indices.data,
                                          d_pos.data,
                                          d_vel.data,
                                          d_tag.data,
                                          d_tmp_pos.data,
                                          d_tmp_vel.data,
                                          first_idx,
                                          first_tag,
                                          num_selected,
                                          m_tuner3->getParam()[0]);
        if (this->m_exec_conf->isCUDAErrorCheckingEnabled())
            CHECK_CUDA_ERROR();
        m_tuner3->end();
        }
    }

namespace detail
    {
//! Export RejectionVirtualParticleFillerGPU to python
template<class Geometry> void export_RejectionVirtualParticleFillerGPU(pybind11::module& m)
    {
    namespace py = pybind11;
    const std::string name = Geometry::getName() + "GeometryFillerGPU";
    py::class_<mpcd::RejectionVirtualParticleFillerGPU<Geometry>,
               mpcd::RejectionVirtualParticleFiller<Geometry>,
               std::shared_ptr<mpcd::RejectionVirtualParticleFillerGPU<Geometry>>>(m, name.c_str())
        .def(pybind11::init<std::shared_ptr<SystemDefinition>,
                            const std::string&,
                            Scalar,
                            std::shared_ptr<Variant>,
                            std::shared_ptr<const Geometry>,
                            unsigned int,
                            unsigned int>());
    }
    } // end namespace detail
    } // end namespace mpcd
    } // namespace hoomd
#endif // MPCD_REJECTION_VIRTUAL_PARTICLE_FILLER_GPU_H_
