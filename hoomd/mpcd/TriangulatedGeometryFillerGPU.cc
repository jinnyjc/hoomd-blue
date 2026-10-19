// Copyright (c) 2009-2026 The Regents of the University of Michigan.
// Part of HOOMD-blue, released under the BSD 3-Clause License.

/*!
 * \file mpcd/TriangulatedGeometryFillerGPU.cc
 * \brief Definition of mpcd::TriangulatedGeometryFillerGPU
 */

#include "TriangulatedGeometryFillerGPU.h"

namespace hoomd
    {
namespace mpcd
    {
TriangulatedGeometryFillerGPU::TriangulatedGeometryFillerGPU(
    std::shared_ptr<SystemDefinition> sysdef,
    const std::string& type,
    Scalar density,
    std::shared_ptr<Variant> T,
    std::shared_ptr<const TriangulatedGeometry> geom,
    unsigned int num_classify_trials)
    : mpcd::TriangulatedGeometryFiller(sysdef, type, density, T, geom, num_classify_trials),
      m_fill_candidates(m_exec_conf), m_keep_particles(m_exec_conf), m_keep_indices(m_exec_conf),
      m_num_keep(m_exec_conf), m_max_observed(m_exec_conf)
    {
    m_traverser.reset(new gpu::CellLBVHTraverserWrapper());

    m_traverse_tuner.reset(new Autotuner<1>({m_traverser->getTunableParameters()},
                                            m_exec_conf,
                                            "mpcd_triangulated_filler_traverse"));
    m_classify_tuner.reset(new Autotuner<1>({AutotunerBase::makeBlockSizeRange(m_exec_conf)},
                                            m_exec_conf,
                                            "mpcd_triangulated_filler_classify_cells"));
    m_draw_tuner.reset(new Autotuner<2>({AutotunerBase::makeBlockSizeRange(m_exec_conf),
                                         AutotunerBase::getTppListPow2(m_exec_conf)},
                                        m_exec_conf,
                                        "mpcd_triangulated_filler_draw_particles"));
    m_tag_tuner.reset(new Autotuner<1>({AutotunerBase::makeBlockSizeRange(m_exec_conf)},
                                       m_exec_conf,
                                       "mpcd_triangulated_filler_tag_particles"));
    m_autotuners.insert(m_autotuners.end(),
                        {m_traverse_tuner, m_classify_tuner, m_draw_tuner, m_tag_tuner});
    }

void TriangulatedGeometryFillerGPU::classifyCells()
    {
    // size the cell list against the current box before reading its dimensions
    m_cl->computeDimensions();

    const BoxDim& global_box = m_pdata->getGlobalBox();
    const Scalar3 global_L = global_box.getL();

    // the cells this rank owns, and where they sit in the global grid
    const uint3 global_dim = m_cl->getGlobalDim();
    const int3 origin = m_cl->getOriginIndex();
    const Index3D& ci = m_cl->getCellIndexer();
    const unsigned int num_cells = ci.getNumElements();

    // fractional width of one cell along each lattice vector
    const Scalar3 inv_dim = make_scalar3(Scalar(1.0) / global_dim.x,
                                         Scalar(1.0) / global_dim.y,
                                         Scalar(1.0) / global_dim.z);

    // trial points are drawn in the swept region of each cell, see TriangulatedGeometryFiller
    const Scalar3 max_shift
        = m_cl->isGridShifting() ? m_cl->getMaxGridShift() : make_scalar3(0, 0, 0);
    const Scalar3 width = make_scalar3(inv_dim.x + Scalar(2.0) * max_shift.x,
                                       inv_dim.y + Scalar(2.0) * max_shift.y,
                                       inv_dim.z + Scalar(2.0) * max_shift.z);
    BoxDim draw_box(make_scalar3(width.x * global_L.x, width.y * global_L.y, width.z * global_L.z));
    draw_box.setTiltFactors(global_box.getTiltFactorXY(),
                            global_box.getTiltFactorXZ(),
                            global_box.getTiltFactorYZ());

    // when not set by the user, use the mean number of solvent particles in the swept region
    // plus eight standard deviations of the Poisson distribution
    const Scalar mean_in_region = m_density * draw_box.getVolume();
    const unsigned int num_classify_trials
        = (m_num_classify_trials > 0)
              ? m_num_classify_trials
              : static_cast<unsigned int>(
                    std::ceil(mean_in_region + Scalar(8.0) * std::sqrt(mean_in_region)));

    // the swept region is queried to classify the cells, and the cell itself is queried for the
    // candidate triangles that are kept for filling
    mpcd::gpu::triangulated_classify_cells_args_t swept_args(global_box,
                                                             draw_box,
                                                             ci,
                                                             global_dim,
                                                             origin,
                                                             inv_dim,
                                                             max_shift,
                                                             num_classify_trials,
                                                             m_sysdef->getSeed(),
                                                             m_filler_id);
    mpcd::gpu::triangulated_classify_cells_args_t cell_args(global_box,
                                                            draw_box,
                                                            ci,
                                                            global_dim,
                                                            origin,
                                                            inv_dim,
                                                            make_scalar3(0, 0, 0),
                                                            num_classify_trials,
                                                            m_sysdef->getSeed(),
                                                            m_filler_id);

    ArrayHandle<ShortReal3> d_vertices(m_geom->getVertices(),
                                       access_location::device,
                                       access_mode::read);
    ArrayHandle<uint3> d_triangles(m_geom->getTriangles(),
                                   access_location::device,
                                   access_mode::read);

    // Step 1: Set up the traverser
    neighbor::LBVH& lbvh = *(m_geom->getTriangleLBVH().get());
    m_traverser->setup(lbvh, 0);
    if (m_exec_conf->isCUDAErrorCheckingEnabled())
        CHECK_CUDA_ERROR();

    // Step 2: Count the triangles in the swept region of every cell. Cells that no triangles
    // can reach under grid shifting are not candidates, and the largest count sizes the
    // per-cell blocks of the candidate triangle lists.
    GPUArray<unsigned int> cell_num_tri(num_cells, m_exec_conf);
    GPUArray<bool> candidate_flags(num_cells, m_exec_conf);
    GPUArray<unsigned int> candidate_cells(num_cells, m_exec_conf);
    unsigned int num_candidate_cells = 0;
    unsigned int max_tri_per_cell = 0;
        {
        ArrayHandle<unsigned int> d_cell_num_tri(cell_num_tri,
                                                 access_location::device,
                                                 access_mode::overwrite);
        ArrayHandle<bool> d_candidate_flags(candidate_flags,
                                            access_location::device,
                                            access_mode::overwrite);
        ArrayHandle<unsigned int> d_candidate_cells(candidate_cells,
                                                    access_location::device,
                                                    access_mode::overwrite);

        m_max_observed.resetFlags(0);
        m_traverse_tuner->begin();
        m_traverser->traverse(swept_args,
                              NULL,
                              num_cells,
                              d_cell_num_tri.data,
                              NULL,
                              0,
                              d_candidate_flags.data,
                              m_max_observed.getDeviceFlags(),
                              lbvh,
                              0,
                              m_traverse_tuner->getParam()[0]);
        if (m_exec_conf->isCUDAErrorCheckingEnabled())
            CHECK_CUDA_ERROR();
        m_traverse_tuner->end();
        max_tri_per_cell = m_max_observed.readFlags();

            {
            // compact the candidate cells into a list of indices with CUB
            void* d_tmp_storage = NULL;
            size_t tmp_storage_bytes = 0;
            gpu::compact_virtual_particle_indices(d_tmp_storage,
                                                  tmp_storage_bytes,
                                                  d_candidate_flags.data,
                                                  num_cells,
                                                  d_candidate_cells.data,
                                                  m_num_keep.getDeviceFlags());
            ScopedAllocation<unsigned char> d_tmp_alloc(m_exec_conf->getCachedAllocator(),
                                                        (tmp_storage_bytes > 0) ? tmp_storage_bytes
                                                                                : 1);
            d_tmp_storage = (void*)d_tmp_alloc();

            // run selection
            gpu::compact_virtual_particle_indices(d_tmp_storage,
                                                  tmp_storage_bytes,
                                                  d_candidate_flags.data,
                                                  num_cells,
                                                  d_candidate_cells.data,
                                                  m_num_keep.getDeviceFlags());
            }
        num_candidate_cells = m_num_keep.readFlags();
        }

    // Step 3: List the candidate triangles in the swept region of the candidate cells and
    // classify them with trial points. The cut cells are compacted and queried again with the
    // cell itself, which overwrites the lists in cut cell order, and a cut cell is a fill cell
    // if any triangles reach into it or its center is outside the fluid.
    m_max_tri_per_cell = max_tri_per_cell;
    if (num_candidate_cells > m_num_tri_per_cell.getNumElements())
        {
        GPUArray<unsigned int> num_tri_per_cell(num_candidate_cells, m_exec_conf);
        m_num_tri_per_cell.swap(num_tri_per_cell);
        GPUArray<unsigned int> fill_cells(num_candidate_cells, m_exec_conf);
        m_fill_cells.swap(fill_cells);
        GPUArray<unsigned int> fill_candidates(num_candidate_cells, m_exec_conf);
        m_fill_candidates.swap(fill_candidates);
        }
    const unsigned int num_tri_max = num_candidate_cells * m_max_tri_per_cell;
    if (num_tri_max > m_tri_list.getNumElements())
        {
        GPUArray<unsigned int> tri_list(num_tri_max, m_exec_conf);
        m_tri_list.swap(tri_list);
        }
    GPUArray<bool> cut_flags(num_candidate_cells, m_exec_conf);
    GPUArray<bool> center_outside(num_candidate_cells, m_exec_conf);
    GPUArray<unsigned int> cut_cells(num_candidate_cells, m_exec_conf);
    GPUArray<unsigned int> cut_candidates(num_candidate_cells, m_exec_conf);
    GPUArray<bool> fill_flags(num_candidate_cells, m_exec_conf);
    m_num_fill_cells = 0;
    if (num_candidate_cells > 0)
        {
        ArrayHandle<unsigned int> d_candidate_cells(candidate_cells,
                                                    access_location::device,
                                                    access_mode::read);
        ArrayHandle<unsigned int> d_num_tri_per_cell(m_num_tri_per_cell,
                                                     access_location::device,
                                                     access_mode::overwrite);
        ArrayHandle<unsigned int> d_tri_list(m_tri_list,
                                             access_location::device,
                                             access_mode::overwrite);
        ArrayHandle<bool> d_cut_flags(cut_flags, access_location::device, access_mode::overwrite);
        ArrayHandle<bool> d_center_outside(center_outside,
                                           access_location::device,
                                           access_mode::overwrite);
        ArrayHandle<unsigned int> d_cut_cells(cut_cells,
                                              access_location::device,
                                              access_mode::overwrite);
        ArrayHandle<unsigned int> d_cut_candidates(cut_candidates,
                                                   access_location::device,
                                                   access_mode::overwrite);
        ArrayHandle<bool> d_fill_flags(fill_flags, access_location::device, access_mode::overwrite);
        ArrayHandle<unsigned int> d_fill_cells(m_fill_cells,
                                               access_location::device,
                                               access_mode::overwrite);
        ArrayHandle<unsigned int> d_fill_candidates(m_fill_candidates,
                                                    access_location::device,
                                                    access_mode::overwrite);

        // candidate triangles in the swept region
        m_traverse_tuner->begin();
        m_traverser->traverse(swept_args,
                              d_candidate_cells.data,
                              num_candidate_cells,
                              d_num_tri_per_cell.data,
                              d_tri_list.data,
                              m_max_tri_per_cell,
                              NULL,
                              NULL,
                              lbvh,
                              0,
                              m_traverse_tuner->getParam()[0]);
        if (m_exec_conf->isCUDAErrorCheckingEnabled())
            CHECK_CUDA_ERROR();
        m_traverse_tuner->end();

        // trial points decide which cells are cut
        m_classify_tuner->begin();
        gpu::triangulated_classify_cells(swept_args,
                                         d_cut_flags.data,
                                         d_center_outside.data,
                                         d_tri_list.data,
                                         d_num_tri_per_cell.data,
                                         m_max_tri_per_cell,
                                         d_candidate_cells.data,
                                         num_candidate_cells,
                                         d_vertices.data,
                                         d_triangles.data,
                                         m_classify_tuner->getParam()[0]);
        if (m_exec_conf->isCUDAErrorCheckingEnabled())
            CHECK_CUDA_ERROR();
        m_classify_tuner->end();

        unsigned int num_cut_cells = 0;
            {
            // compact the cut cells and their index in the candidate arrays with CUB
            void* d_tmp_storage = NULL;
            size_t tmp_storage_bytes = 0;
            gpu::compact_flagged_cells(d_tmp_storage,
                                       tmp_storage_bytes,
                                       d_cut_flags.data,
                                       d_candidate_cells.data,
                                       num_candidate_cells,
                                       d_cut_cells.data,
                                       d_cut_candidates.data,
                                       m_num_keep.getDeviceFlags());
            ScopedAllocation<unsigned char> d_tmp_alloc(m_exec_conf->getCachedAllocator(),
                                                        (tmp_storage_bytes > 0) ? tmp_storage_bytes
                                                                                : 1);
            d_tmp_storage = (void*)d_tmp_alloc();

            // run selection
            gpu::compact_flagged_cells(d_tmp_storage,
                                       tmp_storage_bytes,
                                       d_cut_flags.data,
                                       d_candidate_cells.data,
                                       num_candidate_cells,
                                       d_cut_cells.data,
                                       d_cut_candidates.data,
                                       m_num_keep.getDeviceFlags());
            num_cut_cells = m_num_keep.readFlags();
            }

        if (num_cut_cells > 0)
            {
            // candidate triangles in the cell itself, overwriting the lists in cut cell order
            m_traverse_tuner->begin();
            m_traverser->traverse(cell_args,
                                  d_cut_cells.data,
                                  num_cut_cells,
                                  d_num_tri_per_cell.data,
                                  d_tri_list.data,
                                  m_max_tri_per_cell,
                                  NULL,
                                  NULL,
                                  lbvh,
                                  0,
                                  m_traverse_tuner->getParam()[0]);
            if (m_exec_conf->isCUDAErrorCheckingEnabled())
                CHECK_CUDA_ERROR();
            m_traverse_tuner->end();

            // decide the fill cells
            m_classify_tuner->begin();
            gpu::triangulated_flag_fill_cells(d_fill_flags.data,
                                              d_num_tri_per_cell.data,
                                              d_cut_candidates.data,
                                              d_center_outside.data,
                                              num_cut_cells,
                                              m_classify_tuner->getParam()[0]);
            if (m_exec_conf->isCUDAErrorCheckingEnabled())
                CHECK_CUDA_ERROR();
            m_classify_tuner->end();

                {
                // compact the fill cells and their index in the cut cell arrays with CUB
                void* d_tmp_storage = NULL;
                size_t tmp_storage_bytes = 0;
                gpu::compact_flagged_cells(d_tmp_storage,
                                           tmp_storage_bytes,
                                           d_fill_flags.data,
                                           d_cut_cells.data,
                                           num_cut_cells,
                                           d_fill_cells.data,
                                           d_fill_candidates.data,
                                           m_num_keep.getDeviceFlags());
                ScopedAllocation<unsigned char> d_tmp_alloc(
                    m_exec_conf->getCachedAllocator(),
                    (tmp_storage_bytes > 0) ? tmp_storage_bytes : 1);
                d_tmp_storage = (void*)d_tmp_alloc();

                // run selection
                gpu::compact_flagged_cells(d_tmp_storage,
                                           tmp_storage_bytes,
                                           d_fill_flags.data,
                                           d_cut_cells.data,
                                           num_cut_cells,
                                           d_fill_cells.data,
                                           d_fill_candidates.data,
                                           m_num_keep.getDeviceFlags());
                m_num_fill_cells = m_num_keep.readFlags();
                }
            }
        }

    m_exec_conf->msg->notice(6) << "MPCD TriangulatedGeometryFillerGPU: filling "
                                << m_num_fill_cells << " of " << num_cells << " cells, at most "
                                << m_max_tri_per_cell << " triangles per cell" << std::endl;
    }

void TriangulatedGeometryFillerGPU::fill(uint64_t timestep)
    {
    // size the cell list against the current box before reading its dimensions
    m_cl->computeDimensions();

    if (m_need_classify)
        {
        classifyCells();
        m_need_classify = false;
        }

    const BoxDim& global_box = m_pdata->getGlobalBox();
    const Scalar3 global_L = global_box.getL();

    // the cells this rank owns, and where they sit in the global grid
    const uint3 global_dim = m_cl->getGlobalDim();
    const int3 origin = m_cl->getOriginIndex();
    const Index3D& ci = m_cl->getCellIndexer();

    // fractional width of one cell along each lattice vector
    const Scalar3 inv_dim = make_scalar3(Scalar(1.0) / global_dim.x,
                                         Scalar(1.0) / global_dim.y,
                                         Scalar(1.0) / global_dim.z);

    // particles are drawn in the cell itself rather than in the swept region, so this box carries
    // the cell width and the tilt of the global box. It is centered on the origin, and a particle
    // is placed by translating onto the center of its cell
    BoxDim draw_box(
        make_scalar3(inv_dim.x * global_L.x, inv_dim.y * global_L.y, inv_dim.z * global_L.z));
    draw_box.setTiltFactors(global_box.getTiltFactorXY(),
                            global_box.getTiltFactorXZ(),
                            global_box.getTiltFactorYZ());

    // the number of particles drawn per cell is Poisson distributed with mean equal to the fill
    // density times the cell volume, which is the same for every cell even in a triclinic box
    // because shear preserves volume
    const Scalar cell_volume
        = global_box.getVolume() / (global_dim.x * global_dim.y * global_dim.z);
    const Scalar mean_per_cell = m_density * cell_volume;

    // initial guess for the per-cell capacity of the temporary arrays: eight standard deviations
    // above the Poisson mean
    if (m_alloc_per_cell == 0)
        {
        m_alloc_per_cell = static_cast<unsigned int>(
            std::ceil(mean_per_cell + Scalar(8.0) * std::sqrt(mean_per_cell)));
        }

    const Scalar vel_factor = fast::sqrt((*m_T)(timestep) / m_mpcd_pdata->getMass());

    unsigned int num_selected = 0;
    unsigned int max_observed = m_alloc_per_cell;
    do
        {
        // grow the capacity if the previous attempt overflowed
        m_alloc_per_cell = max_observed;
        m_max_observed.resetFlags(0);

        // Step 1: Size the temporary arrays for the current per-cell capacity.
        const unsigned int num_virtual_max = m_num_fill_cells * m_alloc_per_cell;
        if (num_virtual_max > m_tmp_pos.getNumElements())
            {
            GPUArray<Scalar4> tmp_pos(num_virtual_max, m_exec_conf);
            m_tmp_pos.swap(tmp_pos);
            GPUArray<Scalar4> tmp_vel(num_virtual_max, m_exec_conf);
            m_tmp_vel.swap(tmp_vel);
            GPUArray<bool> keep_particles(num_virtual_max, m_exec_conf);
            m_keep_particles.swap(keep_particles);
            GPUArray<unsigned int> keep_indices(num_virtual_max, m_exec_conf);
            m_keep_indices.swap(keep_indices);
            }

            // Step 2: Draw a Poisson number of particles in each cell that needs filling, then keep
            // only the ones that are outside the fluid
            {
            ArrayHandle<Scalar4> d_tmp_pos(m_tmp_pos,
                                           access_location::device,
                                           access_mode::overwrite);
            ArrayHandle<Scalar4> d_tmp_vel(m_tmp_vel,
                                           access_location::device,
                                           access_mode::overwrite);
            ArrayHandle<bool> d_keep_particles(m_keep_particles,
                                               access_location::device,
                                               access_mode::overwrite);
            ArrayHandle<unsigned int> d_keep_indices(m_keep_indices,
                                                     access_location::device,
                                                     access_mode::overwrite);
            ArrayHandle<unsigned int> d_fill_cells(m_fill_cells,
                                                   access_location::device,
                                                   access_mode::read);
            ArrayHandle<unsigned int> d_fill_candidates(m_fill_candidates,
                                                        access_location::device,
                                                        access_mode::read);
            ArrayHandle<unsigned int> d_tri_list(m_tri_list,
                                                 access_location::device,
                                                 access_mode::read);
            ArrayHandle<unsigned int> d_num_tri_per_cell(m_num_tri_per_cell,
                                                         access_location::device,
                                                         access_mode::read);
            ArrayHandle<ShortReal3> d_vertices(m_geom->getVertices(),
                                               access_location::device,
                                               access_mode::read);
            ArrayHandle<uint3> d_triangles(m_geom->getTriangles(),
                                           access_location::device,
                                           access_mode::read);

            mpcd::gpu::draw_virtual_particles_args_t args(d_tmp_pos.data,
                                                          d_tmp_vel.data,
                                                          d_keep_particles.data,
                                                          d_fill_cells.data,
                                                          m_num_fill_cells,
                                                          global_box,
                                                          draw_box,
                                                          ci,
                                                          global_dim,
                                                          origin,
                                                          inv_dim,
                                                          mean_per_cell,
                                                          m_alloc_per_cell,
                                                          vel_factor,
                                                          m_type,
                                                          timestep,
                                                          m_sysdef->getSeed(),
                                                          m_filler_id,
                                                          m_draw_tuner->getParam()[0],
                                                          m_draw_tuner->getParam()[1],
                                                          m_max_observed.getDeviceFlags());
            m_draw_tuner->begin();
            mpcd::gpu::triangulated_draw_virtual_particles(args,
                                                           d_fill_candidates.data,
                                                           d_tri_list.data,
                                                           d_num_tri_per_cell.data,
                                                           m_max_tri_per_cell,
                                                           d_vertices.data,
                                                           d_triangles.data);
            if (m_exec_conf->isCUDAErrorCheckingEnabled())
                CHECK_CUDA_ERROR();
            m_draw_tuner->end();

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
                ScopedAllocation<unsigned char> d_tmp_alloc(
                    m_exec_conf->getCachedAllocator(),
                    (tmp_storage_bytes > 0) ? tmp_storage_bytes : 1);
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

        // the kernel only records draws that exceeded the capacity, so zero means no overflow
        max_observed = m_max_observed.readFlags();
        } while (max_observed > m_alloc_per_cell);

    // Step 3: Allocate memory for the new virtual particles and copy them in. The tags can only be
    // assigned now because the number that survived rejection is not known in advance.
    const unsigned int first_tag = computeFirstTag(num_selected);
    const unsigned int first_idx = m_mpcd_pdata->addVirtualParticles(num_selected);
        {
        ArrayHandle<Scalar4> d_tmp_pos(m_tmp_pos, access_location::device, access_mode::read);
        ArrayHandle<Scalar4> d_tmp_vel(m_tmp_vel, access_location::device, access_mode::read);
        ArrayHandle<unsigned int> d_keep_indices(m_keep_indices,
                                                 access_location::device,
                                                 access_mode::read);
        ArrayHandle<Scalar4> d_pos(m_mpcd_pdata->getPositions(),
                                   access_location::device,
                                   access_mode::readwrite);
        ArrayHandle<Scalar4> d_vel(m_mpcd_pdata->getVelocities(),
                                   access_location::device,
                                   access_mode::readwrite);
        ArrayHandle<unsigned int> d_tag(m_mpcd_pdata->getTags(),
                                        access_location::device,
                                        access_mode::readwrite);
        m_tag_tuner->begin();
        mpcd::gpu::copy_virtual_particles(d_keep_indices.data,
                                          d_pos.data,
                                          d_vel.data,
                                          d_tag.data,
                                          d_tmp_pos.data,
                                          d_tmp_vel.data,
                                          first_idx,
                                          first_tag,
                                          num_selected,
                                          m_tag_tuner->getParam()[0]);
        if (m_exec_conf->isCUDAErrorCheckingEnabled())
            CHECK_CUDA_ERROR();
        m_tag_tuner->end();
        }
    }

namespace detail
    {
//! Export TriangulatedGeometryFillerGPU to python
void export_TriangulatedGeometryFillerGPU(pybind11::module& m)
    {
    pybind11::class_<mpcd::TriangulatedGeometryFillerGPU,
                     mpcd::TriangulatedGeometryFiller,
                     std::shared_ptr<mpcd::TriangulatedGeometryFillerGPU>>(
        m,
        "TriangulatedGeometryFillerGPU")
        .def(pybind11::init<std::shared_ptr<SystemDefinition>,
                            const std::string&,
                            Scalar,
                            std::shared_ptr<Variant>,
                            std::shared_ptr<const TriangulatedGeometry>,
                            unsigned int>());
    }
    } // end namespace detail
    } // end namespace mpcd
    } // end namespace hoomd
