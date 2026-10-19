// Copyright (c) 2009-2026 The Regents of the University of Michigan.
// Part of HOOMD-blue, released under the BSD 3-Clause License.

/*!
 * \file mpcd/RejectionVirtualParticleFiller.h
 * \brief Declaration and definition of RejectionVirtualParticleFiller
 */

#ifndef MPCD_REJECTION_VIRTUAL_PARTICLE_FILLER_H_
#define MPCD_REJECTION_VIRTUAL_PARTICLE_FILLER_H_

#ifdef __HIPCC__
#error This header cannot be compiled by nvcc
#endif

#include "VirtualParticleFiller.h"
#include "WrappedCellIndex.h"

#include "hoomd/RNGIdentifiers.h"
#include "hoomd/RandomNumbers.h"

#include <pybind11/pybind11.h>

namespace hoomd
    {
namespace mpcd
    {

//! Adds virtual particles to MPCD particle data for a given geometry.
/*!
 * Here we implement a virtual particle filler using rejection sampling method. The filler first
 * identifies the collision cells that are cut by the boundary, i.e., the cells that can contain
 * both fluid and solid, and then only draws particles in those cells.
 *
 * To identify these cells, a number of trial points are drawn per cell in a region that covers
 * every position the cell can occupy under grid shifting. If all trial points lie inside the
 * geometry, the cell is pure fluid and needs no virtual particles. If all lie outside, the cell
 * is entirely solid and can never share a collision cell with fluid particles. Only the cells
 * with trial points on both sides of the boundary are filled. This classification is done once
 * and recomputed only when the simulation box changes.
 *
 * During filling, the number of particles drawn in each cell is Poisson distributed with mean equal
 * to the fill density times the cell volume. Each drawn particle is kept only if it lies outside
 * the confinement defined by the template geometry.
 */
template<class Geometry>
class PYBIND11_EXPORT RejectionVirtualParticleFiller : public mpcd::VirtualParticleFiller
    {
    public:
    //! Constructor
    RejectionVirtualParticleFiller(std::shared_ptr<SystemDefinition> sysdef,
                                   const std::string& type,
                                   Scalar density,
                                   std::shared_ptr<Variant> T,
                                   std::shared_ptr<const Geometry> geom,
                                   unsigned int num_classify_trials)
        : mpcd::VirtualParticleFiller(sysdef, type, density, T), m_geom(geom),
          m_tmp_pos(m_exec_conf), m_tmp_vel(m_exec_conf),
          m_num_classify_trials(num_classify_trials), m_fill_cells(m_exec_conf),
          m_num_fill_cells(0), m_need_classify(true), m_alloc_per_cell(0)
        {
        m_exec_conf->msg->notice(5)
            << "Constructing MPCD RejectionVirtualParticleFiller : " + Geometry::getName()
            << std::endl;

        m_pdata->getBoxChangeSignal()
            .connect<mpcd::RejectionVirtualParticleFiller<Geometry>,
                     &mpcd::RejectionVirtualParticleFiller<Geometry>::setNeedClassify>(this);
        }

    //! Destructor
    virtual ~RejectionVirtualParticleFiller()
        {
        m_exec_conf->msg->notice(5)
            << "Destroying MPCD RejectionVirtualParticleFiller" << std::endl;

        m_pdata->getBoxChangeSignal()
            .disconnect<mpcd::RejectionVirtualParticleFiller<Geometry>,
                        &mpcd::RejectionVirtualParticleFiller<Geometry>::setNeedClassify>(this);
        }

    //! Get the streaming geometry
    std::shared_ptr<const Geometry> getGeometry() const
        {
        return m_geom;
        }

    //! Set the streaming geometry
    void setGeometry(std::shared_ptr<const Geometry> geom)
        {
        m_geom = geom;
        m_need_classify = true;
        }

    //! Get the number of trial points used per cell during classification
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
    std::shared_ptr<const Geometry> m_geom;
    GPUArray<Scalar4> m_tmp_pos;
    GPUArray<Scalar4> m_tmp_vel;

    unsigned int m_num_classify_trials;  //!< Trial points per cell during classification
    GPUArray<unsigned int> m_fill_cells; //!< Local 1D indices of the cells that need filling
    unsigned int m_num_fill_cells;       //!< Number of cells that need filling
    bool m_need_classify;                //!< True if the cell classification is out of date
    unsigned int m_alloc_per_cell;       //!< Per-cell capacity of the temporary arrays

    //! Flag the cell classification as out of date
    void setNeedClassify()
        {
        m_need_classify = true;
        }

    //! Determine which cells need virtual particles
    virtual void classifyCells();
    };

template<class Geometry> void RejectionVirtualParticleFiller<Geometry>::classifyCells()
    {
    // size the cell list against the current box before reading its dimensions
    m_cl->computeDimensions();

    const BoxDim& global_box = m_pdata->getGlobalBox();
    const Scalar3 global_L = global_box.getL();

    // the cells this rank owns, and where they sit in the global grid
    const uint3 local_dim = m_cl->getDim();
    const uint3 global_dim = m_cl->getGlobalDim();
    const int3 origin = m_cl->getOriginIndex();
    const Index3D& ci = m_cl->getCellIndexer();

    // fractional width of one cell along each lattice vector
    const Scalar3 inv_dim = make_scalar3(Scalar(1.0) / global_dim.x,
                                         Scalar(1.0) / global_dim.y,
                                         Scalar(1.0) / global_dim.z);

    // a cell can be shifted by at most this much along each lattice vector, so sampling the cell
    // grown by that amount on both sides covers every position it can take under grid shifting
    const Scalar3 max_shift
        = m_cl->isGridShifting() ? m_cl->getMaxGridShift() : make_scalar3(0, 0, 0);
    const Scalar3 width = make_scalar3(inv_dim.x + Scalar(2.0) * max_shift.x,
                                       inv_dim.y + Scalar(2.0) * max_shift.y,
                                       inv_dim.z + Scalar(2.0) * max_shift.z);

    // trial points are drawn uniformly in an orthorhombic box with the size of the sampling
    // region, wrapped by this tilted box onto the skewed shape a cell actually occupies, and
    // then translated onto each cell by adding its center
    BoxDim draw_box(make_scalar3(width.x * global_L.x, width.y * global_L.y, width.z * global_L.z));
    draw_box.setTiltFactors(global_box.getTiltFactorXY(),
                            global_box.getTiltFactorXZ(),
                            global_box.getTiltFactorYZ());
    const Scalar3 half_w = make_scalar3(Scalar(0.5) * draw_box.getL().x,
                                        Scalar(0.5) * draw_box.getL().y,
                                        Scalar(0.5) * draw_box.getL().z);

    const Scalar mean_per_cell = m_density * draw_box.getVolume();

    // when not set by the user, use the mean number of solvent particles in the sampling region
    // plus eight standard deviations of the Poisson distribution
    const unsigned int num_classify_trials
        = (m_num_classify_trials > 0)
              ? m_num_classify_trials
              : static_cast<unsigned int>(
                    std::ceil(mean_per_cell + Scalar(8.0) * std::sqrt(mean_per_cell)));

    uint16_t seed = m_sysdef->getSeed();

    std::vector<unsigned int> marked;
    for (unsigned int k = 0; k < local_dim.z; ++k)
        {
        for (unsigned int j = 0; j < local_dim.y; ++j)
            {
            for (unsigned int i = 0; i < local_dim.x; ++i)
                {
                const int gi = static_cast<int>(i) + origin.x;
                const int gj = static_cast<int>(j) + origin.y;
                const int gk = static_cast<int>(k) + origin.z;

                // cell center in Cartesian coordinates
                const Scalar3 f_center = make_scalar3((gi + Scalar(0.5)) * inv_dim.x,
                                                      (gj + Scalar(0.5)) * inv_dim.y,
                                                      (gk + Scalar(0.5)) * inv_dim.z);
                const Scalar3 cell_center = global_box.makeCoordinates(f_center);

                hoomd::RandomGenerator rng(
                    hoomd::Seed(hoomd::RNGIdentifier::VirtualParticleFiller, 0, seed),
                    hoomd::Counter(wrappedCellIndex(gi, gj, gk, global_dim), m_filler_id, 0));

                // 0 until the first usable trial point, then -1 if only inside points have been
                // seen so far and +1 if only outside points have been seen
                int code = 0;
                for (unsigned int n = 0; n < num_classify_trials; ++n)
                    {
                    Scalar3 point = make_scalar3(
                        hoomd::UniformDistribution<Scalar>(-half_w.x, half_w.x)(rng),
                        hoomd::UniformDistribution<Scalar>(-half_w.y, half_w.y)(rng),
                        hoomd::UniformDistribution<Scalar>(-half_w.z, half_w.z)(rng));

                    int3 img = make_int3(0, 0, 0);
                    draw_box.wrap(point, img);
                    point += cell_center;

                    // cells at the edge of the box have a sampling region that overhangs the
                    // global box, so trial points can land outside it. The geometry is assumed
                    // to lie inside the global box, so these points are discarded.
                    const Scalar3 f = global_box.makeFraction(point);
                    if (f.x < Scalar(0.0) || f.x >= Scalar(1.0) || f.y < Scalar(0.0)
                        || f.y >= Scalar(1.0) || f.z < Scalar(0.0) || f.z >= Scalar(1.0))
                        {
                        continue;
                        }

                    const bool is_outside = m_geom->isOutside(point);

                    if (code == 0)
                        {
                        code = is_outside ? 1 : -1;
                        continue;
                        }

                    // one point on each side is enough to decide, so stop as soon as a trial point
                    // contradicts the ones before it
                    if ((code == -1 && is_outside) || (code == 1 && !is_outside))
                        {
                        marked.push_back(ci(i, j, k));
                        break;
                        }
                    }
                }
            }
        }

    m_num_fill_cells = static_cast<unsigned int>(marked.size());
    if (m_num_fill_cells > m_fill_cells.getNumElements())
        {
        GPUArray<unsigned int> fill_cells(m_num_fill_cells, m_exec_conf);
        m_fill_cells.swap(fill_cells);
        }

    if (m_num_fill_cells > 0)
        {
        ArrayHandle<unsigned int> h_fill_cells(m_fill_cells,
                                               access_location::host,
                                               access_mode::overwrite);
        std::copy(marked.begin(), marked.end(), h_fill_cells.data);
        }

    m_exec_conf->msg->notice(6) << "MPCD RejectionVirtualParticleFiller: filling "
                                << m_num_fill_cells << " of "
                                << (local_dim.x * local_dim.y * local_dim.z) << " cells"
                                << std::endl;
    }

template<class Geometry> void RejectionVirtualParticleFiller<Geometry>::fill(uint64_t timestep)
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

    // particles are drawn in the cell itself rather than in the region swept out by grid shifting,
    // so this box carries the cell width and the tilt of the global box. It is centered on the
    // origin, and a particle is placed by translating onto the center of its cell
    BoxDim draw_box(
        make_scalar3(inv_dim.x * global_L.x, inv_dim.y * global_L.y, inv_dim.z * global_L.z));
    draw_box.setTiltFactors(global_box.getTiltFactorXY(),
                            global_box.getTiltFactorXZ(),
                            global_box.getTiltFactorYZ());
    const Scalar3 half_w = make_scalar3(Scalar(0.5) * draw_box.getL().x,
                                        Scalar(0.5) * draw_box.getL().y,
                                        Scalar(0.5) * draw_box.getL().z);

    // the number of particles drawn per cell is Poisson distributed with mean equal to the fill
    // density times the cell volume, which is the same for every cell even in a triclinic box
    // because shear preserves volume.
    const Scalar cell_volume
        = global_box.getVolume() / (global_dim.x * global_dim.y * global_dim.z);
    const Scalar mean_per_cell = m_density * cell_volume;

    // initial guess for the per-cell capacity of the temporary arrays: eight standard deviations
    // above the Poisson mean.
    if (m_alloc_per_cell == 0)
        {
        m_alloc_per_cell = static_cast<unsigned int>(
            std::ceil(mean_per_cell + Scalar(8.0) * std::sqrt(mean_per_cell)));
        }

    unsigned int num_selected = 0;
    unsigned int max_observed = m_alloc_per_cell;
    do
        {
        // grow the capacity if the previous attempt overflowed
        m_alloc_per_cell = max_observed;

        // Step 1: Size the temporary arrays for the current per-cell capacity.
        const unsigned int num_virtual_max = m_num_fill_cells * m_alloc_per_cell;
        if (num_virtual_max > m_tmp_pos.getNumElements())
            {
            GPUArray<Scalar4> tmp_pos(num_virtual_max, m_exec_conf);
            GPUArray<Scalar4> tmp_vel(num_virtual_max, m_exec_conf);
            m_tmp_pos.swap(tmp_pos);
            m_tmp_vel.swap(tmp_vel);
            }

        // Step 2: Draw the particles and assign velocities simultaneously by using temporary
        // memory. Only keep the ones that are outside the geometry.
        uint16_t seed = m_sysdef->getSeed();
        const Scalar vel_factor = fast::sqrt((*m_T)(timestep) / m_mpcd_pdata->getMass());
        ArrayHandle<Scalar4> h_tmp_pos(m_tmp_pos, access_location::host, access_mode::overwrite);
        ArrayHandle<Scalar4> h_tmp_vel(m_tmp_vel, access_location::host, access_mode::overwrite);
        ArrayHandle<unsigned int> h_fill_cells(m_fill_cells,
                                               access_location::host,
                                               access_mode::read);

        num_selected = 0;
        for (unsigned int n = 0; n < m_num_fill_cells; ++n)
            {
            const unsigned int cell_idx = h_fill_cells.data[n];
            const uint3 cell_ijk = ci.getTriple(cell_idx);
            const int gi = static_cast<int>(cell_ijk.x) + origin.x;
            const int gj = static_cast<int>(cell_ijk.y) + origin.y;
            const int gk = static_cast<int>(cell_ijk.z) + origin.z;

            // cell center in Cartesian coordinates
            Scalar3 f_center = make_scalar3((gi + Scalar(0.5)) * inv_dim.x,
                                            (gj + Scalar(0.5)) * inv_dim.y,
                                            (gk + Scalar(0.5)) * inv_dim.z);
            const Scalar3 cell_center = global_box.makeCoordinates(f_center);

            hoomd::RandomGenerator rng(
                hoomd::Seed(hoomd::RNGIdentifier::VirtualParticleFiller, timestep, seed),
                hoomd::Counter(wrappedCellIndex(gi, gj, gk, global_dim), m_filler_id, 1));

            unsigned int num_in_cell = hoomd::PoissonDistribution<Scalar>(mean_per_cell)(rng);
            max_observed = std::max(max_observed, num_in_cell);

            // once any cell has drawn more than the capacity, the whole fill is redone with
            // larger arrays, so the rest of this attempt only tracks the largest draw.
            if (max_observed > m_alloc_per_cell)
                continue;

            for (unsigned int p = 0; p < num_in_cell; ++p)
                {
                Scalar3 particle
                    = make_scalar3(hoomd::UniformDistribution<Scalar>(-half_w.x, half_w.x)(rng),
                                   hoomd::UniformDistribution<Scalar>(-half_w.y, half_w.y)(rng),
                                   hoomd::UniformDistribution<Scalar>(-half_w.z, half_w.z)(rng));
                int3 img = make_int3(0, 0, 0);
                draw_box.wrap(particle, img);
                particle += cell_center;

                if (m_geom->isOutside(particle))
                    {
                    h_tmp_pos.data[num_selected]
                        = make_scalar4(particle.x, particle.y, particle.z, __int_as_scalar(m_type));

                    hoomd::NormalDistribution<Scalar> gen(vel_factor, 0.0);
                    Scalar3 vel;
                    gen(vel.x, vel.y, rng);
                    vel.z = gen(rng);
                    m_geom->addToVirtualParticleVelocity(vel, particle);
                    h_tmp_vel.data[num_selected]
                        = make_scalar4(vel.x, vel.y, vel.z, __int_as_scalar(mpcd::detail::NO_CELL));
                    ++num_selected;
                    }
                }
            }
        } while (max_observed > m_alloc_per_cell);

    // Step 3: Allocate memory for the new virtual particles and copy them in. The tags can only be
    // assigned now because the number that survived rejection is not known in advance.
    const unsigned int first_tag = computeFirstTag(num_selected);
    const unsigned int first_idx = m_mpcd_pdata->addVirtualParticles(num_selected);
    ArrayHandle<Scalar4> h_tmp_pos(m_tmp_pos, access_location::host, access_mode::read);
    ArrayHandle<Scalar4> h_tmp_vel(m_tmp_vel, access_location::host, access_mode::read);
    ArrayHandle<Scalar4> h_pos(m_mpcd_pdata->getPositions(),
                               access_location::host,
                               access_mode::readwrite);
    ArrayHandle<Scalar4> h_vel(m_mpcd_pdata->getVelocities(),
                               access_location::host,
                               access_mode::readwrite);
    ArrayHandle<unsigned int> h_tag(m_mpcd_pdata->getTags(),
                                    access_location::host,
                                    access_mode::readwrite);
    for (unsigned int i = 0; i < num_selected; ++i)
        {
        const unsigned int idx = first_idx + i;
        h_pos.data[idx] = h_tmp_pos.data[i];
        h_vel.data[idx] = h_tmp_vel.data[i];
        h_tag.data[idx] = first_tag + i;
        }
    }

namespace detail
    {
//! Export RejectionVirtualParticleFiller to python
template<class Geometry> void export_RejectionVirtualParticleFiller(pybind11::module& m)
    {
    namespace py = pybind11;
    const std::string name = Geometry::getName() + "GeometryFiller";
    py::class_<mpcd::RejectionVirtualParticleFiller<Geometry>,
               mpcd::VirtualParticleFiller,
               std::shared_ptr<mpcd::RejectionVirtualParticleFiller<Geometry>>>(m, name.c_str())
        .def(pybind11::init<std::shared_ptr<SystemDefinition>,
                            const std::string&,
                            Scalar,
                            std::shared_ptr<Variant>,
                            std::shared_ptr<const Geometry>,
                            unsigned int>())
        .def_property_readonly("geometry",
                               &mpcd::RejectionVirtualParticleFiller<Geometry>::getGeometry)
        .def_property_readonly(
            "num_classify_trials",
            &mpcd::RejectionVirtualParticleFiller<Geometry>::getNumClassifyTrials);
    }
    } // end namespace detail
    } // end namespace mpcd
    } // end namespace hoomd
#endif // MPCD_REJECTION_VIRTUAL_PARTICLE_FILLER_H_
