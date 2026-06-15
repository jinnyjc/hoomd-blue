// Copyright (c) 2009-2025 The Regents of the University of Michigan.
// Part of HOOMD-blue, released under the BSD 3-Clause License.

/*!
 * \file mpcd/TriangulatedGeometryStreamingMethodGPU.h
 * \brief Declaration of mpcd::TriangulatedGeometryStreamingMethodGPU
 */

#ifndef MPCD_TRIANGULATED_GEOMETRY_STREAMING_METHOD_GPU_H_
#define MPCD_TRIANGULATED_GEOMETRY_STREAMING_METHOD_GPU_H_

#ifdef __HIPCC__
#error This header cannot be compiled by nvcc
#endif

#include "TriangulatedGeometryStreamingMethod.h"
#include "TriangulatedGeometryStreamingMethodGPU.cuh"
#include "hoomd/Autotuner.h"
#include "hoomd/CachedAllocator.h"

namespace hoomd
    {
namespace mpcd
    {
//! MPCD triangulated geometry streaming method
/*!
 * This method implements the GPU version of ballistic propagation of MPCD
 * particles in confined triangulated geometry.
 */
template<class Force>
class PYBIND11_EXPORT TriangulatedGeometryStreamingMethodGPU
    : public mpcd::TriangulatedGeometryStreamingMethod<Force>
    {
    public:
    //! Constructor
    /*!
     * \param sysdef System definition
     * \param cur_timestep Current system timestep
     * \param period Number of timesteps between collisions
     * \param phase Phase shift for periodic updates
     * \param geom Triangulated geometry
     * \param force Solvent force
     */
    TriangulatedGeometryStreamingMethodGPU(std::shared_ptr<SystemDefinition> sysdef,
                                           unsigned int cur_timestep,
                                           unsigned int period,
                                           int phase,
                                           std::shared_ptr<TriangulatedGeometry> geom,
                                           std::shared_ptr<Force> force,
                                           unsigned int max_bounce)
        : mpcd::TriangulatedGeometryStreamingMethod<Force>(sysdef,
                                                           cur_timestep,
                                                           period,
                                                           phase,
                                                           geom,
                                                           force,
                                                           max_bounce),
          m_lbvh_valid(false), m_num_active(this->m_exec_conf)
        {
        m_lbvh.reset(new gpu::TriangleLBVHWrapper());
        m_traverser.reset(new gpu::TriangleLBVHTraverserWrapper());

        m_build_tuner.reset(new Autotuner<1>({m_lbvh->getTunableParameters()},
                                             this->m_exec_conf,
                                             "mpcd_triangulated_lbvh_build"));

        m_init_tuner.reset(new Autotuner<1>({AutotunerBase::makeBlockSizeRange(this->m_exec_conf)},
                                            this->m_exec_conf,
                                            "mpcd_triangulated_stream_init"));

        m_traverse_tuner.reset(new Autotuner<1>({m_traverser->getTunableParameters()},
                                                this->m_exec_conf,
                                                "mpcd_triangulated_stream_traverse"));

        m_finalize_tuner.reset(
            new Autotuner<1>({AutotunerBase::makeBlockSizeRange(this->m_exec_conf)},
                             this->m_exec_conf,
                             "mpcd_triangulated_stream_finalize"));

        this->m_autotuners.insert(
            this->m_autotuners.end(),
            {m_build_tuner, m_init_tuner, m_traverse_tuner, m_finalize_tuner});
        }

    //! Implementation of the streaming rule
    void stream(uint64_t timestep) override;

    protected:
    // Autotuners
    std::shared_ptr<Autotuner<1>> m_build_tuner;    //!< Tuner for LBVH build
    std::shared_ptr<Autotuner<1>> m_init_tuner;     //!< Tuner for init kernel
    std::shared_ptr<Autotuner<1>> m_traverse_tuner; //!< Tuner for traversal
    std::shared_ptr<Autotuner<1>> m_finalize_tuner; //!< Tuner for finalize kernel

    // LBVH
    std::unique_ptr<gpu::TriangleLBVHWrapper> m_lbvh;               //!< LBVH
    std::unique_ptr<gpu::TriangleLBVHTraverserWrapper> m_traverser; //!< LBVH traverser
    bool m_lbvh_valid; //!< whether the LBVH has been built

    // Per-particle loop arrays
    GPUArray<Scalar> m_dt_remain;             //!< Remaining timesteps
    GPUArray<unsigned int> m_active_idx;      //!< Active particle indices
    GPUArray<unsigned int> m_next_active_idx; //!< Next active particle indices
    GPUArray<unsigned int> m_flags;           //!< Particle flag: 1 if dt_remain > 0
    GPUFlags<unsigned int> m_num_active;      //!< Number of active particles
    };

/*!
 * \param timestep Current time to stream
 *
 * This implements the three-phase GPU streaming algorithm:
 *  1. Init kernel: apply first half-kick from the solvent force and initialize
 *        dt_remain and active list
 *  2. Traverse loop: for every active particle whose \a dt_remain > 0, query the LBVH to
 *        find the nearest triangle hit along its remaining path.
 *         - If a hit is found, advance the particle to the collision point, subtract the
 *           time consumed, and reflect the velocity based on boundary conditions.
 *         - Otherwise, the particle is advanced by \a dt_remain * \a vel and \a dt_remain
 *           becomes zero.
 *        Particles with  \a dt_remain > 0 are then compacted into the next active list
 *        with CUB::DevicePartition::Flagged, and the loop repeats until every particle has
 *        consumed its full timestep.
 *  3. Finalize kernel: apply final half-kick and wrap the particles back into the
 *        simulation box.
 */
template<class Force> void TriangulatedGeometryStreamingMethodGPU<Force>::stream(uint64_t timestep)
    {
    if (!this->shouldStream(timestep))
        return;

    if (!this->m_cl)
        {
        throw std::runtime_error("Cell list has not been set");
        }

    if (!this->m_geom)
        throw std::runtime_error("Triangulated geometry has not been set");

    const unsigned int N = this->m_mpcd_pdata->getN();

    // Allocate per-particle arrays
    if (m_dt_remain.getNumElements() != N)
        {
        GPUArray<Scalar> dt_remain(N, this->m_exec_conf);
        m_dt_remain.swap(dt_remain);

        GPUArray<unsigned int> active_idx(N, this->m_exec_conf);
        m_active_idx.swap(active_idx);

        GPUArray<unsigned int> next_active_idx(N, this->m_exec_conf);
        m_next_active_idx.swap(next_active_idx);

        GPUArray<unsigned int> flags(N, this->m_exec_conf);
        m_flags.swap(flags);
        }

    ArrayHandle<Scalar4> d_pos(this->m_mpcd_pdata->getPositions(),
                               access_location::device,
                               access_mode::readwrite);
    ArrayHandle<Scalar4> d_vel(this->m_mpcd_pdata->getVelocities(),
                               access_location::device,
                               access_mode::readwrite);
    ArrayHandle<Scalar3> d_vertices(this->m_geom->getVertices(),
                                    access_location::device,
                                    access_mode::read);
    ArrayHandle<uint3> d_triangles(this->m_geom->getTriangles(),
                                   access_location::device,
                                   access_mode::read);
    ArrayHandle<Scalar> d_dt_remain(m_dt_remain, access_location::device, access_mode::readwrite);
    ArrayHandle<unsigned int> d_active_idx(m_active_idx,
                                           access_location::device,
                                           access_mode::readwrite);
    ArrayHandle<unsigned int> d_next_active_idx(m_next_active_idx,
                                                access_location::device,
                                                access_mode::readwrite);
    ArrayHandle<unsigned int> d_flags(m_flags, access_location::device, access_mode::readwrite);

    // default construct a force if one is not set
    const Force force = (this->m_force) ? *(this->m_force) : Force();

    const Scalar mass = this->m_mpcd_pdata->getMass();
    const BoxDim box = this->m_cl->getCoverageBox();
    const Scalar unwrap_distance = this->m_geom->getUnwrapDistance();
    const Scalar3 lo
        = box.getLo() - make_scalar3(unwrap_distance, unwrap_distance, unwrap_distance);
    const Scalar3 hi
        = box.getHi() + make_scalar3(unwrap_distance, unwrap_distance, unwrap_distance);

    // build LBVH once; geometry is static
    if (!m_lbvh_valid)
        {
        m_lbvh->setup(d_vertices.data, d_triangles.data, this->m_geom->getNumTotalTriangles(), 0);

        m_build_tuner->begin();
        m_lbvh->build(d_vertices.data,
                      d_triangles.data,
                      this->m_geom->getNumTotalTriangles(),
                      lo,
                      hi,
                      0,
                      m_build_tuner->getParam()[0]);
        m_build_tuner->end();

        m_traverser->setup(*(m_lbvh->get()), 0);

        if (this->m_exec_conf->isCUDAErrorCheckingEnabled())
            CHECK_CUDA_ERROR();

        m_lbvh_valid = true;
        }

    mpcd::gpu::triangulated_stream_args_t args(d_pos.data,
                                               d_vel.data,
                                               mass,
                                               box,
                                               this->m_mpcd_dt,
                                               N);

    // Step 1: init
    m_init_tuner->begin();
    mpcd::gpu::triangulated_stream_init(args,
                                        force,
                                        d_dt_remain.data,
                                        d_active_idx.data,
                                        m_init_tuner->getParam()[0]);
    if (this->m_exec_conf->isCUDAErrorCheckingEnabled())
        CHECK_CUDA_ERROR();
    m_init_tuner->end();

    // Step 2: traverse loop
    unsigned int* cur_active = d_active_idx.data;
    unsigned int* next_active = d_next_active_idx.data;
    unsigned int h_num_active = N;
    unsigned int bounce = 0;

    do
        {
        m_traverse_tuner->begin();
        m_traverser->traverse(args,
                              d_vertices.data,
                              d_triangles.data,
                              d_dt_remain.data,
                              cur_active,
                              d_flags.data,
                              h_num_active,
                              this->m_geom->getNoSlip(),
                              *(m_lbvh->get()),
                              0,
                              m_traverse_tuner->getParam()[0]);
        if (this->m_exec_conf->isCUDAErrorCheckingEnabled())
            CHECK_CUDA_ERROR();
        m_traverse_tuner->end();

            {
            void* d_tmp = NULL;
            size_t tmp_bytes = 0;
            mpcd::gpu::triangulated_compact_active(d_tmp,
                                                   tmp_bytes,
                                                   cur_active,
                                                   next_active,
                                                   d_flags.data,
                                                   m_num_active.getDeviceFlags(),
                                                   h_num_active);

            // make requested temporary allocation
            ScopedAllocation<unsigned char> d_tmp_alloc(this->m_exec_conf->getCachedAllocator(),
                                                        (tmp_bytes > 0) ? tmp_bytes : 1);
            d_tmp = (void*)d_tmp_alloc();

            // partition particles to keep
            mpcd::gpu::triangulated_compact_active(d_tmp,
                                                   tmp_bytes,
                                                   cur_active,
                                                   next_active,
                                                   d_flags.data,
                                                   m_num_active.getDeviceFlags(),
                                                   h_num_active);
            }
        h_num_active = m_num_active.readFlags();

        std::swap(cur_active, next_active);
        ++bounce;
        } while (h_num_active > 0 && bounce < this->m_max_bounce);

    if (h_num_active > 0)
        {
        throw std::runtime_error("Particle did not finish collision after "
                                 + std::to_string(this->m_max_bounce)
                                 + " bounces. Check the triangulated geometry");
        }

    // Step 3: finalize
    m_finalize_tuner->begin();
    mpcd::gpu::triangulated_stream_finalize(args, force, m_finalize_tuner->getParam()[0]);
    if (this->m_exec_conf->isCUDAErrorCheckingEnabled())
        CHECK_CUDA_ERROR();
    m_finalize_tuner->end();

    // particles have moved, so the cell cache is no longer valid
    this->m_mpcd_pdata->invalidateCellCache();
    }

namespace detail
    {
//! Export mpcd::TriangulatedGeometryStreamingMethodGPU to python
/*!
 * \param m Python module to export to
 */
template<class Force> void export_TriangulatedGeometryStreamingMethodGPU(pybind11::module& m)
    {
    const std::string name = "TriangulatedGeometryStreamingMethod" + Force::getName() + "GPU";
    pybind11::class_<mpcd::TriangulatedGeometryStreamingMethodGPU<Force>,
                     mpcd::TriangulatedGeometryStreamingMethod<Force>,
                     std::shared_ptr<mpcd::TriangulatedGeometryStreamingMethodGPU<Force>>>(
        m,
        name.c_str())
        .def(pybind11::init<std::shared_ptr<SystemDefinition>,
                            unsigned int,
                            unsigned int,
                            int,
                            std::shared_ptr<TriangulatedGeometry>,
                            std::shared_ptr<Force>,
                            unsigned int>());
    }
    } // end namespace detail
    } // end namespace mpcd
    } // end namespace hoomd
#endif // MPCD_TRIANGULATED_GEOMETRY_STREAMING_METHOD_GPU_H_
