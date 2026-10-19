// Copyright (c) 2009-2026 The Regents of the University of Michigan.
// Part of HOOMD-blue, released under the BSD 3-Clause License.

/*!
 * \file mpcd/TriangulatedGeometryFillerGPU.cu
 * \brief Definitions of CUDA kernels for mpcd::TriangulatedGeometryFillerGPU
 */

#include "TriangulatedGeometryFillerGPU.cuh"

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wconversion"
#include <hipcub/hipcub.hpp>
#pragma GCC diagnostic pop
#include <thrust/iterator/counting_iterator.h>

#include <neighbor/neighbor.h>

namespace hoomd
    {
namespace mpcd
    {
namespace gpu
    {
#define DEVICE __device__ __forceinline__

//! Query operation for cells
/*!
 * One traversal thread handles one cell. The query bounding volume is the AABB of the cell grown
 * by \a args.max_shift on both sides, which is the swept region of the cell when the maximum grid
 * shift is passed and the cell itself when zero is passed. Any triangle that overlaps it is
 * forwarded to TriangleOutputOp::process as a candidate triangle.
 */
struct CellQueryOp
    {
    //! Constructor
    /*!
     * \param args_ Common arguments for classification
     * \param d_cells_ Local 1D indices of the cells to query, or NULL for all local cells
     * \param N_ Number of cells to query
     */
    CellQueryOp(const triangulated_classify_cells_args_t& args_,
                const unsigned int* d_cells_,
                unsigned int N_)
        : args(args_), d_cells(d_cells_), N(N_)
        {
        }

    //! Data stored per thread for traversal
    struct ThreadData
        {
        DEVICE ThreadData(unsigned int cell_idx_) : cell_idx(cell_idx_) { }

        unsigned int cell_idx; //!< Local 1D index of the cell
        };

    // Traversal volume type used by neighbor::LBVHTraverser
    typedef neighbor::BoundingBox Volume;

    //! Loads the per-thread data
    /*!
     * \param i Traversal thread index in [0, N)
     * \returns The ThreadData required for traversal
     */
    DEVICE ThreadData setup(unsigned int i) const
        {
        return ThreadData((d_cells) ? d_cells[i] : i);
        }

    //! Return the traversal volume
    /*!
     * \param q The current thread data
     * \returns The traversal bounding volume
     *
     * The search volume is the AABB of the cell grown by \a args.max_shift on both sides. The
     * image argument is ignored as triangles are already unwrapped.
     */
    DEVICE Volume get(const ThreadData& q, const float3&) const
        {
        const uint3 cell_ijk = args.ci.getTriple(q.cell_idx);
        const int gi = (int)cell_ijk.x + args.origin.x;
        const int gj = (int)cell_ijk.y + args.origin.y;
        const int gk = (int)cell_ijk.z + args.origin.z;

        // fractional bounds of the cell grown by the margin on both sides
        const Scalar3 f_lo = make_scalar3(gi * args.inv_dim.x - args.max_shift.x,
                                          gj * args.inv_dim.y - args.max_shift.y,
                                          gk * args.inv_dim.z - args.max_shift.z);
        const Scalar3 f_hi = make_scalar3((gi + 1) * args.inv_dim.x + args.max_shift.x,
                                          (gj + 1) * args.inv_dim.y + args.max_shift.y,
                                          (gk + 1) * args.inv_dim.z + args.max_shift.z);

        // take the bounding box of the 8 corners so that tilted boxes are handled
        Scalar3 lower = args.global_box.makeCoordinates(f_lo);
        Scalar3 upper = lower;
        for (int cx = 0; cx < 2; ++cx)
            {
            for (int cy = 0; cy < 2; ++cy)
                {
                for (int cz = 0; cz < 2; ++cz)
                    {
                    const Scalar3 f = make_scalar3(cx ? f_hi.x : f_lo.x,
                                                   cy ? f_hi.y : f_lo.y,
                                                   cz ? f_hi.z : f_lo.z);
                    const Scalar3 r = args.global_box.makeCoordinates(f);
                    lower.x = min(lower.x, r.x);
                    lower.y = min(lower.y, r.y);
                    lower.z = min(lower.z, r.z);
                    upper.x = max(upper.x, r.x);
                    upper.y = max(upper.y, r.y);
                    upper.z = max(upper.z, r.z);
                    }
                }
            }
        return neighbor::BoundingBox(lower, upper);
        }

    //! Perform the overlap test with the LBVH
    /*!
     * \param v Traversal volume
     * \param box Box in LBVH to intersect with
     * \returns True if the volume and box overlap
     *
     * This performs broad-phase search through the overlap test between the query box of the
     * cell and the triangle AABB.
     */
    DEVICE bool overlap(const Volume& v, const neighbor::BoundingBox& box) const
        {
        return v.overlap(box);
        }

    //! Refine the rough overlap test with a primitive
    /*!
     * \param q The current thread data
     * \param primitive Index of the intersected primitive
     * \returns True if the volumes still overlap after refinement
     *
     * Every overlapping triangle is a candidate, so this always returns \a True.
     */
    DEVICE bool refine(const ThreadData&, const int) const
        {
        return true;
        }

    //! Get the number of query threads
    __host__ DEVICE unsigned int size() const
        {
        return N;
        }

    const triangulated_classify_cells_args_t args; //!< Common arguments for classification
    const unsigned int* d_cells;                   //!< Cells to query, or NULL for all local cells
    const unsigned int N;                          //!< Number of cells to query
    };

//! Output operation for the triangles found for a cell
/*!
 * The triangles found for a cell are counted, and listed as well when \a d_tri_list is set.
 * The first traversal only counts so that the lists can be sized, and the later ones list the
 * candidate triangles.
 */
struct TriangleOutputOp
    {
    //! Constructor
    /*!
     * \param d_num_tri_per_cell_ Number of triangles per cell
     * \param d_tri_list_ Triangles, one block of \a max_tri_per_cell_ per cell,
     *                    or NULL to count only
     * \param max_tri_per_cell_ Uniform per-cell capacity of \a d_tri_list_
     * \param d_candidate_flags_ Flag set if the cell has any triangles, or NULL
     * \param d_max_tri_per_cell_ Largest count over the cells, or NULL
     */
    TriangleOutputOp(unsigned int* d_num_tri_per_cell_,
                     unsigned int* d_tri_list_,
                     unsigned int max_tri_per_cell_,
                     bool* d_candidate_flags_,
                     unsigned int* d_max_tri_per_cell_)
        : d_num_tri_per_cell(d_num_tri_per_cell_), d_tri_list(d_tri_list_),
          max_tri_per_cell(max_tri_per_cell_), d_candidate_flags(d_candidate_flags_),
          d_max_tri_per_cell(d_max_tri_per_cell_)
        {
        }

    //! Thread-local data
    struct ThreadData
        {
        DEVICE ThreadData(unsigned int i_) : i(i_), count(0) { }

        unsigned int i;     //!< Thread index
        unsigned int count; //!< Number of triangles found
        };

    //! Setup the thread data
    /*!
     * \param i Thread index
     * \returns The ThreadData for output
     *
     * \tparam Type of QueryData.
     */
    template<class QueryDataT> DEVICE ThreadData setup(unsigned int i, const QueryDataT&) const
        {
        return ThreadData(i);
        }

    //! Processes a newly intersected primitive
    /*!
     * \param t The thread's output data
     * \param primitive The index of the candidate triangle
     *
     * When listing, a counting traversal with a box at least as large already sized the blocks,
     * so the capacity is never exceeded.
     */
    DEVICE void process(ThreadData& t, const int primitive) const
        {
        if (d_tri_list && t.count < max_tri_per_cell)
            {
            d_tri_list[t.i * max_tri_per_cell + t.count] = (unsigned int)primitive;
            }
        ++t.count;
        }

    //! Finish the output job once the thread is ready to terminate
    /*!
     * \param t The thread's output data
     *
     * Cells that no triangles can reach under grid shifting are not candidates. The largest
     * count sizes the per-cell blocks of the triangle lists.
     */
    DEVICE void finalize(ThreadData& t) const
        {
        d_num_tri_per_cell[t.i] = (d_tri_list) ? min(t.count, max_tri_per_cell) : t.count;
        if (d_candidate_flags)
            {
            d_candidate_flags[t.i] = (t.count > 0);
            }
        if (d_max_tri_per_cell && t.count > 0)
            {
            atomicMax(d_max_tri_per_cell, t.count);
            }
        }

    unsigned int* d_num_tri_per_cell;    //!< Number of triangles per cell
    unsigned int* d_tri_list;            //!< Triangles per cell, or NULL to count only
    const unsigned int max_tri_per_cell; //!< Uniform per-cell capacity of d_tri_list
    bool* d_candidate_flags;             //!< Flag set if the cell has any triangles, or NULL
    unsigned int* d_max_tri_per_cell;    //!< Largest count over the cells, or NULL
    };

/*!
 * Initializes the underlying LBVHTraverser.
 */
CellLBVHTraverserWrapper::CellLBVHTraverserWrapper()
    {
    trav_ = new neighbor::LBVHTraverser();
    }

CellLBVHTraverserWrapper::~CellLBVHTraverserWrapper()
    {
    delete trav_;
    }

/*!
 * \param lbvh LBVH to traverse
 * \param stream CUDA stream for execution
 *
 * A NullTransformOp is used as the LBVH primitive order is the triangle order.
 */
void CellLBVHTraverserWrapper::setup(neighbor::LBVH& lbvh, hipStream_t stream)
    {
    neighbor::NullTransformOp transform;
    trav_->setup(stream, lbvh, transform);
    }

/*!
 * \param args Common arguments for classification, whose max_shift sets the query box
 * \param d_cells Local 1D indices of the cells to query, or NULL for all local cells
 * \param num_cells Number of cells to query
 * \param d_num_tri_per_cell Number of triangles per cell
 * \param d_tri_list Triangles, one block of \a max_tri_per_cell per cell, or NULL to count only
 * \param max_tri_per_cell Uniform per-cell capacity of \a d_tri_list
 * \param d_candidate_flags Flag set if the cell has any triangles, or NULL
 * \param d_max_tri_per_cell Largest count over the cells, or NULL
 * \param lbvh LBVH to traverse
 * \param stream CUDA stream for execution
 * \param block_size CUDA block size for execution
 *
 * Each cell is queried with the bounding box of the cell grown by \a args.max_shift. A single
 * zero image (neighbor::SelfOp) is used as the triangles are already unwrapped given the unwrap
 * distance. A NullTransformOp is used as no primitive mapping is needed.
 */
void CellLBVHTraverserWrapper::traverse(const triangulated_classify_cells_args_t& args,
                                        const unsigned int* d_cells,
                                        unsigned int num_cells,
                                        unsigned int* d_num_tri_per_cell,
                                        unsigned int* d_tri_list,
                                        unsigned int max_tri_per_cell,
                                        bool* d_candidate_flags,
                                        unsigned int* d_max_tri_per_cell,
                                        neighbor::LBVH& lbvh,
                                        hipStream_t stream,
                                        unsigned int block_size)
    {
    CellQueryOp query(args, d_cells, num_cells);
    TriangleOutputOp output(d_num_tri_per_cell,
                            d_tri_list,
                            max_tri_per_cell,
                            d_candidate_flags,
                            d_max_tri_per_cell);

    neighbor::NullTransformOp transform;
    neighbor::SelfOp translate;

    trav_->traverse(neighbor::LBVHTraverser::LaunchParameters(block_size, stream),
                    lbvh,
                    query,
                    output,
                    translate,
                    transform);
    }

std::vector<unsigned int> CellLBVHTraverserWrapper::getTunableParameters() const
    {
    return trav_->getTunableParameters();
    }

namespace kernel
    {
//! Kernel to flag the candidate cells that are cut by the surface
/*!
 * \param d_cut_flags Flags marking the candidate cells that are cut by the surface
 * \param d_center_outside Flags marking the candidate cells whose center is outside the fluid
 * \param d_tri_list Candidate triangles, one block of \a max_tri_per_cell per candidate cell
 * \param d_num_tri_per_cell Number of candidate triangles per candidate cell
 * \param max_tri_per_cell Uniform per-cell capacity of \a d_tri_list
 * \param d_candidate_cells Local 1D indices of the candidate cells
 * \param num_candidate_cells Number of candidate cells
 * \param d_vertices Triangle vertices
 * \param d_triangles Triangle vertex indices
 * \param global_box Global simulation box
 * \param draw_box Box with the extent of the swept region and the tilt of the global box
 * \param ci Indexer for the local cells
 * \param global_dim Global cell dimensions
 * \param origin Index of the first local cell in the global grid
 * \param inv_dim Fractional width of one cell along each lattice vector
 * \param num_classify_trials Number of trial points drawn per cell
 * \param seed User seed for RNG
 * \param filler_id Identifier for the filler (rng argument)
 *
 * \b implementation
 * We assign one thread per candidate cell. Trial points are drawn in the swept region of the
 * cell and tested against the candidate triangles found in that region, and the cell is cut as
 * soon as points are found on both sides of the surface. The center of a cut cell is tested as
 * well so that the cell can be decided later if no triangles reach into the cell itself.
 */
__global__ void triangulated_classify_cells(bool* d_cut_flags,
                                            bool* d_center_outside,
                                            const unsigned int* d_tri_list,
                                            const unsigned int* d_num_tri_per_cell,
                                            const unsigned int max_tri_per_cell,
                                            const unsigned int* d_candidate_cells,
                                            const unsigned int num_candidate_cells,
                                            const ShortReal3* d_vertices,
                                            const uint3* d_triangles,
                                            const BoxDim global_box,
                                            const BoxDim draw_box,
                                            const Index3D ci,
                                            const uint3 global_dim,
                                            const int3 origin,
                                            const Scalar3 inv_dim,
                                            const unsigned int num_classify_trials,
                                            const unsigned int seed,
                                            const unsigned int filler_id)
    {
    // one thread per candidate cell
    const unsigned int n = blockIdx.x * blockDim.x + threadIdx.x;
    if (n >= num_candidate_cells)
        return;

    const unsigned int cell_idx = d_candidate_cells[n];
    const uint3 cell_ijk = ci.getTriple(cell_idx);
    const int gi = (int)cell_ijk.x + origin.x;
    const int gj = (int)cell_ijk.y + origin.y;
    const int gk = (int)cell_ijk.z + origin.z;

    const Scalar3 f_center = make_scalar3((gi + Scalar(0.5)) * inv_dim.x,
                                          (gj + Scalar(0.5)) * inv_dim.y,
                                          (gk + Scalar(0.5)) * inv_dim.z);
    const Scalar3 cell_center = global_box.makeCoordinates(f_center);

    const Scalar3 half_w = make_scalar3(Scalar(0.5) * draw_box.getL().x,
                                        Scalar(0.5) * draw_box.getL().y,
                                        Scalar(0.5) * draw_box.getL().z);

    // candidate triangles in the swept region
    const unsigned int* cell_tris = d_tri_list + n * max_tri_per_cell;
    const unsigned int num_cell_tris = d_num_tri_per_cell[n];

    hoomd::RandomGenerator rng(
        hoomd::Seed(hoomd::RNGIdentifier::VirtualParticleFiller, 0, seed),
        hoomd::Counter(wrappedCellIndex(gi, gj, gk, global_dim), filler_id, 0));

    // 0 until the first usable trial point, then -1 if only inside points have been seen so far
    // and +1 if only outside points have been seen
    int code = 0;
    bool is_cut = false;
    for (unsigned int t = 0; t < num_classify_trials; ++t)
        {
        Scalar3 point = make_scalar3(hoomd::UniformDistribution<Scalar>(-half_w.x, half_w.x)(rng),
                                     hoomd::UniformDistribution<Scalar>(-half_w.y, half_w.y)(rng),
                                     hoomd::UniformDistribution<Scalar>(-half_w.z, half_w.z)(rng));

        int3 img = make_int3(0, 0, 0);
        draw_box.wrap(point, img);
        point += cell_center;

        // cells at the edge of the box have a swept region that overhangs the global box, so
        // trial points can land outside it. The geometry is assumed to lie inside the global
        // box, so these points are discarded.
        const Scalar3 f = global_box.makeFraction(point);
        if (f.x < Scalar(0.0) || f.x >= Scalar(1.0) || f.y < Scalar(0.0) || f.y >= Scalar(1.0)
            || f.z < Scalar(0.0) || f.z >= Scalar(1.0))
            continue;

        const bool is_outside = isOutside(point, cell_tris, num_cell_tris, d_vertices, d_triangles);

        if (code == 0)
            {
            code = is_outside ? 1 : -1;
            continue;
            }

        // one point on each side is enough to decide, so stop as soon as a trial point
        // contradicts the ones before it
        if ((code == -1 && is_outside) || (code == 1 && !is_outside))
            {
            is_cut = true;
            break;
            }
        }

    d_cut_flags[n] = is_cut;

    // the side of the cell is decided with the swept region candidates, which always hold the
    // nearest triangle to the center
    d_center_outside[n]
        = is_cut && isOutside(cell_center, cell_tris, num_cell_tris, d_vertices, d_triangles);
    }

//! Kernel to flag the cut cells that need filling
/*!
 * \param d_fill_flags Flags marking the cut cells that need filling
 * \param d_num_tri_per_cell Number of candidate triangles per cut cell
 * \param d_cut_candidates Index of each cut cell in the candidate arrays
 * \param d_center_outside Flags marking the candidate cells whose center is outside the fluid
 * \param num_cut_cells Number of cut cells
 *
 * \b implementation
 * We assign one thread per cut cell. A cut cell is filled if any triangles reach into it, or
 * without candidates if its center is outside the fluid, and dropped otherwise.
 */
__global__ void triangulated_flag_fill_cells(bool* d_fill_flags,
                                             const unsigned int* d_num_tri_per_cell,
                                             const unsigned int* d_cut_candidates,
                                             const bool* d_center_outside,
                                             const unsigned int num_cut_cells)
    {
    // one thread per cut cell
    const unsigned int n = blockIdx.x * blockDim.x + threadIdx.x;
    if (n >= num_cut_cells)
        return;

    d_fill_flags[n] = (d_num_tri_per_cell[n] > 0) || d_center_outside[d_cut_candidates[n]];
    }

//! Kernel to draw a Poisson number of virtual particles in each cell that needs filling
/*!
 * \param d_tmp_pos Temporary positions
 * \param d_tmp_vel Temporary velocities
 * \param d_keep_particles Flags for the particles that lie outside the fluid
 * \param d_fill_cells Local 1D indices of the cells that need filling
 * \param num_fill_cells Number of cells that need filling
 * \param global_box Global simulation box
 * \param draw_box Box with the extent of one cell and the tilt of the global box
 * \param ci Indexer for the local cells
 * \param global_dim Global cell dimensions
 * \param origin Index of the first local cell in the global grid
 * \param inv_dim Fractional width of one cell along each lattice vector
 * \param mean_per_cell Mean of the Poisson distribution for the number of particles per cell
 * \param max_per_cell Per-cell capacity of the temporary arrays
 * \param vel_factor Scale factor for uniform normal velocities consistent with particle mass /
 * temperature
 * \param type Particle type for filling
 * \param timestep Current timestep
 * \param seed User seed for RNG
 * \param filler_id Identifier for the filler (rng argument)
 * \param threads_per_cell Number of threads that draw for each cell
 * \param d_max_observed Largest per-cell Poisson draw that exceeded the capacity
 * \param d_fill_candidates Index of each fill cell in the cut cell arrays
 * \param d_tri_list Candidate triangles, one block of \a max_tri_per_cell per cut cell
 * \param d_num_tri_per_cell Number of candidate triangles per cut cell
 * \param max_tri_per_cell Uniform per-cell capacity of \a d_tri_list
 * \param d_vertices Triangle vertices
 * \param d_triangles Triangle vertex indices
 *
 * \b implementation
 * Same layout as mpcd::gpu::kernel::draw_virtual_particles: cells are processed by groups of
 * \a threads_per_cell threads, each cell owns a fixed block of \a max_per_cell entries, and a
 * cell whose draw exceeds the capacity records it in \a d_max_observed. The candidate triangles
 * of a cell are looked up in the cut cell arrays through \a d_fill_candidates. A particle is
 * kept if it lies outside the fluid according to those triangles, or without a test if the cell
 * has none because classification found it to be solid.
 */
__global__ void triangulated_draw_virtual_particles(Scalar4* d_tmp_pos,
                                                    Scalar4* d_tmp_vel,
                                                    bool* d_keep_particles,
                                                    const unsigned int* d_fill_cells,
                                                    const unsigned int num_fill_cells,
                                                    const BoxDim global_box,
                                                    const BoxDim draw_box,
                                                    const Index3D ci,
                                                    const uint3 global_dim,
                                                    const int3 origin,
                                                    const Scalar3 inv_dim,
                                                    const Scalar mean_per_cell,
                                                    const unsigned int max_per_cell,
                                                    const Scalar vel_factor,
                                                    const unsigned int type,
                                                    const uint64_t timestep,
                                                    const unsigned int seed,
                                                    const unsigned int filler_id,
                                                    const unsigned int threads_per_cell,
                                                    unsigned int* d_max_observed,
                                                    const unsigned int* d_fill_candidates,
                                                    const unsigned int* d_tri_list,
                                                    const unsigned int* d_num_tri_per_cell,
                                                    const unsigned int max_tri_per_cell,
                                                    const ShortReal3* d_vertices,
                                                    const uint3* d_triangles)
    {
    // multiple threads per cell
    const unsigned int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= num_fill_cells * threads_per_cell)
        return;
    const unsigned int n = idx / threads_per_cell;
    const unsigned int offset = idx % threads_per_cell;

    const unsigned int cell_idx = d_fill_cells[n];
    const uint3 cell_ijk = ci.getTriple(cell_idx);
    const int gi = (int)cell_ijk.x + origin.x;
    const int gj = (int)cell_ijk.y + origin.y;
    const int gk = (int)cell_ijk.z + origin.z;

    const Scalar3 f_center = make_scalar3((gi + Scalar(0.5)) * inv_dim.x,
                                          (gj + Scalar(0.5)) * inv_dim.y,
                                          (gk + Scalar(0.5)) * inv_dim.z);
    const Scalar3 cell_center = global_box.makeCoordinates(f_center);

    const Scalar3 half_w = make_scalar3(Scalar(0.5) * draw_box.getL().x,
                                        Scalar(0.5) * draw_box.getL().y,
                                        Scalar(0.5) * draw_box.getL().z);

    // candidate triangles stored for this cell during classification
    const unsigned int c = d_fill_candidates[n];
    const unsigned int* cell_tris = d_tri_list + c * max_tri_per_cell;
    const unsigned int num_cell_tris = d_num_tri_per_cell[c];

    hoomd::RandomGenerator rng_cell(
        hoomd::Seed(hoomd::RNGIdentifier::VirtualParticleFiller, timestep, seed),
        hoomd::Counter(wrappedCellIndex(gi, gj, gk, global_dim), filler_id, 1));

    unsigned int num_in_cell = hoomd::PoissonDistribution<Scalar>(mean_per_cell)(rng_cell);
    if (num_in_cell > max_per_cell)
        {
        // report the overflow for reallocation
        if (offset == 0)
            {
            atomicMax(d_max_observed, num_in_cell);
            }
        num_in_cell = 0;
        }

    hoomd::RandomGenerator rng(
        hoomd::Seed(hoomd::RNGIdentifier::VirtualParticleFiller, timestep, seed),
        hoomd::Counter(wrappedCellIndex(gi, gj, gk, global_dim), filler_id, 2 + offset));

    const unsigned int base = n * max_per_cell;
    for (unsigned int p = offset; p < max_per_cell; p += threads_per_cell)
        {
        if (p >= num_in_cell)
            {
            d_keep_particles[base + p] = false;
            continue;
            }
        Scalar3 pos = make_scalar3(hoomd::UniformDistribution<Scalar>(-half_w.x, half_w.x)(rng),
                                   hoomd::UniformDistribution<Scalar>(-half_w.y, half_w.y)(rng),
                                   hoomd::UniformDistribution<Scalar>(-half_w.z, half_w.z)(rng));

        int3 img = make_int3(0, 0, 0);
        draw_box.wrap(pos, img);
        pos += cell_center;

        // a cell with no candidate triangles was classified as solid, so every particle is kept
        const bool is_outside
            = (num_cell_tris == 0)
              || isOutside(pos, cell_tris, num_cell_tris, d_vertices, d_triangles);
        d_keep_particles[base + p] = is_outside;

        if (!is_outside)
            continue;

        d_tmp_pos[base + p] = make_scalar4(pos.x, pos.y, pos.z, __int_as_scalar(type));

        hoomd::NormalDistribution<Scalar> gen(vel_factor, 0.0);
        Scalar3 vel;
        gen(vel.x, vel.y, rng);
        vel.z = gen(rng);
        d_tmp_vel[base + p]
            = make_scalar4(vel.x, vel.y, vel.z, __int_as_scalar(mpcd::detail::NO_CELL));
        }
    }
    } // end namespace kernel

/*!
 * \param args Common arguments for classification
 * \param d_cut_flags Flags marking the candidate cells that are cut by the surface
 * \param d_center_outside Flags marking the candidate cells whose center is outside the fluid
 * \param d_tri_list Candidate triangles, one block of \a max_tri_per_cell per candidate cell
 * \param d_num_tri_per_cell Number of candidate triangles per candidate cell
 * \param max_tri_per_cell Uniform per-cell capacity of \a d_tri_list
 * \param d_candidate_cells Local 1D indices of the candidate cells
 * \param num_candidate_cells Number of candidate cells
 * \param d_vertices Triangle vertices
 * \param d_triangles Triangle vertex indices
 * \param block_size CUDA block size for execution
 *
 * \sa mpcd::gpu::kernel::triangulated_classify_cells
 */
hipError_t triangulated_classify_cells(const triangulated_classify_cells_args_t& args,
                                       bool* d_cut_flags,
                                       bool* d_center_outside,
                                       const unsigned int* d_tri_list,
                                       const unsigned int* d_num_tri_per_cell,
                                       unsigned int max_tri_per_cell,
                                       const unsigned int* d_candidate_cells,
                                       unsigned int num_candidate_cells,
                                       const ShortReal3* d_vertices,
                                       const uint3* d_triangles,
                                       unsigned int block_size)
    {
    if (num_candidate_cells == 0)
        return hipSuccess;

    hipFuncAttributes attr;
    hipFuncGetAttributes(&attr, reinterpret_cast<const void*>(kernel::triangulated_classify_cells));
    const unsigned int max_block_size = attr.maxThreadsPerBlock;

    const unsigned int run_block_size = min(block_size, max_block_size);
    hipLaunchKernelGGL(kernel::triangulated_classify_cells,
                       dim3(num_candidate_cells / run_block_size + 1),
                       dim3(run_block_size),
                       0,
                       0,
                       d_cut_flags,
                       d_center_outside,
                       d_tri_list,
                       d_num_tri_per_cell,
                       max_tri_per_cell,
                       d_candidate_cells,
                       num_candidate_cells,
                       d_vertices,
                       d_triangles,
                       args.global_box,
                       args.draw_box,
                       args.ci,
                       args.global_dim,
                       args.origin,
                       args.inv_dim,
                       args.num_classify_trials,
                       args.seed,
                       args.filler_id);
    return hipSuccess;
    }

/*!
 * \param d_fill_flags Flags marking the cut cells that need filling
 * \param d_num_tri_per_cell Number of candidate triangles per cut cell
 * \param d_cut_candidates Index of each cut cell in the candidate arrays
 * \param d_center_outside Flags marking the candidate cells whose center is outside the fluid
 * \param num_cut_cells Number of cut cells
 * \param block_size CUDA block size for execution
 *
 * \sa mpcd::gpu::kernel::triangulated_flag_fill_cells
 */
hipError_t triangulated_flag_fill_cells(bool* d_fill_flags,
                                        const unsigned int* d_num_tri_per_cell,
                                        const unsigned int* d_cut_candidates,
                                        const bool* d_center_outside,
                                        unsigned int num_cut_cells,
                                        unsigned int block_size)
    {
    if (num_cut_cells == 0)
        return hipSuccess;

    hipFuncAttributes attr;
    hipFuncGetAttributes(&attr,
                         reinterpret_cast<const void*>(kernel::triangulated_flag_fill_cells));
    const unsigned int max_block_size = attr.maxThreadsPerBlock;

    const unsigned int run_block_size = min(block_size, max_block_size);
    hipLaunchKernelGGL(kernel::triangulated_flag_fill_cells,
                       dim3(num_cut_cells / run_block_size + 1),
                       dim3(run_block_size),
                       0,
                       0,
                       d_fill_flags,
                       d_num_tri_per_cell,
                       d_cut_candidates,
                       d_center_outside,
                       num_cut_cells);
    return hipSuccess;
    }

/*!
 * \param d_tmp Temporary storage
 * \param tmp_bytes Number of bytes in temporary storage
 * \param d_flags Flags marking the cells to keep
 * \param d_cells Local 1D indices of the cells
 * \param num_cells Number of cells
 * \param d_out_cells Local 1D indices of the flagged cells
 * \param d_out_indices Index of each flagged cell in \a d_cells
 * \param d_num_out Number of flagged cells
 */
hipError_t compact_flagged_cells(void* d_tmp,
                                 size_t& tmp_bytes,
                                 const bool* d_flags,
                                 const unsigned int* d_cells,
                                 unsigned int num_cells,
                                 unsigned int* d_out_cells,
                                 unsigned int* d_out_indices,
                                 unsigned int* d_num_out)
    {
    thrust::counting_iterator<unsigned int> itr(0);
    if (d_tmp == NULL)
        {
        size_t cells_bytes = 0;
        size_t indices_bytes = 0;
        hipcub::DeviceSelect::Flagged(d_tmp,
                                      cells_bytes,
                                      d_cells,
                                      d_flags,
                                      d_out_cells,
                                      d_num_out,
                                      num_cells);
        hipcub::DeviceSelect::Flagged(d_tmp,
                                      indices_bytes,
                                      itr,
                                      d_flags,
                                      d_out_indices,
                                      d_num_out,
                                      num_cells);
        tmp_bytes = max(cells_bytes, indices_bytes);
        return hipSuccess;
        }

    hipcub::DeviceSelect::Flagged(d_tmp,
                                  tmp_bytes,
                                  d_cells,
                                  d_flags,
                                  d_out_cells,
                                  d_num_out,
                                  num_cells);
    hipcub::DeviceSelect::Flagged(d_tmp,
                                  tmp_bytes,
                                  itr,
                                  d_flags,
                                  d_out_indices,
                                  d_num_out,
                                  num_cells);
    return hipSuccess;
    }

/*!
 * \param args Common arguments for drawing virtual particles cell by cell
 * \param d_fill_candidates Index of each fill cell in the cut cell arrays
 * \param d_tri_list Candidate triangles, one block of \a max_tri_per_cell per cut cell
 * \param d_num_tri_per_cell Number of candidate triangles per cut cell
 * \param max_tri_per_cell Uniform per-cell capacity of \a d_tri_list
 * \param d_vertices Triangle vertices
 * \param d_triangles Triangle vertex indices
 *
 * \sa mpcd::gpu::kernel::triangulated_draw_virtual_particles
 */
hipError_t triangulated_draw_virtual_particles(const draw_virtual_particles_args_t& args,
                                               const unsigned int* d_fill_candidates,
                                               const unsigned int* d_tri_list,
                                               const unsigned int* d_num_tri_per_cell,
                                               unsigned int max_tri_per_cell,
                                               const ShortReal3* d_vertices,
                                               const uint3* d_triangles)
    {
    if (args.num_fill_cells == 0)
        return hipSuccess;

    hipFuncAttributes attr;
    hipFuncGetAttributes(
        &attr,
        reinterpret_cast<const void*>(kernel::triangulated_draw_virtual_particles));
    const unsigned int max_block_size = attr.maxThreadsPerBlock;

    const unsigned int run_block_size = min(args.block_size, max_block_size);
    const unsigned int num_threads = args.num_fill_cells * args.threads_per_cell;
    hipLaunchKernelGGL(kernel::triangulated_draw_virtual_particles,
                       dim3(num_threads / run_block_size + 1),
                       dim3(run_block_size),
                       0,
                       0,
                       args.d_tmp_pos,
                       args.d_tmp_vel,
                       args.d_keep_particles,
                       args.d_fill_cells,
                       args.num_fill_cells,
                       args.global_box,
                       args.draw_box,
                       args.ci,
                       args.global_dim,
                       args.origin,
                       args.inv_dim,
                       args.mean_per_cell,
                       args.max_per_cell,
                       args.vel_factor,
                       args.type,
                       args.timestep,
                       args.seed,
                       args.filler_id,
                       args.threads_per_cell,
                       args.d_max_observed,
                       d_fill_candidates,
                       d_tri_list,
                       d_num_tri_per_cell,
                       max_tri_per_cell,
                       d_vertices,
                       d_triangles);
    return hipSuccess;
    }

#undef DEVICE

    } // end namespace gpu
    } // end namespace mpcd
    } // end namespace hoomd
