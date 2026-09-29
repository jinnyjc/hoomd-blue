// Copyright (c) 2009-2026 The Regents of the University of Michigan.
// Part of HOOMD-blue, released under the BSD 3-Clause License.

/*!
 * \file mpcd/TriangulatedGeometryFiller.h
 * \brief Declaration of TriangulatedGeometryFiller
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

//! Closest point on a triangle to a point
/*!
 * \param p Query point
 * \param a First triangle vertex
 * \param b Second triangle vertex
 * \param c Third triangle vertex
 * \returns The point on triangle abc (including its edges and vertices) closest to p
 *
 * The plane containing the triangle is partitioned into seven Voronoi regions (three
 * vertex regions, three edge regions, and the face). The region containing p is found
 * from the signs of six dot products, and p is projected accordingly: onto a vertex,
 * onto an edge segment, or onto the face using barycentric coordinates.
 *
 * This implementation follows the method described in:
 * Ericson, C. (2005). Real-Time Collision Detection. Section 5.1.5.
 */
inline Scalar3
closestPointOnTriangle(const Scalar3& p, const Scalar3& a, const Scalar3& b, const Scalar3& c)
    {
    const Scalar3 ab = b - a;
    const Scalar3 ac = c - a;
    const Scalar3 ap = p - a;

    // vertex region a
    const Scalar d1 = dot(ab, ap);
    const Scalar d2 = dot(ac, ap);
    if (d1 <= Scalar(0.0) && d2 <= Scalar(0.0))
        return a;

    // vertex region b
    const Scalar3 bp = p - b;
    const Scalar d3 = dot(ab, bp);
    const Scalar d4 = dot(ac, bp);
    if (d3 >= Scalar(0.0) && d4 <= d3)
        return b;

    // vertex region c
    const Scalar3 cp = p - c;
    const Scalar d5 = dot(ab, cp);
    const Scalar d6 = dot(ac, cp);
    if (d6 >= Scalar(0.0) && d5 <= d6)
        return c;

    // edge region ab
    const Scalar vc = d1 * d4 - d3 * d2;
    if (vc <= Scalar(0.0) && d1 >= Scalar(0.0) && d3 <= Scalar(0.0))
        {
        const Scalar v = d1 / (d1 - d3);
        return a + v * ab;
        }

    // edge region ac
    const Scalar vb = d5 * d2 - d1 * d6;
    if (vb <= Scalar(0.0) && d2 >= Scalar(0.0) && d6 <= Scalar(0.0))
        {
        const Scalar w = d2 / (d2 - d6);
        return a + w * ac;
        }

    // edge region bc
    const Scalar va = d3 * d6 - d5 * d4;
    if (va <= Scalar(0.0) && (d4 - d3) >= Scalar(0.0) && (d5 - d6) >= Scalar(0.0))
        {
        const Scalar w = (d4 - d3) / ((d4 - d3) + (d5 - d6));
        return b + w * (c - b);
        }

    // face region: interpolate with normalized barycentric coordinates
    const Scalar denom = Scalar(1.0) / (va + vb + vc);
    const Scalar v = vb * denom;
    const Scalar w = vc * denom;
    return a + v * ab + w * ac;
    }

//! Test whether a point lies outside the fluid using a list of candidate triangles
/*!
 * \param p Query point
 * \param tris Candidate triangle indices
 * \param num_tris Number of candidates
 * \param verts List of vertices
 * \param tri_data List of triangles
 * \returns True if p lies on the side that the normal of the nearest candidate points to
 *
 * The sign of the dot product between (p - q) and the triangle normal decides the side,
 * where q is the closest point on the nearest candidate. Normals point out of the fluid,
 * and a point exactly on the surface is treated as fluid.
 */
inline bool isOutside(const Scalar3& p,
                      const unsigned int* tris,
                      unsigned int num_tris,
                      const ShortReal3* verts,
                      const uint3* tri_data)
    {
    Scalar best_d2(0);
    Scalar best_sign(0);
    for (unsigned int t = 0; t < num_tris; ++t)
        {
        const uint3 tri = tri_data[tris[t]];

        const ShortReal3 va = verts[tri.x];
        const ShortReal3 vb = verts[tri.y];
        const ShortReal3 vc = verts[tri.z];
        const Scalar3 a = make_scalar3(va.x, va.y, va.z);
        const Scalar3 b = make_scalar3(vb.x, vb.y, vb.z);
        const Scalar3 c = make_scalar3(vc.x, vc.y, vc.z);

        const Scalar3 q = closestPointOnTriangle(p, a, b, c);
        const Scalar3 d = p - q;
        const Scalar d2 = dot(d, d);
        if (t == 0 || d2 < best_d2)
            {
            best_d2 = d2;
            const Scalar3 n = cross(b - a, c - a);
            best_sign = dot(d, n);
            }
        }
    return best_sign > Scalar(0.0);
    }

//! Adds virtual particles to MPCD particle data for a triangulated geometry
/*!
 * The filler first identifies the collision cells that are cut by the mesh boundary and
 * only draws particles in those cells. Because the triangles are static, the candidate
 * triangles for each cell are found once during classification and stored, then reused
 * every fill step. Classification builds the candidate list of each cut cell and stores
 * it in a uniform per-cell allocation sized by the largest count.
 *
 * A point is classified against the mesh by finding the nearest candidate triangle and
 * checking the sign of the dot product between (p - q) and the triangle normal, where q
 * is the closest point on that triangle. Normals point out of the fluid, so a strictly
 * positive sign means the point is outside the fluid and should be filled.
 *
 * During filling, a Poisson number of particles with mean equal to the fill density
 * times the cell volume is drawn uniformly in each marked cell. Each particle is tested
 * against the candidate triangles of its cell and kept only if it lies outside. The
 * per-cell capacity of temporary arrays is grown whenever a draw exceeds it.
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

    //! Get the streaming geometry
    std::shared_ptr<const TriangulatedGeometry> getGeometry() const
        {
        return m_geom;
        }

    //! Set the streaming geometry
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

    GPUArray<unsigned int> m_tri_list;         //!< Candidate triangles
    GPUArray<unsigned int> m_num_tri_per_cell; //!< Number of candidate triangles per fill cell
    unsigned int m_max_tri_per_cell;           //!< Uniform per-cell capacity of m_tri_list

    //! Flag the cell classification as out of date
    void setNeedClassify()
        {
        m_need_classify = true;
        }

    //! Determine which cells need virtual particles and store their candidate triangles
    virtual void classifyCells();

    //! Compute the Cartesian AABB of the region a cell sweeps under grid shifting
    /*!
     * \param gi Global cell index along x
     * \param gj Global cell index along y
     * \param gk Global cell index along z
     * \param inv_dim Fractional width of one cell along each lattice vector
     * \param max_shift Maximum fractional grid shift
     * \param global_box Global simulation box
     * \returns Axis-aligned bounding box of the swept region
     *
     * The swept region is the cell grown by the maximum grid shift on both sides. Its
     * bounding box is taken over the 8 corners so that tilted boxes are handled.
     */
    hoomd::detail::AABB computeSweptCellAABB(int gi,
                                             int gj,
                                             int gk,
                                             const Scalar3& inv_dim,
                                             const Scalar3& max_shift,
                                             const BoxDim& global_box) const
        {
        const Scalar3 f_lo = make_scalar3(gi * inv_dim.x - max_shift.x,
                                          gj * inv_dim.y - max_shift.y,
                                          gk * inv_dim.z - max_shift.z);
        const Scalar3 f_hi = make_scalar3((gi + 1) * inv_dim.x + max_shift.x,
                                          (gj + 1) * inv_dim.y + max_shift.y,
                                          (gk + 1) * inv_dim.z + max_shift.z);

        // take the bounding box of the 8 corners so that tilted boxes are handled
        Scalar3 lower = global_box.makeCoordinates(f_lo);
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
                    const Scalar3 r = global_box.makeCoordinates(f);
                    lower.x = std::min(lower.x, r.x);
                    lower.y = std::min(lower.y, r.y);
                    lower.z = std::min(lower.z, r.z);
                    upper.x = std::max(upper.x, r.x);
                    upper.y = std::max(upper.y, r.y);
                    upper.z = std::max(upper.z, r.z);
                    }
                }
            }

        return hoomd::detail::AABB(vec3<Scalar>(lower.x, lower.y, lower.z),
                                   vec3<Scalar>(upper.x, upper.y, upper.z));
        }
    };

namespace detail
    {
//! Export TriangulatedGeometryFiller to python
void export_TriangulatedGeometryFiller(pybind11::module& m);
    } // end namespace detail
    } // end namespace mpcd
    } // end namespace hoomd
#endif // MPCD_TRIANGULATED_GEOMETRY_FILLER_H_
