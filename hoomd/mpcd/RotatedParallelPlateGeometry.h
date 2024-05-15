// Copyright (c) 2009-2024 The Regents of the University of Michigan.
// Part of HOOMD-blue, released under the BSD 3-Clause License.

/*!
 * \file mpcd/RotatedParallelPlateGeometry.h
 * \brief Definition of the MPCD rotated slit channel geometry
 */

#ifndef MPCD_ROTATED_PARALLEL_PLATE_GEOMETRY_H_
#define MPCD_ROTATED_PARALLEL_PLATE_GEOMETRY_H_

#include "hoomd/BoxDim.h"
#include "hoomd/HOOMDMath.h"

#ifdef __HIPCC__
#define HOSTDEVICE __host__ __device__ inline
#else
#define HOSTDEVICE inline __attribute__((always_inline))
#include <string>
#endif // __HIPCC__

namespace hoomd
    {
namespace mpcd
    {
//! Rotated parallel plate (slit) geometry
/*!
 */
class __attribute__((visibility("default"))) RotatedParallelPlateGeometry
    {
    public:
    //! Constructor
    /*!
     * \param H Channel half-width
     * \param angle Rotation angle (degrees)
     * \param V Velocity of the wall
     * \param no_slip Boundary condition at the wall (slip or no-slip)
     */
    HOSTDEVICE RotatedParallelPlateGeometry(Scalar H, Scalar angle, Scalar V, bool no_slip)
        : m_H(H), m_angle(angle * M_PI / 180.), m_V(V), m_no_slip(no_slip)
        {
        }

    //! Detect collision between the particle and the boundary
    /*!
     * \param pos Proposed particle position
     * \param vel Proposed particle velocity
     * \param dt Integration time remaining
     *
     * \returns True if a collision occurred, and false otherwise
     *
     * \post The particle position \a pos is moved to the point of reflection, the velocity \a vel
     * is updated according to the appropriate bounce back rule, and the integration time \a dt is
     * decreased to the amount of time remaining.
     */
    HOSTDEVICE bool detectCollision(Scalar3& pos, Scalar3& vel, Scalar& dt) const
        {
        // TODO: fill this in
        const Scalar y1 = pos.x * tan(m_angle * M_PI / 180) + m_H;
        const Scalar y2 = pos.x * tan(m_angle * M_PI / 180) - m_H;

        const signed char sign = (char)((pos.y > y1) - (pos.y < y2));

        if (sign == 0 || (vel.x * tan(m_angle * M_PI) - vel.y) == Scalar(0))
            {
            dt = Scalar(0);
            return false;
            }

        if (sign == +1)
            {
            dt = (pos.y - y1) / (vel.x * tan(m_angle * M_PI) - vel.y);
            }

        if (sign == -1)
            {
            dt = (y2 - pos.y) / (vel.x * tan(m_angle * M_PI) - vel.y);
            }

        pos.x -= vel.x * dt;
        pos.y = sign * m_H;
        pos.z -= vel.z * dt;

        if (m_no_slip)
            {
            vel.x = -vel.x + Scalar(sign * 2) * m_V;
            vel.z = -vel.z;
            }
        // both slip and no-slip have no penetration of the surface
        vel.y = -vel.y;

        return true;
        }

    //! Check if a particle is out of bounds
    /*!
     * \param pos Current particle position
     * \returns True if particle is out of bounds, and false otherwise
     */
    HOSTDEVICE bool isOutside(const Scalar3& pos) const
        {
        // TODO: fill this in
        const Scalar y1 = pos.x * tan(m_angle * M_PI / 180.) + m_H;
        const Scalar y2 = pos.x * tan(m_angle * M_PI / 180.) - m_H;

        return (pos.y > y1 || pos.y < y2);
        }

    //! Get channel half width
    /*!
     * \returns Channel half width
     */
    HOSTDEVICE Scalar getH() const
        {
        return m_H;
        }

    //! Get the rotation angle of the plates
    HOSTDEVICE Scalar getAngle() const
        {
        return m_angle * (180. / M_PI);
        }

    //! Get the wall velocity
    /*!
     * \returns Wall velocity
     */
    HOSTDEVICE Scalar getVelocity() const
        {
        return m_V;
        }

    //! Get the wall boundary condition
    /*!
     * \returns Boundary condition at wall
     */
    HOSTDEVICE bool getNoSlip() const
        {
        return m_no_slip;
        }

#ifndef __HIPCC__
    //! Get the unique name of this geometry
    static std::string getName()
        {
        return std::string("RotatedParallelPlates");
        }
#endif // __HIPCC__

    private:
    const Scalar m_H;     //!< Half of the channel width
    const Scalar m_angle; //!< Rotation angle of plates
    const Scalar m_V;     //!< Velocity of the wall
    const bool m_no_slip; //!< Boundary condition
    };

    } // end namespace mpcd
    } // end namespace hoomd
#undef HOSTDEVICE

#endif // MPCD_ROTATED_PARALLEL_PLATE_GEOMETRY_H_
