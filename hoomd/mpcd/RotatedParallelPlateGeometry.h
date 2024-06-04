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

class __attribute__((visibility("default"))) RotatedParallelPlateGeometry
    {
    public:
    //! Constructor
    /*!
     * \param separation Distance between plates
     * \param angle Rotation angle (degrees)
     * \param speed Speed of the wall
     * \param no_slip Boundary condition at the wall (slip or no-slip)
     */
    HOSTDEVICE
    RotatedParallelPlateGeometry(Scalar separation, Scalar angle, Scalar speed, bool no_slip)
        : m_V(speed), m_no_slip(no_slip)
        {
        m_angle = angle * M_PI / Scalar(180.);
        m_H = separation / (Scalar(2) * slow::cos(m_angle));
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
        const Scalar sin_angle = slow::sin(m_angle);
        const Scalar cos_angle = slow::cos(m_angle);
        const Scalar tan_angle = sin_angle / cos_angle;
        const Scalar y_center = pos.x * tan_angle;
        const Scalar y1 = y_center + m_H;
        const Scalar y2 = y_center - m_H;

        /* Check if particle is in bounds or could not have left because
         * normal speed is zero, and exit early. The normal is
         * n = [-sin(angle), cos(angle)] with no z component.
         */
        const signed char sign = (char)((pos.y > y1) - (pos.y < y2));
        const Scalar normal_speed = vel.x * -sin_angle + vel.y * cos_angle;
        if (sign == 0 || normal_speed == Scalar(0))
            {
            dt = Scalar(0);
            return false;
            }

        /* Find the time remaining when the particle is collided with wall. This time is computed
         * using the distance difference between pos.y and y_wall in the normal direction divided by
         * normal speed. The distance difference = (pos.y - y_wall) * cos_angle.
         */
        const Scalar y_wall = (sign == 1) ? y1 : y2;
        dt = (pos.y - y_wall) * cos_angle / normal_speed;

        // backtrack the particle for dt to get to point of contact
        pos -= vel * dt;

        // update velocity according to bounday conditions
        // no-slip requires reflection of the tangential components, so vel = -vel + 2V.
        if (m_no_slip)
            {
            const Scalar dV = Scalar(sign * 2) * m_V;
            vel.x = -vel.x + dV * cos_angle;
            vel.y = -vel.y + dV * sin_angle;
            vel.z = -vel.z;
            }
        // slip requires reflection of the normal components, so vel = vel - 2n,
        else
            {
            const Scalar dv = Scalar(2) * normal_speed;
            vel.x -= dv * -sin_angle;
            vel.y -= dv * cos_angle;
            }

        return true;
        }

    //! Check if a particle is out of bounds
    /*!
     * \param pos Current particle position
     * \returns True if particle is out of bounds, and false otherwise
     */
    HOSTDEVICE bool isOutside(const Scalar3& pos) const
        {
        const Scalar tan_angle = slow::tan(m_angle);
        const Scalar y_center = pos.x * tan_angle;
        const Scalar y1 = y_center + m_H;
        const Scalar y2 = y_center - m_H;
        return (pos.y > y1 || pos.y < y2);
        }

    //! Get distance between plates
    /*!
     * \returns Distance between plates
     */
    HOSTDEVICE Scalar getSeparation() const
        {
        return Scalar(2) * slow::cos(m_angle) * m_H;
        }

    //! Get the rotation angle of the plates
    HOSTDEVICE Scalar getAngle() const
        {
        return m_angle * Scalar(180.) / M_PI;
        }

    //! Get the wall speed
    /*!
     * \returns Wall speed
     */
    HOSTDEVICE Scalar getSpeed() const
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
    Scalar m_H;           //!< Half of the channel width
    Scalar m_angle;       //!< Rotation angle of plates
    const Scalar m_V;     //!< Velocity of the wall
    const bool m_no_slip; //!< Boundary condition
    };

    } // end namespace mpcd
    } // end namespace hoomd
#undef HOSTDEVICE

#endif // MPCD_ROTATED_PARALLEL_PLATE_GEOMETRY_H_
