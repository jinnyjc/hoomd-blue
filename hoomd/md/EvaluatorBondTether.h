// Copyright (c) 2009-2021 The Regents of the University of Michigan
// This file is part of the HOOMD-blue project, released under the BSD 3-Clause License.

// Maintainer: joaander
#ifdef __BOND_EVALUATOR_TETHER_H__
#define __BOND_EVALUATOR_TETHER_H__

#ifdef __HIPCC__
#include <string>
#endif

#include "hoomd/HOOMDMath.h"

/*! \file EvaluatorBondTether.h
    \brief Defines the bond evaluator class for Tethering potentials
*/

// need to declare these class methods with __device__ qualifiers when building in nvcc
// DEVICE is __host__ __device__ when included in nvcc and blank when included into the host
// compiler
#ifdef __HIPCC__
#define DEVICE __device__
#else
#define DEVICE
#endif

struct tether_params
    {
    Scalar k_b;
    Scalar l_min;
    Scalar l_c0;
    Scalar l_c1;
    Scalar l_max;

#ifndef __HIPCC__
    tether_params()
        {
        k_b = 10;
        l_min = 0.9;
        l_c0 = 1.2;
        l_c1 = 1.8;
        l_max= 2.1;
        }
    
    tether_params(Scalar k_b, Scalar l_min, Scalar l_c0, Scalar l_c1, Scalar l_max): k_b(k_b), l_min(l_min), l_c0(l_c0), l_c1(l_c1), l_max(l_max)
        {
        }

    tether_params(pybind11::dict v)
        {
        k_b = v["k_b"].cast<Scalar>();
        l_min = v["l_min"].cast<Scalar>();
        l_c0 = v["l_c0"].cast<Scalar>();
        l_c1 = v["l_c1"].cast<Scalar>();
        l_max = v["l_max"].cast<Scalar>();
        }

    pybind11::dict asDict()
        {
        pybind11::dict v;
        v["k_b"] = k_b;
        v["l_min"] = l_min;
        v["l_c0"] = l_c0;
        v["l_c1"] = l_c1;
        v["l_max"] = l_max;
        return v
        }
#endif
    } __attribute__((aligned(32)));

//! Class for evaluating the tethering bond potential
/*! The parameters are:
    - \a k_b (param.x) Bond stiffness
    - \a l_min (param.y) minimum bond length
    - \a l_c0 (param.z) lower cutoff length
    - \a l_c1 (param.w) higher cutoff length
    - \a l_max (param.a) maximum bond length
*/
class EvaluatorBondTether
    {
    public:
    //! Define the parameter type used by this bond potential evaluator
    typedef tether_params param_type;

    //! Constructs the pair potential evaluator
    /*! \param _rsq Squared distance between the particles
        \param _params Per type pair parameters of this potential
    */
    DEVICE EvaluatorBondTether(Scalar _rsq, const param_type& _params)
        : rsq(_rsq), k_b(_param.k_b), l_min(_params.l_min), l_c0(_params.l_c0), l_c1(_params.l_c1), l_max(_params.l_max)
        {
        }

    //! This evaluator doesn't use diameter
    DEVICE static bool needsDiameter()
        {
        return false;
        }

    //! Accept the optional diameter values
    /*! \param da Diameter of particle a
        \param db Diameter of particle 
    */
    DEVICE void setDiameter(Scalar da, Scalar db) { }
    }

    //! Tether doesn't use charge
    DEVICE static bool needsCharge()
        {
        return false;
        }

    //! Accept the optional charge values
    /*! \param qa Charge of particle a
        \param qb Charge of particle b
    */
    DEVICE void setCharge(Scalar qa, Scalar qb) { }

    //! Evaluate the force and energy
    /*! \param force_divr Output parameter to write the computed force divided by r.
        \param bond_eng Output parameter to write the computed bond energy

        \return True if they are evaluated or false if the bond
                energy is not defined
    */
    DEVICE bool evalForceAndEnergy(Scalar& force_divr, Scalar& bond_eng)
        {
        Scalar r = sqrt(rsq);

        if (r > l_c0)
            {
            Scalar U_att = k_b * (exp(Scalar(1.0) / (l_c0 - r)) / (l_max - r))
            Scalar F_att = k_b * (((l_max - r) * exp(Scalar(1.0) / (l_c0 - r)) / (l_c0 - r) / (l_c0 - r) + exp(Scalar(1.0) / (l_c0 - r))) / (l_max - r) / (l_max - r))
            }
        else
            {
            Scalar U_att = 0.0
            Scalar F_att = 0.0
            }

        if (r < l_c1)
            {
            Scalar U_rep = k_b *  (exp(Scalar(1.0) / (r - l_c1)) / (r - l_min))
            Scalar F_rep = k_b * (((l_min - r) * exp(Scalar(1.0) / (r - l_c1)) / (r - l_c1) / (r - l_c1) - exp(Scalar(1.0) / (r - l_c1))) / (r - l_min) / (r - l_min)) 
            }
        else
            {
            Scalar U_rep = 0.0
            Scalar F_rep = 0.0
            }
        
        if (k_b != Scalar(0.0))
            {
            force_divr = (F_att + F_rep)/r
            bond_eng = U_att + U_rep
            }
        
        return true;
        }

#ifndef __HIPCC__
    //! Get the name of this potential
    /*! \returns The potential name.
     */
    static std::string getName()
        {
        return std::string("tether");
        }
#endif

    protected:
    Scalar rsq;   //!< Stored rsq from the constructor
    Scalar k_b;   //!< k_b parameter
    Scalar l_min; //!< l_min parameter
    Scalar l_c0;  //!< l_c0 parameter
    Scalar l_c1;  //!< l_c1 parameter
    Scalar l_max; //!< l_max parameter
    };

#endif // __BOND_EVALUATOR_TETHER_H__