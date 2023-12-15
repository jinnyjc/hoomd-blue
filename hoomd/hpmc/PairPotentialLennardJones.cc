// Copyright (c) 2009-2023 The Regents of the University of Michigan.
// Part of HOOMD-blue, released under the BSD 3-Clause License.

#include "PairPotentialLennardJones.h"

namespace hoomd
    {
namespace hpmc
    {

PairPotentialLennardJones::PairPotentialLennardJones(std::shared_ptr<SystemDefinition> sysdef)
    : PairPotential(sysdef), m_type_param_index(sysdef->getParticleData()->getNTypes()),
      m_params(sysdef->getParticleData()->getNTypes())
    {
    }

LongReal PairPotentialLennardJones::getRCut()
    {
    LongReal r_cut = 0;
    for (const auto& param : m_params)
        {
        r_cut = std::max(r_cut, slow::sqrt(param.r_cut_squared));
        }
    return r_cut;
    }

LongReal PairPotentialLennardJones::energy(const vec3<LongReal>& r_ij,
                                            unsigned int type_i,
                                            const quat<LongReal>& q_i,
                                            LongReal charge_i,
                                            unsigned int type_j,
                                            const quat<LongReal>& q_j,
                                            LongReal charge_j)
    {
    LongReal r_squared = dot(r_ij, r_ij);

    unsigned int param_index = m_type_param_index(type_i, type_j);
    const auto& param = m_params[param_index];
    if (r_squared > param.r_cut_squared)
        return 0;

    LongReal lj1 = param.epsilon_x_4 * param.sigma_6 * param.sigma_6;
    LongReal lj2 = param.epsilon_x_4 * param.sigma_6;

    LongReal r_2_inverse = LongReal(1.0) / r_squared;
    LongReal r_6_inverse = r_2_inverse * r_2_inverse * r_2_inverse;

    LongReal energy = r_6_inverse * (lj1 * r_6_inverse - lj2);

    if (param.mode == shift)
        {
        LongReal r_cut_2_inverse = LongReal(1.0) / param.r_cut_squared;
        LongReal r_cut_6_inverse = r_cut_2_inverse * r_cut_2_inverse * r_cut_2_inverse;
        energy -= r_cut_6_inverse * (lj1 * r_cut_6_inverse - lj2);
        }

    if (param.mode == xplor && r_squared > param.r_on_squared)
        {
        LongReal a = param.r_cut_squared - param.r_on_squared;
        LongReal denominator = a * a * a;

        LongReal b = param.r_cut_squared - r_squared;
        LongReal numerator = b * b
                              * (param.r_cut_squared + LongReal(2.0) * r_squared
                                 - LongReal(3.0) * param.r_on_squared);
        energy *= numerator / denominator;
        }

    return energy;
    }

void PairPotentialLennardJones::setParamsPython(pybind11::tuple typ, pybind11::dict params)
    {
    auto pdata = m_sysdef->getParticleData();
    auto type_i = pdata->getTypeByName(typ[0].cast<std::string>());
    auto type_j = pdata->getTypeByName(typ[1].cast<std::string>());
    unsigned int param_index_1 = m_type_param_index(type_i, type_j);
    m_params[param_index_1] = ParamType(params);
    unsigned int param_index_2= m_type_param_index(type_j, type_i);
    m_params[param_index_2] = ParamType(params);
    }

pybind11::dict PairPotentialLennardJones::getParamsPython(pybind11::tuple typ)
    {
    auto pdata = m_sysdef->getParticleData();
    auto type_i = pdata->getTypeByName(typ[0].cast<std::string>());
    auto type_j = pdata->getTypeByName(typ[1].cast<std::string>());
    unsigned int param_index = m_type_param_index(type_i, type_j);
    return m_params[param_index].asDict();
    }

namespace detail
    {
void export_PairPotentialLennardJones(pybind11::module& m)
    {
    pybind11::class_<PairPotentialLennardJones,
                     PairPotential,
                     std::shared_ptr<PairPotentialLennardJones>>(m, "PairPotentialLennardJones")
        .def(pybind11::init<std::shared_ptr<SystemDefinition>>())
        .def("setParams", &PairPotentialLennardJones::setParamsPython)
        .def("getParams", &PairPotentialLennardJones::getParamsPython);
    }
    } // end namespace detail
    } // end namespace hpmc
    } // end namespace hoomd
