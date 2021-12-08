# Copyright (c) 2009-2021 The Regents of the University of Michigan
# This file is part of the HOOMD-blue project, released under the BSD 3-Clause
# License.

"""Test hoomd.hpmc.field"""

import hoomd
import pytest
import numpy as np


valid_constructor_args = [
        dict(
            position=[[0, 0, 0]],
            orientation=[[1, 0, 0, 0]],
            k_translational=1.0,
            k_rotational=1.0,
            symmetries=[[1, 0, 0, 0]]
            ),
]


@pytest.mark.cpu
@pytest.mark.parametrize("constructor_args", valid_constructor_args)
def test_valid_construction_latticefield(device, constructor_args):
    """Test that LatticeField can be constructed with valid arguments."""
    field = hoomd.hpmc.field.LatticeField(**constructor_args)

    # validate the params were set properly
    translator = {
            'position': 'reference_positions',
            'orientation': 'reference_orientations'
    }
    for attr, value in constructor_args.items():
        assert np.all(getattr(field, translator.get(attr, attr)) == value)


@pytest.mark.cpu
def test_attaching(device, simulation_factory, two_particle_snapshot_factory):
    mc = hoomd.hpmc.integrate.Sphere()
    mc.shape['A'] = dict(diameter=0)

    # create simulation & attach objects
    sim = simulation_factory(two_particle_snapshot_factory())
    sim.operations.integrator = mc

    # create lattice field
    lattice = hoomd.hpmc.field.LatticeField(
            position=sim.state.get_snapshot().particles.position,
            orientation=sim.state.get_snapshot().particles.orientation,
            k_translational=1.0,
            k_rotational=1.0,
            symmetries=[[1, 0, 0, 0]])
    mc.external_potential = lattice

    # create C++ mirror classes and set parameters
    sim.run(0)

    # make sure objecst are attached
    assert mc._attached
    assert lattice._attached


@pytest.mark.cpu
def test_detaching(device, simulation_factory, two_particle_snapshot_factory):
    mc = hoomd.hpmc.integrate.Sphere()
    mc.shape['A'] = dict(diameter=0)

    # create simulation & attach objects
    sim = simulation_factory(two_particle_snapshot_factory())
    sim.operations.integrator = mc

    # create lattice field
    lattice = hoomd.hpmc.field.LatticeField(
            position=sim.state.get_snapshot().particles.position,
            orientation=sim.state.get_snapshot().particles.orientation,
            k_translational=1.0,
            k_rotational=1.0,
            symmetries=[[1, 0, 0, 0]])
    mc.external_potential = lattice

    # create C++ mirror classes and set parameters
    sim.run(0)

    # make sure objecst are attached
    sim.operations.remove(mc)
    assert not mc._attached
    assert not lattice._attached

@pytest.mark.cpu
def test_lattice_displacement(device, simulation_factory,
        two_particle_snapshot_factory):
    """Ensure lattice displacements result in expected energy"""
    mc = hoomd.hpmc.integrate.Sphere()
    mc.shape['A'] = dict(diameter=0)

    # create simulation & attach objects
    sim = simulation_factory(two_particle_snapshot_factory())
    sim.operations.integrator = mc

    # create lattice field
    k_trans = 1.0
    lattice = hoomd.hpmc.field.LatticeField(
            position=sim.state.get_snapshot().particles.position,
            orientation=sim.state.get_snapshot().particles.orientation,
            k_translational=k_trans,
            k_rotational=1.0,
            symmetries=[[1, 0, 0, 0]])
    mc.external_potential = lattice

    dx = 0.01
    disp = np.array([dx, 0, 0])
    lattice.reference_positions = lattice.reference_positions + disp

    # create C++ mirror classes and set parameters
    sim.run(0)
    assert np.allclose(
            lattice.energy,
            dx**2 * k_trans * sim.state.N_particles
    )
