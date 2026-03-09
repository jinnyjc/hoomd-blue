# Copyright (c) 2009-2025 The Regents of the University of Michigan.
# Part of HOOMD-blue, released under the BSD 3-Clause License.

import numpy as np
import pytest

import hoomd
from hoomd.conftest import pickling_check


@pytest.fixture
def snap():
    snap_ = hoomd.Snapshot()
    if snap_.communicator.rank == 0:
        snap_.configuration.box = [10, 10, 10, 0, 0, 0]
        snap_.particles.N = 0
        snap_.particles.types = ["A"]
        snap_.mpcd.N = 1
        snap_.mpcd.types = ["A"]
    return snap_


class TestConcentricCylinders:
    def test_default_init(self, simulation_factory, snap):
        geom = hoomd.mpcd.geometry.ConcentricCylinders(
            inner_radius=2.0, outer_radius=4.0
        )
        assert geom.inner_radius == 2.0
        assert geom.outer_radius == 4.0
        assert geom.angular_speed == 0.0
        assert geom.no_slip

        sim = simulation_factory(snap)
        geom._attach(sim)
        assert geom.inner_radius == 2.0
        assert geom.outer_radius == 4.0
        assert geom.angular_speed == 0.0
        assert geom.no_slip

    def test_nondefault_init(self, simulation_factory, snap):
        geom = hoomd.mpcd.geometry.ConcentricCylinders(
            inner_radius=2.0, outer_radius=5.0, angular_speed=1.0, no_slip=False
        )
        assert geom.inner_radius == 2.0
        assert geom.outer_radius == 5.0
        assert geom.angular_speed == 1.0
        assert not geom.no_slip

        sim = simulation_factory(snap)
        geom._attach(sim)
        assert geom.inner_radius == 2.0
        assert geom.outer_radius == 5.0
        assert geom.angular_speed == 1.0
        assert not geom.no_slip

    def test_pickling(self, simulation_factory, snap):
        geom = hoomd.mpcd.geometry.ConcentricCylinders(
            inner_radius=2.0, outer_radius=4.0
        )
        pickling_check(geom)

        sim = simulation_factory(snap)
        geom._attach(sim)
        pickling_check(geom)


class TestCosineChannel:
    def test_default_init(self, simulation_factory, snap):
        geom = hoomd.mpcd.geometry.CosineChannel(
            amplitude=4.0, repeat_length=10.0, separation=4.0
        )
        assert geom.amplitude == 4.0
        assert geom.repeat_length == 10.0
        assert geom.separation == 4.0
        assert geom.no_slip

        sim = simulation_factory(snap)
        geom._attach(sim)
        assert geom.amplitude == 4.0
        assert geom.repeat_length == 10.0
        assert geom.separation == 4.0
        assert geom.no_slip

    def test_nondefault_init(self, simulation_factory, snap):
        geom = hoomd.mpcd.geometry.CosineChannel(
            amplitude=4.0, repeat_length=10.0, separation=4.0, no_slip=False
        )
        assert geom.amplitude == 4.0
        assert geom.repeat_length == 10.0
        assert geom.separation == 4.0
        assert not geom.no_slip

        sim = simulation_factory(snap)
        geom._attach(sim)
        assert geom.amplitude == 4.0
        assert geom.repeat_length == 10.0
        assert geom.separation == 4.0
        assert not geom.no_slip

    def test_pickling(self, simulation_factory, snap):
        geom = hoomd.mpcd.geometry.CosineChannel(
            amplitude=4.0, repeat_length=10.0, separation=4.0
        )
        pickling_check(geom)

        sim = simulation_factory(snap)
        geom._attach(sim)
        pickling_check(geom)


class TestCosineExpansionContraction:
    def test_default_init(self, simulation_factory, snap):
        geom = hoomd.mpcd.geometry.CosineExpansionContraction(
            expansion_separation=4, contraction_separation=2, repeat_length=10.0
        )

        assert geom.expansion_separation == 4.0
        assert geom.contraction_separation == 2.0
        assert geom.repeat_length == 10.0
        assert geom.no_slip

        sim = simulation_factory(snap)
        geom._attach(sim)
        assert geom.expansion_separation == 4.0
        assert geom.contraction_separation == 2.0
        assert geom.repeat_length == 10.0
        assert geom.no_slip

    def test_nondefault_init(self, simulation_factory, snap):
        geom = hoomd.mpcd.geometry.CosineExpansionContraction(
            expansion_separation=4,
            contraction_separation=2,
            repeat_length=10.0,
            no_slip=False,
        )
        assert geom.expansion_separation == 4.0
        assert geom.contraction_separation == 2.0
        assert geom.repeat_length == 10.0
        assert not geom.no_slip

        sim = simulation_factory(snap)
        geom._attach(sim)
        assert geom.expansion_separation == 4.0
        assert geom.contraction_separation == 2.0
        assert geom.repeat_length == 10.0
        assert not geom.no_slip

    def test_pickling(self, simulation_factory, snap):
        geom = hoomd.mpcd.geometry.CosineExpansionContraction(
            expansion_separation=4, contraction_separation=2, repeat_length=10.0
        )

        pickling_check(geom)

        sim = simulation_factory(snap)
        geom._attach(sim)
        pickling_check(geom)


class TestParallelPlates:
    def test_default_init(self, simulation_factory, snap):
        geom = hoomd.mpcd.geometry.ParallelPlates(separation=8.0)
        assert geom.separation == 8.0
        assert geom.speed == 0.0
        assert geom.no_slip

        sim = simulation_factory(snap)
        geom._attach(sim)
        assert geom.separation == 8.0
        assert geom.speed == 0.0
        assert geom.no_slip

    def test_nondefault_init(self, simulation_factory, snap):
        geom = hoomd.mpcd.geometry.ParallelPlates(
            separation=10.0, speed=1.0, no_slip=False
        )
        assert geom.separation == 10.0
        assert geom.speed == 1.0
        assert not geom.no_slip

        sim = simulation_factory(snap)
        geom._attach(sim)
        assert geom.separation == 10.0
        assert geom.speed == 1.0
        assert not geom.no_slip

    def test_pickling(self, simulation_factory, snap):
        geom = hoomd.mpcd.geometry.ParallelPlates(separation=8.0)
        pickling_check(geom)

        sim = simulation_factory(snap)
        geom._attach(sim)
        pickling_check(geom)


class TestPlanarPore:
    def test_default_init(self, simulation_factory, snap):
        geom = hoomd.mpcd.geometry.PlanarPore(separation=8.0, length=10.0)
        assert geom.separation == 8.0
        assert geom.length == 10.0
        assert geom.no_slip

        sim = simulation_factory(snap)
        geom._attach(sim)
        assert geom.separation == 8.0
        assert geom.length == 10.0
        assert geom.no_slip

    def test_nondefault_init(self, simulation_factory, snap):
        geom = hoomd.mpcd.geometry.PlanarPore(
            separation=10.0, length=8.0, no_slip=False
        )
        assert geom.separation == 10.0
        assert geom.length == 8.0
        assert not geom.no_slip

        sim = simulation_factory(snap)
        geom._attach(sim)
        assert geom.separation == 10.0
        assert geom.length == 8.0
        assert not geom.no_slip

    def test_pickling(self, simulation_factory, snap):
        geom = hoomd.mpcd.geometry.PlanarPore(separation=8.0, length=10.0)
        pickling_check(geom)

        sim = simulation_factory(snap)
        geom._attach(sim)
        pickling_check(geom)


class TestSphere:
    def test_default_init(self, simulation_factory, snap):
        geom = hoomd.mpcd.geometry.Sphere(radius=4.0)
        assert geom.radius == 4.0
        assert geom.no_slip

        sim = simulation_factory(snap)
        geom._attach(sim)
        assert geom.radius == 4.0
        assert geom.no_slip

    def test_nondefault_init(self, simulation_factory, snap):
        geom = hoomd.mpcd.geometry.Sphere(radius=5.0, no_slip=False)
        assert geom.radius == 5.0
        assert not geom.no_slip

        sim = simulation_factory(snap)
        geom._attach(sim)
        assert geom.radius == 5.0
        assert not geom.no_slip

    def test_pickling(self, simulation_factory, snap):
        geom = hoomd.mpcd.geometry.Sphere(radius=4.0)
        pickling_check(geom)

        sim = simulation_factory(snap)
        geom._attach(sim)
        pickling_check(geom)


class TestTriangulatedGeometry:
    def test_triangulated_geometry(self, simulation_factory, snap):
        """Test that TriangulatedGeometry correctly stores vertex/triangle numbers."""
        sim = simulation_factory(snap)

        vertices = np.array(
            [
                [-10, 10, -10],
                [-10, 10, 10],
                [10, 10, -10],
                [10, 10, 10],
                [-10, -10, 10],
                [-10, -10, -10],
                [10, -10, 10],
                [10, -10, -10],
            ]
        )
        triangles = np.array([[0, 1, 3], [0, 2, 3], [4, 5, 7], [4, 6, 7]], dtype=int)

        geom = hoomd.mpcd.geometry.TriangulatedGeometry(sim, vertices, triangles, 
                                                        unwrap_distance=(0.0, 0.0, 0.0),
                                                        no_slip=True)

        with geom.cpu_view as data:
            assert data.vertices.shape[0] == vertices.shape[0]
            assert data.triangles.shape[0] == triangles.shape[0]

            np.testing.assert_array_almost_equal(data.vertices, vertices)
            np.testing.assert_array_almost_equal(data.triangles, triangles)