# Copyright (c) 2009-2026 The Regents of the University of Michigan.
# Part of HOOMD-blue, released under the BSD 3-Clause License.

r"""Virtual particles are MPCD particles that are added to ensure MPCD
collision cells that are sliced by solid boundaries do not become "underfilled".
From the perspective of the MPCD algorithm, the number density of particles in
these sliced cells is lower than the average density, and so the transport
properties may differ. In practice, this usually means that the boundary
conditions do not appear to be properly enforced.

.. invisible-code-block: python

    simulation = hoomd.util.make_example_simulation(mpcd_types=["A"])
    simulation.operations.integrator = hoomd.mpcd.Integrator(dt=0.1)

"""

import hoomd
from hoomd.data.parameterdicts import ParameterDict
from hoomd.mpcd import _mpcd
from hoomd.mpcd.geometry import Geometry, TriangulatedGeometry
from hoomd.operation import Operation
import inspect


class VirtualParticleFiller(Operation):
    """Base virtual-particle filler.

    Args:
        type (str): Type of particles to fill.
        density (float): Particle number density.
        kT (hoomd.variant.variant_like): Temperature of particles.

    Virtual particles will be added with the specified `type` and `density`.
    Their velocities will be drawn from a Maxwell--Boltzmann distribution
    consistent with `kT`.

    .. invisible-code-block: python

        filler = hoomd.mpcd.fill.VirtualParticleFiller(
            type="A",
            density=5.0,
            kT=1.0)
        simulation.operations.integrator.virtual_particle_fillers = [filler]

    {inherited}

    **Members defined in** `VirtualParticleFiller`:

    Attributes:
        density (float): Particle number density.

            .. rubric:: Example:

            .. code-block:: python

                filler.density = 5.0

        kT (hoomd.variant.variant_like): Temperature of particles.

            .. rubric:: Examples:

            Constant temperature.

            .. code-block:: python

                filler.kT = 1.0

            Variable temperature.

            .. code-block:: python

                filler.kT = hoomd.variant.Ramp(1.0, 2.0, 0, 100)

        type (str): Type of particles to fill.

            .. rubric:: Example:

            .. code-block:: python

                filler.type = "A"

    """

    __doc__ = inspect.cleandoc(__doc__).replace(
        "{inherited}", inspect.cleandoc(Operation._doc_inherited)
    )
    _doc_inherited = (
        Operation._doc_inherited
        + """

    **Members inherited from**
    `VirtualParticleFiller <hoomd.mpcd.fill.VirtualParticleFiller>`:

    .. py:attribute:: density

        Particle number density.
        `Read more... <hoomd.mpcd.fill.VirtualParticleFiller.density>`

    .. py:attribute:: kT

        Temperature of particles.
        `Read more... <hoomd.mpcd.fill.VirtualParticleFiller.kT>`

    .. py:attribute:: type

        Type of particles to fill.
        `Read more... <hoomd.mpcd.fill.VirtualParticleFiller.type>`
    """
    )

    def __init__(self, type, density, kT):
        super().__init__()

        param_dict = ParameterDict(
            type=str(type),
            density=float(density),
            kT=hoomd.variant.Variant,
        )
        param_dict["kT"] = kT
        self._param_dict.update(param_dict)


class GeometryFiller(VirtualParticleFiller):
    """Virtual-particle filler for a bounce-back geometry.

    Args:
        type (str): Type of particles to fill.
        density (float): Particle number density.
        kT (hoomd.variant.variant_like): Temperature of particles.
        geometry (hoomd.mpcd.geometry.Geometry): Surface to fill around.
        num_classify_trials (int): Trial points drawn per collision cell when
            determining which cells need to be filled, or `None` to determine
            it from `density`.

    Virtual particles are inserted in cells whose volume is sliced by the
    specified `geometry`. The algorithm for doing the filling depends on the
    specific `geometry`.

    Those fillers first determine which collision cells can possibly contain
    both fluid and solid, then draw virtual particles only in those cells.
    The classification is performed once and repeated only if the box changes.

    .. rubric:: Example:

    Filler for parallel plate geometry.

    .. code-block:: python

        plates = hoomd.mpcd.geometry.ParallelPlates(separation=6.0)
        filler = hoomd.mpcd.fill.GeometryFiller(
            type="A", density=5.0, kT=1.0, geometry=plates
        )
        simulation.operations.integrator.virtual_particle_fillers = [filler]

    {inherited}

    **Members defined in** `GeometryFiller`:

    Attributes:
        geometry (hoomd.mpcd.geometry.Geometry): Surface to fill around
            (*read only*).

        num_classify_trials (int): Trial points drawn per collision cell when
            determining which cells need to be filled (*read only*).

            If `None`, the number is chosen as the mean number of solvent particles
            in the classification region plus eight standard deviations.

    """

    __doc__ = inspect.cleandoc(__doc__).replace(
        "{inherited}", inspect.cleandoc(VirtualParticleFiller._doc_inherited)
    )
    _cpp_class_map = {}

    def __init__(self, type, density, kT, geometry, num_classify_trials=None):
        super().__init__(type, density, kT)

        param_dict = ParameterDict(
            geometry=Geometry,
        )
        param_dict["geometry"] = geometry
        self._param_dict.update(param_dict)
        self._num_classify_trials = int(
            num_classify_trials if num_classify_trials is not None else 0
        )

    @property
    def num_classify_trials(self):
        """int: Trial points drawn per collision cell during classification."""
        return self._num_classify_trials

    def _attach_hook(self):
        sim = self._simulation
        sim._warn_if_seed_unset()

        self.geometry._attach(sim)

        # try to find class in map, otherwise default to internal MPCD module
        geom_type = type(self.geometry)
        try:
            class_info = self._cpp_class_map[geom_type]
        except KeyError:
            class_info = (_mpcd, geom_type.__name__ + "GeometryFiller")
        class_info = list(class_info)
        if isinstance(sim.device, hoomd.device.GPU):
            class_info[1] += "GPU"
        class_ = getattr(*class_info, None)
        assert class_ is not None, "Virtual particle filler for geometry not found"

        self._cpp_obj = class_(
            sim.state._cpp_sys_def,
            self.type,
            self.density,
            self.kT,
            self.geometry._cpp_obj,
            self.num_classify_trials,
        )

        super()._attach_hook()

    def _detach_hook(self):
        self.geometry._detach()
        super()._detach_hook()

    @classmethod
    def _register_cpp_class(cls, geometry, module, cpp_class_name):
        cls._cpp_class_map[geometry] = (module, cpp_class_name)


class TriangulatedGeometryFiller(VirtualParticleFiller):
    """Virtual-particle filler for a triangulated geometry.

    Args:
        type (str): Type of particles to fill.
        density (float): Particle number density.
        kT (hoomd.variant.variant_like): Temperature of particles.
        geometry (hoomd.mpcd.geometry.TriangulatedGeometry): Surface to fill around.
        num_classify_trials (int): Trial points drawn per collision cell when
            determining which cells need to be filled, or `None` to determine
            it from `density`.

    Virtual particles are inserted in cells whose volume is sliced by the
    triangulated mesh. The filler first determines which cells can contain both
    fluid and solid and stores the triangles near each of them, then draws
    virtual particles only in those cells. A drawn particle is kept if it lies
    on the side of the nearest triangle that its normal points to, so the
    triangles must be wound so that their normals point out of the fluid. The
    classification is performed once and repeated only if the box changes.

    .. rubric:: Example:

    Filler for triangulated plate geometry.

    .. code-block:: python

        vertices = numpy.array(
            [[-5, 2.5, -5], [-5, 2.5, 5], [5, 2.5, -5], [5, 2.5, 5]]
        )
        triangles = numpy.array([[0, 1, 3], [0, 3, 2]])
        plate = hoomd.mpcd.geometry.TriangulatedGeometry(
            simulation, vertices, triangles, unwrap_distance=0.0
        )
        filler = hoomd.mpcd.fill.TriangulatedGeometryFiller(
            type="A", density=5.0, kT=1.0, geometry=plate
        )
        simulation.operations.integrator.virtual_particle_fillers = [filler]

    {inherited}

    **Members defined in** `TriangulatedGeometryFiller`:

    Attributes:
        geometry (hoomd.mpcd.geometry.TriangulatedGeometry): Surface to fill around
            (*read only*).

        num_classify_trials (int): Trial points drawn per collision cell when
            determining which cells need to be filled (*read only*).

            If `None`, the number is chosen as the mean number of solvent particles
            in the classification region plus eight standard deviations.

    """

    __doc__ = inspect.cleandoc(__doc__).replace(
        "{inherited}", inspect.cleandoc(VirtualParticleFiller._doc_inherited)
    )

    def __init__(self, type, density, kT, geometry, num_classify_trials=None):
        super().__init__(type, density, kT)

        if not isinstance(geometry, TriangulatedGeometry):
            raise TypeError("Geometry must be a TriangulatedGeometry")

        self._geometry = geometry
        self._num_classify_trials = int(
            num_classify_trials if num_classify_trials is not None else 0
        )

    @property
    def geometry(self):
        """hoomd.mpcd.geometry.TriangulatedGeometry: Surface to fill around."""
        return self._geometry

    @property
    def num_classify_trials(self):
        """int: Trial points drawn per collision cell during classification."""
        return self._num_classify_trials

    def _attach_hook(self):
        sim = self._simulation
        sim._warn_if_seed_unset()

        self._cpp_obj = _mpcd.TriangulatedGeometryFiller(
            sim.state._cpp_sys_def,
            self.type,
            self.density,
            self.kT,
            self.geometry._cpp_obj,
            self.num_classify_trials,
        )

        super()._attach_hook()


__all__ = [
    "GeometryFiller",
    "TriangulatedGeometryFiller",
    "VirtualParticleFiller",
]
