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
from hoomd.mpcd.geometry import Geometry
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
        num_trials (int): Trial points drawn per collision cell when determining
            which cells need to be filled.
        max_per_cell (int): Maximum number of particles drawn per collision
            cell, or `None` to determine it from `density`.

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

        max_per_cell (int): Maximum number of particles drawn per collision
            cell, or 0 to determine it from `density` (*read only*).

            The number drawn per cell is Poisson distributed with mean equal to
            ``density`` times the cell volume, and this sets how much memory is
            reserved per cell. If the number drawn would exceed it, the count is
            clamped. The automatic bound is eight standard deviations above the
            mean, which is exceeded with probability below 1e-7.

        num_trials (int): Trial points drawn per collision cell when determining
            which cells need to be filled (*read only*).

            A cell is filled only if at least one trial point falls on each side
            of the surface, so this must be large enough to resolve features
            that are thin compared to a collision cell.

    """

    __doc__ = inspect.cleandoc(__doc__).replace(
        "{inherited}", inspect.cleandoc(VirtualParticleFiller._doc_inherited)
    )
    _cpp_class_map = {}

    def __init__(self, type, density, kT, geometry, num_trials=1000, max_per_cell=None):
        super().__init__(type, density, kT)

        param_dict = ParameterDict(
            geometry=Geometry,
        )
        param_dict["geometry"] = geometry
        self._param_dict.update(param_dict)
        self._num_trials = int(num_trials)
        self._max_per_cell = 0 if max_per_cell is None else int(max_per_cell)

    @property
    def num_trials(self):
        """int: Trial points drawn per collision cell during classification."""
        return self._num_trials

    @property
    def max_per_cell(self):
        """int: Maximum particles drawn per collision cell, 0 for automatic."""
        return self._max_per_cell

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
            self.num_trials,
            self.max_per_cell,
        )

        super()._attach_hook()

    def _detach_hook(self):
        self.geometry._detach()
        super()._detach_hook()

    @classmethod
    def _register_cpp_class(cls, geometry, module, cpp_class_name):
        cls._cpp_class_map[geometry] = (module, cpp_class_name)


__all__ = [
    "GeometryFiller",
    "VirtualParticleFiller",
]
