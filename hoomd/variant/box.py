# Copyright (c) 2009-2023 The Regents of the University of Michigan.
# Part of HOOMD-blue, released under the BSD 3-Clause License.

"""Implement variants that return box parameters as a function of time."""

from hoomd import _hoomd
from hoomd.data.typeconverter import box_preprocessing


class Box(_hoomd.VectorVariantBox):
    """Box-like vector variant base class.

    Provides an interface to length-6 vector variants. The return value of the
    ``__call__`` method returns a length-6 array of scalar values representing
    the ``Box`` attributes ``Lx``, ``Ly``, ``Lz``, ``xy``, ``xz``, and ``yz``.
    """
    pass


class Constant(_hoomd.VectorVariantBoxConstant):
    """A constant box variant.

    Args:
        box (hoomd.box.box_like): The box.
    """

    def __init__(self, box):
        box = box_preprocessing(box)
        _hoomd.VectorVariantBoxConstant.__init__(self, box._cpp_obj)


class Ramp(_hoomd.VectorVariantBoxLinear):
    """Interpolate between two boxes linearly in time.

    Args:
        initial_box (hoomd.box.box_like): The initial box.
        final_box (hoomd.box.box_like): The final box.
        t_start (int): The time step at the start of the ramp.
        t_ramp (int): The length of the ramp.

    Ramp returns the array corresponding to *initial_box* for
    :math:`t \\leq t_{\\mathrm{start}}` and *final_box* for
    :math:`t \\geq t_{\\mathrm{start}} + t_{\\mathrm{ramp}}`.
    """

    def __init__(self, initial_box, final_box, t_start, t_ramp):
        box1 = box_preprocessing(initial_box)
        box2 = box_preprocessing(final_box)
        _hoomd.VectorVariantBoxLinear.__init__(self, box1._cpp_obj,
                                               box2._cpp_obj, t_start, t_ramp)


class LinearInverseVolume(_hoomd.VectorVariantBoxInverseVolumeRamp):
    """Produce dependent box arrays whose inverse volume is linear in time.

    ``LinearInverseVolume`` produces box arrays that correspond to a box whose
    **inverse volume** (i.e., density for a constant number of particles) varies
    linearly with time. The shape of the box remains constant, that is, the
    ratios of the lengths of the box vectors (:math:`L_y / L_x` and
    :math:`L_z / L_x` and the tilt factors (:math:`xy`, :math:`xz`, :math:`yz`)
    remain constant.

    Args:
        initial_box (hoomd.box.box_like): The initial box.
        final_volume (float): The final volume of the box.
        t_start (int): The time step at the start of the ramp.
        t_ramp (int): The length of the ramp.
    """

    def __init__(self, initial_box, final_volume, t_start, t_ramp):
        box = box_preprocessing(initial_box)
        _hoomd.VectorVariantBoxInverseVolumeRamp.__init__(
            self, box._cpp_obj, final_volume, t_start, t_ramp)
