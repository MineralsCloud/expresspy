#!/usr/bin/env python
# -*- coding: utf-8 -*-


import re
import textwrap
from typing import List, Optional

import attr
import numpy as np
from attr import attrib, attrs
from crystals import Lattice
from expresspy.cards.base import Card
from expresspy.typeconversion import to_fortran
from singleton_decorator import singleton

__all__ = [
    'LatticeParameters',
    'SpecialKPoint',
    'MonkhorstPackGrid',
    'KPointsCard',
    'AtomicSpeciesCard',
    'AtomicPositionCard',
    'CellParametersCard'
]


@attrs(frozen=True)
class AtomicSpecies(object):
    atom: str = attrib(converter=str)
    mass: float = attrib(converter=float)
    pseudopotential: str = attrib(converter=str)

    def to_qe(self) -> str:
        return f"{self.atom}  {to_fortran(self.mass)}  {self.pseudopotential}"


@attrs
class AtomicSpeciesCard(Card):
    _name: str = "ATOMIC_SPECIES"
    option = None
    data: List[AtomicSpecies] = attrib(factory=list)

    def __getitem__(self, item):
        return self.data.__getitem__(item)

    def __len__(self):
        return self.data.__len__()

    def to_qe(self):
        return "\n".join(f"{x}" for x in self.data)


@attrs(frozen=True)
class AtomicPosition(object):
    atom: str = attrib(converter=str)
    position: List[float] = attrib(factory=list, validator=lambda x: len(x) == 3)

    def to_qe(self) -> str:
        return f"{self.atom}  " + "  ".join(map(to_fortran, self.position))


@attrs
class AtomicPositionCard(Card):
    _name: str = "ATOMIC_POSITIONS"
    _allowed_options = ("alat", "bohr", "angstrom", "crystal", "crystal_sg")
    option: str = attrib(default="alat", validator=attr.validators.in_(_allowed_options))
    data: List[AtomicPosition] = attrib(factory=list)

    def to_qe(self):
        return "\n".join(f"{x}" for x in self.data)


@attrs(frozen=True)
class CellParameters(object):
    a = attrib(converter=float)
    b = attrib(converter=float)
    c = attrib(converter=float)
    alpha = attrib(converter=float)
    beta = attrib(converter=float)
    gamma = attrib(converter=float)

    @property
    def edges(self):
        return [self.a, self.b, self.c]

    @property
    def angles(self):
        return [self.alpha, self.beta, self.gamma]

    def to_tuple(self):
        return attr.astuple(self)

    def __getitem__(self, item):
        return self.to_tuple().__getitem__(item)

    def __len__(self):
        return self.to_tuple().__len__()


@attrs
class CellParametersCard(Card):
    _name: str = "CELL_PARAMETERS"
    _allowed_options = ("alat", "bohr", "angstrom")
    option: str = attrib(default="alat", validator=attr.validators.in_(_allowed_options))
    lattice: Lattice = attrib(default=Lattice(np.diag([1, 1, 1])))

    @classmethod
    def from_parameters(cls, option: str, a: float, b: float, c: float, alpha: float, beta: float, gamma: float):
        return cls(option, Lattice.from_parameters(a, b, c, alpha, beta, gamma))

    @classmethod
    def from_array(cls, option: str, array: np.ndarray):
        if not array.shape == (3, 3):
            raise ValueError(f"Expected array of shape (3, 3), given {array.shape}!")
        return cls(option, Lattice(array))

    @property
    def lattice_parameters(self):
        return CellParameters(*self.lattice.lattice_parameters)

    @property
    def lattice_system(self):
        return Lattice.from_parameters(*self.lattice_parameters.to_tuple()).lattice_system

    @property
    def lattice_vectors(self):
        return Lattice.from_parameters(*self.lattice_parameters.to_tuple()).lattice_vectors

    @property
    def periodicity(self):
        return Lattice.from_parameters(*self.lattice_parameters.to_tuple()).periodicity

    @property
    def reciprocal_lattice(self):
        return Lattice.from_parameters(*self.lattice_parameters.to_tuple()).reciprocal

    @property
    def reciprocal_vectors(self):
        return Lattice.from_parameters(*self.lattice_parameters.to_tuple()).reciprocal_vectors

    @property
    def volume(self):
        return Lattice.from_parameters(*self.lattice_parameters.to_tuple()).volume

    def to_qe(self):
        return re.sub("[\[\]]", ' ',
                      np.array2string(self.lattice_vectors, formatter={'float_kind': lambda x: "{:20.10f}".format(x)}))


@attrs
class MonkhorstPackGrid(object):
    grid = attrib(validator=attr.validators.deep_iterable(
        member_validator=attr.validators.instance_of(int),
        iterable_validator=attr.validators.instance_of(list)
    ))
    offsets = attrib(validator=attr.validators.deep_iterable(
        member_validator=attr.validators.instance_of(int),
        iterable_validator=attr.validators.instance_of(list)
    ))

    @property
    def nk1(self):
        return self.grid[0]

    @property
    def nk2(self):
        return self.grid[1]

    @property
    def nk3(self):
        return self.grid[2]

    @property
    def sk1(self):
        return self.offsets[0]

    @property
    def sk2(self):
        return self.offsets[1]

    @property
    def sk3(self):
        return self.offsets[2]


@singleton
class GammaPoint(object):
    ...


@attrs
class SpecialKPoint(object):
    coordinates = attrib(validator=attr.validators.deep_iterable(
        member_validator=attr.validators.instance_of(float),
        iterable_validator=attr.validators.instance_of(list)
    ))
    weight = attrib(validator=attr.validators.instance_of(float))

    @property
    def x(self):
        return self.coordinates[0]

    @property
    def y(self):
        return self.coordinates[1]

    @property
    def z(self):
        return self.coordinates[2]

    @property
    def w(self):
        return self.weight


@attrs
class KPointsCard(Card):
    _name: str = "K_POINTS"
    _allowed_options = ("tpiba", "automatic", "crystal", "gamma", "tpiba_b", "crystal_b", "tpiba_c", "crystal_c")
    option: str = attrib(default="tpiba", validator=attr.validators.in_(_allowed_options))
    points = attrib(factory=list)

    @points.validator
    def _check_points(self, attribute, value):
        if self.option == "automatic":
            if not isinstance(value, MonkhorstPackGrid):
                raise TypeError
        if self.option == "gamma":
            if not isinstance(value, GammaPoint):
                raise TypeError
        else:  # `self.option` in `("tpiba", "crystal", "tpiba_b", "crystal_b", "tpiba_c", "crystal_c")`
            if not all(map(lambda x: isinstance(x, SpecialKPoint), value)):
                raise TypeError

    @classmethod
    def from_monkhorst_pack_grid(cls, mp: MonkhorstPackGrid):
        return cls("automatic", mp)

    @classmethod
    def from_gamma_point(cls, g: GammaPoint):
        return cls("gamma", g)

    def __len__(self) -> Optional[int]:
        if self.option == "automatic":
            return None
        if self.option == "gamma":
            return 1
        else:  # `self.option` in `("tpiba", "crystal", "tpiba_b", "crystal_b", "tpiba_c", "crystal_c")`
            return self.points.__len__()

    def to_qe(self):
        if self.option == "gamma":
            return ""
        if self.option == "automatic":
            return "{} {}".format(self.points.grid, self.points.offsets)
        else:
            return textwrap.dedent(f"""\
            {self._name} {self.option}
            {len(self)}
            """ + f"{point}\n" for point in self.points)
