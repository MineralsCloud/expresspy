#!/usr/bin/env python
# -*- coding: utf-8 -*-


import re
import textwrap
from collections import namedtuple
from typing import List, Optional

import attr
import numpy as np
from attr import attrs, attrib
from crystals import Element, Atom, Lattice

__all__ = [
    'LatticeParameters',
    'Card',
    'AtomicSpecies',
    'AtomicPosition',
    'CellParameters',
    'KPoints'
]

LatticeParameters = namedtuple('LatticeParameters', ['a', 'b', 'c', 'alpha', 'beta', 'gamma'])


@attrs
class Card(object):
    name: str = attrib(attr.validators.instance_of(str))
    option: Optional[str] = attrib(default=None, validator=attr.validators.optional(attr.validators.instance_of(str)))


@attrs
class AtomicSpecies(Card):
    name: str = attrib(default="ATOMIC_SPECIES")
    atoms: List[Element] = attrib(factory=list)
    pseudopotentials: List[str] = attrib(factory=list)

    @pseudopotentials.validator
    def _check_length_equal(self, attribute, value):
        if not len(value) == len(self.atoms):
            raise ValueError("Length mismatch!")

    @property
    def data(self):
        return list(zip(self.atoms, self.pseudopotentials))

    def __iter__(self):
        yield self.atoms
        yield self.pseudopotentials

    def __len__(self):
        return self.data.__len__()

    def to_qe(self):
        return textwrap.dedent(f"""\
        {self.name}
        """ + f"{x}" for x in self.data)


@attrs
class AtomicPosition(Card):
    name: str = attrib(default="ATOMIC_POSITIONS")
    _allowed_options = ("alat", "bohr", "angstrom", "crystal", "crystal_sg")
    option: str = attrib(default="alat", validator=attr.validators.in_(_allowed_options))
    atoms: List[Atom] = attrib(factory=list)

    def to_qe(self):
        return textwrap.dedent(f"""\
        {self.name} {self.option}
        """ + f"{x}" for x in self.atoms)


@attrs
class CellParameters(Card):
    name: str = attrib(default="CELL_PARAMETERS")
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
        return LatticeParameters(*self.lattice.lattice_parameters)

    @property
    def lattice_system(self):
        return Lattice.from_parameters(*self.lattice_parameters).lattice_system

    @property
    def lattice_vectors(self):
        return Lattice.from_parameters(*self.lattice_parameters).lattice_vectors

    @property
    def periodicity(self):
        return Lattice.from_parameters(*self.lattice_parameters).periodicity

    @property
    def reciprocal_lattice(self):
        return Lattice.from_parameters(*self.lattice_parameters).reciprocal

    @property
    def reciprocal_vectors(self):
        return Lattice.from_parameters(*self.lattice_parameters).reciprocal_vectors

    @property
    def volume(self):
        return Lattice.from_parameters(*self.lattice_parameters).volume

    def to_qe(self):
        return re.sub("[\[\]]", ' ',
                      np.array2string(self.lattice_vectors, formatter={'float_kind': lambda x: "{:20.10f}".format(x)}))


@attrs
class KPoints(Card):
    name: str = attrib(default="K_POINTS")
    _allowed_options = ("tpiba", "automatic", "crystal", "gamma", "tpiba_b", "crystal_b", "tpiba_c", "crystal_c")
    option: str = attrib(default="tpiba", validator=attr.validators.in_(_allowed_options))
    points: List = attrib(factory=list)

    @points.validator
    def _check_points(self, attribute, value):
        if self.option == "automatic":
            if not len(value) == 6:
                raise ValueError

    @property
    def mesh(self):
        if not self.option == "automatic":
            raise ValueError
        return self.points[0:3]

    @property
    def shift(self):
        if not self.option == "automatic":
            raise ValueError
        return self.points[3:6]

    def __len__(self) -> Optional[int]:
        if self.option == "automatic":
            return None
        if self.option in ("tpiba", "crystal", "tpiba_b", "crystal_b", "tpiba_c", "crystal_c"):
            return self.points.__len__()
        if self.option == "gamma":
            return 1

    def to_qe(self):
        if self.option == "gamma":
            return ""
        if self.option == "automatic":
            return "{} {}".format(self.mesh, self.shift)
        else:
            return textwrap.dedent(f"""\
            {self.name} {self.option}
            {len(self)}
            """ + f"{point}\n" for point in self.points)
