#!/usr/bin/env python
# -*- coding: utf-8 -*-


import re
import textwrap
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
    'CellParameters'
]


@attrs(frozen=True)
class LatticeParameters(object):
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

    def __getitem__(self, i):
        return self.to_tuple()[i]


@attrs(frozen=True)
class Card(object):
    _name: str
    option: Optional[str] = attrib(default=None, validator=attr.validators.optional(attr.validators.instance_of(str)))

    def evolve(self, **changes):
        return attr.evolve(self, **changes)

    def to_dict(self):
        return attr.asdict(self)

    def to_qe(self) -> str:
        ...

    def write(self, filename: str):
        with open(filename, "r+") as f:
            f.write(self.to_qe())

    def dump(self, filename: str):
        d = self.to_dict()
        with open(filename, "r+") as f:
            if filename.endswith(".json"):
                import json
                json.dump(d, f)
            if filename.endswith(".yaml|.yml"):
                import yaml
                try:
                    from yaml import CDumper as Dumper
                except ImportError:
                    from yaml import Dumper
                yaml.dump(d, f, Dumper=Dumper)


@attrs
class AtomicSpecies(Card):
    _name: str = "ATOMIC_SPECIES"
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
        {self._name}
        """ + f"{x}" for x in self.data)


@attrs
class AtomicPosition(Card):
    _name: str = "ATOMIC_POSITIONS"
    _allowed_options = ("alat", "bohr", "angstrom", "crystal", "crystal_sg")
    option: str = attrib(default="alat", validator=attr.validators.in_(_allowed_options))
    atoms: List[Atom] = attrib(factory=list)

    def to_qe(self):
        return textwrap.dedent(f"""\
        {self._name} {self.option}
        """ + f"{x}" for x in self.atoms)


@attrs
class CellParameters(Card):
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
        return LatticeParameters(*self.lattice.lattice_parameters)

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
