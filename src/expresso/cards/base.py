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
    'CellParameters'
]

LatticeParameters = namedtuple('LatticeParameters', ['a', 'b', 'c', 'alpha', 'beta', 'gamma'])


@attrs(frozen=True)
class Card(object):
    _name: str
    option: Optional[str] = attrib(default=None, validator=attr.validators.optional(attr.validators.instance_of(str)))

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
