#!/usr/bin/env python
# -*- coding: utf-8 -*-


from typing import List, Optional

import attr
import numpy as np
from attr import attrs, attrib
from crystals import Element, Atom, Lattice


@attrs
class Card(object):
    option: Optional[str] = attrib(default=None, validator=attr.validators.optional(attr.validators.instance_of(str)))


@attrs
class AtomicSpecies(Card):
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


@attrs
class AtomicPosition(Card):
    _allowed_options = ("alat", "bohr", "angstrom", "crystal", "crystal_sg")
    option: str = attrib(default="alat", validator=attr.validators.in_(_allowed_options))
    atoms: List[Atom] = attrib(factory=list)


@attrs
class CellParameters(Card):
    _allowed_options = ("alat", "bohr", "angstrom")
    option: str = attrib(default="alat", validator=attr.validators.in_(_allowed_options))
    lattice: Lattice = attrib(default=Lattice(np.diag([1, 1, 1])))


@attrs
class KPoints(Card):
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
