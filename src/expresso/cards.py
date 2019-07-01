#!/usr/bin/env python
# -*- coding: utf-8 -*-


from typing import List, Optional

import attr
from attr import attrs, attrib
from crystals import Element


@attrs
class Card(object):
    option: Optional[str] = attrib(validator=attr.validators.optional(attr.validators.instance_of(str)))


@attrs
class AtomicSpecies(Card):
    option = attrib(default=None)
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
    option = attrib(default="alat", validator=attr.validators.optional(attr.validators.instance_of(str)))


@attrs
class CellParameters(Card):
    option = attrib(default="alat", validator=attr.validators.optional(attr.validators.instance_of(str)))
