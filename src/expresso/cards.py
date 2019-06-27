#!/usr/bin/env python
# -*- coding: utf-8 -*-


from typing import List

from attr import attrs, attrib
from crystals import Element


class Card(object):
    ...


@attrs
class AtomicSpecies(Card):
    atoms: List[Element] = attrib(factory=list)
    pseudopotentials: List[str] = attrib(factory=list)

    @pseudopotentials.validator
    def _check_length(self, attribute, value):
        if not len(value) == len(self.atoms):
            raise ValueError("Length mismatch!")

    def __iter__(self):
        yield self.atoms
        yield self.pseudopotentials


class AtomicPosition(Card):
    ...


class CellParameters(Card):
    ...
