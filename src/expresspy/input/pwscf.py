#!/usr/bin/env python
# -*- coding: utf-8 -*-

from typing import Dict

import attr
from attr import attrib, attrs
from expresspy.cards import AtomicPosition, AtomicSpecies, Card, KPoints
from expresspy.namelists import (ControlNamelist, ElectronsNamelist, Namelist,
                                 SystemNamelist)

__all__ = [
    'PWscfInput'
]


@attrs
class PWscfInput(object):
    namelists: Dict[str, Namelist] = attrib(factory=dict)
    cards: Dict[str, Card] = attrib(factory=dict)

    @namelists.validator
    def _check_must_contain(self, attribute, value):
        if {"CONTROL", "SYSTEM", "ELECTRONS"}.difference(set(value.keys())) != set():
            raise ValueError
        if not all([isinstance(value["CONTROL"], ControlNamelist),
                    isinstance(value["SYSTEM"], SystemNamelist),
                    isinstance(value["ELECTRONS"], ElectronsNamelist)]):
            raise ValueError

    @cards.validator
    def _check_must_contain(self, attribute, value):
        if {"ATOMIC_SPECIES", "ATOMIC_POSITIONS", "K_POINTS"}.difference(set(value.keys())) != set():
            raise ValueError
        if not all([isinstance(value["ATOMIC_SPECIES"], AtomicSpecies),
                    isinstance(value["ATOMIC_POSITIONS"], AtomicPosition),
                    isinstance(value["K_POINTS"], KPoints)]):
            raise ValueError

    def to_dict(self):
        return attr.asdict(self)
