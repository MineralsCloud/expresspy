#!/usr/bin/env python
# -*- coding: utf-8 -*-

from typing import Dict

from attr import attrib, attrs

from src.expresso.cards import *
from src.expresso.namelists.pwscf import *

__all__ = [
    'PWscfInput'
]


@attrs
class PWscfInput(object):
    namelists: Dict[str, Namelist] = attrib(factory=dict)
    cards: Dict[str, Card] = attrib(factory=dict)

    @namelists.validator
    def _check_contains(self, attribute, value):
        if {"CONTROL", "SYSTEM", "ELECTRONS"}.difference(set(value.keys())) is not {}:
            raise ValueError
        if not all([isinstance(value["CONTROL"], ControlNamelist), isinstance(value["SYSTEM"], SystemNamelist),
                    isinstance(value["ELECTRONS"], ElectronsNamelist)]):
            raise ValueError
