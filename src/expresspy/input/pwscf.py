#!/usr/bin/env python
# -*- coding: utf-8 -*-

import attr
from attr import attrib, attrs
from expresspy.cards import AtomicPositionCard, AtomicSpeciesCard, KPointsCard, CellParametersCard
from expresspy.namelists import ControlNamelist, ElectronsNamelist, SystemNamelist, IonsNamelist, CellNamelist

__all__ = [
    'PWscfInput'
]


@attrs
class PWscfInput(object):
    control: ControlNamelist = attrib(default=ControlNamelist())
    system: SystemNamelist = attrib(default=SystemNamelist())
    electrons: ElectronsNamelist = attrib(default=ElectronsNamelist())
    ions: IonsNamelist = attrib(default=IonsNamelist())
    cell: CellNamelist = attrib(default=CellNamelist())
    atomicspecies: AtomicSpeciesCard = attrib(default=None)
    atomicpositions: AtomicPositionCard = attrib(default=None)
    kpoints: KPointsCard = attrib(default=None)
    cellparameters: CellParametersCard = attrib(default=None)

    def to_dict(self):
        return attr.asdict(self)
