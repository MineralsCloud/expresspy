#!/usr/bin/env python
# -*- coding: utf-8 -*-


import textwrap
from collections import namedtuple
from typing import List, Optional

import attr
from attr import attrs, attrib
from src.expresso.cards.base import Card

__all__ = [
    'MonkhorstPackGrid',
    'KPoints'
]

MonkhorstPackGrid = namedtuple('MonkhorstPackGrid', ['grid', 'offsets'])


@attrs
class KPoints(Card):
    _name: str = "K_POINTS"
    _allowed_options = ("tpiba", "automatic", "crystal", "gamma", "tpiba_b", "crystal_b", "tpiba_c", "crystal_c")
    option: str = attrib(default="tpiba", validator=attr.validators.in_(_allowed_options))
    points: List = attrib(factory=list)

    @points.validator
    def _check_points(self, attribute, value):
        if self.option == "automatic":
            if not len(value) == 6:
                raise ValueError

    @classmethod
    def from_monkhorst_pack_grid(cls, tuple: MonkhorstPackGrid):
        return cls("automatic", tuple)

    @property
    def grid(self):
        if not self.option == "automatic":
            raise ValueError
        return self.points[0:3]

    @property
    def offsets(self):
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
            return "{} {}".format(self.grid, self.offsets)
        else:
            return textwrap.dedent(f"""\
            {self._name} {self.option}
            {len(self)}
            """ + f"{point}\n" for point in self.points)
