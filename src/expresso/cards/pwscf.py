#!/usr/bin/env python
# -*- coding: utf-8 -*-


import textwrap
from typing import Optional

import attr
from attr import attrib, attrs
from expresso.cards.base import Card
from singleton_decorator import singleton

__all__ = [
    'SpecialKPoint',
    'MonkhorstPackGrid',
    'KPoints'
]


@attrs
class MonkhorstPackGrid(object):
    grid = attrib(validator=attr.validators.deep_iterable(
        member_validator=attr.validators.instance_of(int),
        iterable_validator=attr.validators.instance_of(list)
    ))
    offsets = attrib(validator=attr.validators.deep_iterable(
        member_validator=attr.validators.instance_of(int),
        iterable_validator=attr.validators.instance_of(list)
    ))

    @property
    def nk1(self):
        return self.grid[0]

    @property
    def nk2(self):
        return self.grid[1]

    @property
    def nk3(self):
        return self.grid[2]

    @property
    def sk1(self):
        return self.offsets[0]

    @property
    def sk2(self):
        return self.offsets[1]

    @property
    def sk3(self):
        return self.offsets[2]


@singleton
class GammaPoint(object):
    ...


@attrs
class SpecialKPoint(object):
    coordinates = attrib(validator=attr.validators.deep_iterable(
        member_validator=attr.validators.instance_of(float),
        iterable_validator=attr.validators.instance_of(list)
    ))
    weight = attrib(validator=attr.validators.instance_of(float))

    @property
    def x(self):
        return self.coordinates[0]

    @property
    def y(self):
        return self.coordinates[1]

    @property
    def z(self):
        return self.coordinates[2]

    @property
    def w(self):
        return self.weight


@attrs
class KPoints(Card):
    _name: str = "K_POINTS"
    _allowed_options = ("tpiba", "automatic", "crystal", "gamma", "tpiba_b", "crystal_b", "tpiba_c", "crystal_c")
    option: str = attrib(default="tpiba", validator=attr.validators.in_(_allowed_options))
    points = attrib(factory=list)

    @points.validator
    def _check_points(self, attribute, value):
        if self.option == "automatic":
            if not isinstance(value, MonkhorstPackGrid):
                raise TypeError
        if self.option == "gamma":
            if not isinstance(value, GammaPoint):
                raise TypeError
        else:  # `self.option` in `("tpiba", "crystal", "tpiba_b", "crystal_b", "tpiba_c", "crystal_c")`
            if not all(map(lambda x: isinstance(x, SpecialKPoint), value)):
                raise TypeError

    @classmethod
    def from_monkhorst_pack_grid(cls, mp: MonkhorstPackGrid):
        return cls("automatic", mp)

    @classmethod
    def from_gamma_point(cls, g: GammaPoint):
        return cls("gamma", g)

    def __len__(self) -> Optional[int]:
        if self.option == "automatic":
            return None
        if self.option == "gamma":
            return 1
        else:  # `self.option` in `("tpiba", "crystal", "tpiba_b", "crystal_b", "tpiba_c", "crystal_c")`
            return self.points.__len__()

    def to_qe(self):
        if self.option == "gamma":
            return ""
        if self.option == "automatic":
            return "{} {}".format(self.points.grid, self.points.offsets)
        else:
            return textwrap.dedent(f"""\
            {self._name} {self.option}
            {len(self)}
            """ + f"{point}\n" for point in self.points)
