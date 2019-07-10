#!/usr/bin/env python
# -*- coding: utf-8 -*-


from typing import Optional

import attr
from attr import attrs, attrib

__all__ = [
    'Card'
]


@attrs(frozen=True)
class Card(object):
    _name: str
    option: Optional[str]

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
