#!/usr/bin/env python
# -*- coding: utf-8 -*-

from typing import List, Dict, Any

import attr
from attr import attrs

from ..typeconversion import to_fortran

__all__ = [
    'Namelist'
]


@attrs(frozen=True)
class Namelist(object):
    _name: str

    @classmethod
    def from_dict(cls, d: Dict[str, Any]):
        return cls(**d)

    def evolve(self, **changes):
        return attr.evolve(self, **changes)

    def to_dict(self):
        return attr.asdict(self)

    @property
    def names(self) -> List[str]:
        return list(self.to_dict().keys())

    def to_qe(self) -> str:
        entries = {key: to_fortran(value) for (key, value) in self.to_dict().items()}
        import textwrap
        return textwrap.dedent(
            f"&{self._name}\n" +
            "\n".join(f"    {key} = {value}" for (key, value) in entries.items()) +
            "\n/\n"
        )

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
