#!/usr/bin/env python
# -*- coding: utf-8 -*-

from typing import List

import attr
from attr import attrib, attrs
from src.expresso.typeconversion import to_fortran

__all__ = [
    'Namelist'
]


def not_private(attribute, value):
    return not attribute.name.startswith("_")


@attrs
class Namelist(object):
    _name: str = attrib(validator=attr.validators.instance_of(str))

    def asdict(self):
        return attr.asdict(self, filter=not_private)

    @property
    def names(self) -> List[str]:
        return list(self.asdict().keys())

    def to_qe(self) -> str:
        entries = {key: to_fortran(value) for (key, value) in self.asdict().items()}
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
        d = attr.asdict(self)
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
