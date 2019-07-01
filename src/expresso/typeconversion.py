#!/usr/bin/env python
# -*- coding: utf-8 -*-

from functools import singledispatch
from typing import List

__all__ = [
    'to_fortran'
]


@singledispatch
def to_fortran(v):
    raise TypeError("Unsupported type!")


@to_fortran.register(int)
def _(v: int) -> str:
    return str(v)


@to_fortran.register(float)
def _(v: float) -> str:
    return "{:.16E}".format(v)


@to_fortran.register(complex)
def _(v: complex) -> str:
    return "CMPLX({}, {})".format(to_fortran(v.real), to_fortran(v.imag))


@to_fortran.register(bool)
def _(v: bool) -> str:
    if v is True:
        return ".true."
    return ".false."


@to_fortran.register(str)
def _(v: str) -> str:
    return v


@to_fortran.register(list)
def _(v: list) -> List[str]:
    return list(map(to_fortran, v))
