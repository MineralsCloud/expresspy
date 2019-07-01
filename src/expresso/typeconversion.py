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
    if not isinstance(v, int):
        raise TypeError("Expected input to be a int!")
    return str(v)


@to_fortran.register(float)
def _(v: float) -> str:
    if not isinstance(v, float):
        raise TypeError("Expected input to be a float!")
    return "{:.16E}".format(v)


@to_fortran.register(complex)
def _(v: complex) -> str:
    if not isinstance(v, complex):
        raise TypeError("Expected input to be a complex!")
    return "CMPLX({}, {})".format(to_fortran(v.real), to_fortran(v.imag))


@to_fortran.register(bool)
def _(v: bool) -> str:
    if not isinstance(v, bool):
        raise TypeError("Expected input to be a boolean!")
    if v is True:
        return ".true."
    return ".false."


@to_fortran.register(str)
def _(v: str) -> str:
    if not isinstance(v, str):
        raise TypeError("Expected input to be a string!")
    return v


@to_fortran.register(list)
def _(v: list) -> List[str]:
    return list(map(to_fortran, v))
