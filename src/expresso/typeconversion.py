#!/usr/bin/env python
# -*- coding: utf-8 -*-


__all__ = [
    'fortranfloat',
    'fortrancomplex',
    'fortranbool',
    'fortranstr'
]


def fortranfloat(v: float) -> str:
    if not isinstance(v, float):
        raise TypeError("Expected input to be a float!")
    return "{:.16E}".format(v)


def fortrancomplex(v: complex) -> str:
    return "CMPLX({}, {})".format(fortranfloat(v.real), fortranfloat(v.imag))


def fortranbool(v: bool) -> str:
    if not isinstance(v, bool):
        raise TypeError("Expected input to be a boolean!")
    if v is True:
        return ".true."
    return ".false."


def fortranstr(v: str) -> str:
    return '{}'.format(v)
