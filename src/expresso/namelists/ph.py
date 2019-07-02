#!/usr/bin/env python
# -*- coding: utf-8 -*-

from typing import List

from attr import attrs, attrib
from src.expresso.namelists import Namelist

__all__ = [
    'INPUTPHNamelist'
]


@attrs
class INPUTPHNamelist(Namelist):
    name: str = attrib("INPUTPH")
    amass: List[float] = attrib(converter=lambda x: list(map(float, x)), factory=list)
    outdir: str = attrib(converter=str, default="./")
    prefix: str = attrib(converter=str, default="pwscf")
    niter_ph: int = attrib(converter=int, default=100)
    tr2_ph: float = attrib(converter=float, default=1e-12)
    alpha_mix: float = attrib(converter=float, default=0.7)
    nmix_ph: int = attrib(converter=int, default=4)
    verbosity: str = attrib(converter=str, default="default")
    reduce_io: bool = attrib(converter=bool, default=False)
    max_seconds: float = attrib(converter=float, default=1e7)
    fildyn: str = attrib(converter=str, default="matdyn")
    fildrho: str = attrib(converter=str, default="")
    fildvscf: str = attrib(converter=str, default="")
    epsil: bool = attrib(converter=bool, default=False)
    lrpa: bool = attrib(converter=bool, default=False)
    lnoloc: bool = attrib(converter=bool, default=False)
    trans: bool = attrib(converter=bool, default=True)
    lraman: bool = attrib(converter=bool, default=False)
    eth_rps: float = attrib(converter=float, default=1e-9)
    eth_ns: float = attrib(converter=float, default=1e-12)
    dek: float = attrib(converter=float, default=1e-3)
    recover: bool = attrib(converter=bool, default=False)
    low_directory_check: bool = attrib(converter=bool, default=False)
    only_init: bool = attrib(converter=bool, default=False)
    qplot: bool = attrib(converter=bool, default=False)
    q2d: bool = attrib(converter=bool, default=False)
    q_in_band_form: bool = attrib(converter=bool, default=False)
    electron_phonon: str = attrib(converter=str, default="")
    lshift_q: bool = attrib(converter=bool, default=False)
    zeu: bool = attrib(converter=bool, default=epsil)
    zue: bool = attrib(converter=bool, default=False)
    elop: bool = attrib(converter=bool, default=False)
    fpol: bool = attrib(converter=bool, default=False)
    ldisp: bool = attrib(converter=bool, default=False)
    nogg: bool = attrib(converter=bool, default=False)
    asr: bool = attrib(converter=bool, default=False)
    ldiag: bool = attrib(converter=bool, default=False)
    lqdir: bool = attrib(converter=bool, default=False)
    search_sym: bool = attrib(converter=bool, default=True)
    nq1: int = attrib(converter=int, default=0)
    nq2: int = attrib(converter=int, default=0)
    nq3: int = attrib(converter=int, default=0)
    nk1: int = attrib(converter=int, default=0)
    nk2: int = attrib(converter=int, default=0)
    nk3: int = attrib(converter=int, default=0)
    k1: int = attrib(converter=int, default=0)
    k2: int = attrib(converter=int, default=0)
    k3: int = attrib(converter=int, default=0)
    start_irr: int = attrib(converter=int, default=1)
    last_irr: int = attrib(converter=int, default=3)
    nat_todo: int = attrib(converter=int, default=0)
    modenum: int = attrib(converter=int, default=0)
    start_q: int = attrib(converter=int, default=1)
    last_q: int = attrib(converter=int, default=10)
    # dvscf_star: str = attrib(converter=str, default=1)
    # drho_star: str = attrib(converter=str, default=1)
