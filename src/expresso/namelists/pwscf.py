#!/usr/bin/env python
# -*- coding: utf-8 -*-

from textwrap import dedent
from typing import List

import attr
from attr import attrib, attrs
from expresso.typeconversion import to_fortran

__all__ = [
    'Namelist',
    'ControlNamelist',
    'SystemNamelist',
    'ElectronsNamelist',
    'IonsNamelist',
    'CellNamelist'
]


@attrs
class Namelist(object):
    name: str = attrib()

    def to_fortran(self) -> str:
        entries = {key: to_fortran(value) for (key, value) in attr.asdict(self).items()}
        return dedent("""\
            &{}
                {}
            /
            """.format(self.name, {f"{key} = {value}" for (key, value) in entries.items()}))

    def write(self, filename: str):
        with open(filename, "r+") as f:
            f.write(self.to_fortran())

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


@attrs
class ControlNamelist(Namelist):
    name: str = attrib("CONTROL")
    calculation: str = attrib(converter=str, default='scf')
    title: str = attrib(converter=str, default=' ')
    verbosity: str = attrib(converter=str, default='low')
    restart_mode: str = attrib(converter=str, default='from_scratch')
    wf_collect: bool = attrib(converter=bool, default=True)
    nstep: int = attrib(converter=int, default=1)
    iprint: int = attrib(converter=int, default=1)
    tstress: bool = attrib(converter=bool, default=False)
    tprnfor: bool = attrib(converter=bool, default=False)
    dt: float = attrib(converter=float, default=20.0)
    outdir: str = attrib(converter=str, default='./')
    wfcdir: str = attrib(converter=str, default='./')
    prefix: str = attrib(converter=str, default='pwscf')
    lkpoint_dir: bool = attrib(converter=bool, default=True)
    max_seconds: float = attrib(converter=float, default=10000000.0)
    etot_conv_thr: float = attrib(converter=float, default=0.0001)
    forc_conv_thr: float = attrib(converter=float, default=0.001)
    disk_io: str = attrib(converter=str, default='medium')
    pseudo_dir: str = attrib(converter=str, default='$ESPRESSO_PSEUDO')
    tefield: bool = attrib(converter=bool, default=False)
    dipfield: bool = attrib(converter=bool, default=False)
    lelfield: bool = attrib(converter=bool, default=False)
    nberrycyc: int = attrib(converter=int, default=1)
    lorbm: bool = attrib(converter=bool, default=False)
    lberry: bool = attrib(converter=bool, default=False)
    gdir: int = attrib(converter=int, default=1)
    nppstr: int = attrib(converter=int, default=1)
    lfcpopt: bool = attrib(converter=bool, default=False)
    gate: bool = attrib(converter=bool, default=False)


@attrs
class SystemNamelist(Namelist):
    name: str = attrib("SYSTEM")
    ibrav: int = attrib(converter=int, default=0)
    celldm: List[float] = attrib(converter=lambda x: list(map(float, x)), factory=list)
    A: float = attrib(converter=float, default=0.0)
    B: float = attrib(converter=float, default=0.0)
    C: float = attrib(converter=float, default=0.0)
    cosAB: float = attrib(converter=float, default=0.0)
    cosAC: float = attrib(converter=float, default=0.0)
    cosBC: float = attrib(converter=float, default=0.0)
    nat: int = attrib(converter=int, default=1)
    ntyp: int = attrib(converter=int, default=1)
    nbnd: int = attrib(converter=int, default=20)
    tot_charge: float = attrib(converter=float, default=0.0)
    starting_charge: List[float] = attrib(converter=lambda x: list(map(float, x)), factory=list)
    tot_magnetization: float = attrib(converter=float, default=-1.0)
    starting_magnetization: List[float] = attrib(converter=lambda x: list(map(float, x)), factory=list)
    ecutwfc: float = attrib(converter=float, default=90.0)
    ecutrho: float = attrib(converter=float, default=360.0)
    ecutfock: float = attrib(converter=float, default=120.0)
    nr1: int = attrib(converter=int, default=24)
    nr2: int = attrib(converter=int, default=24)
    nr3: int = attrib(converter=int, default=24)
    nr1s: int = attrib(converter=int, default=24)
    nr2s: int = attrib(converter=int, default=24)
    nr3s: int = attrib(converter=int, default=24)
    nosym: bool = attrib(converter=bool, default=False)
    nosym_evc: bool = attrib(converter=bool, default=False)
    noinv: bool = attrib(converter=bool, default=False)
    no_t_rev: bool = attrib(converter=bool, default=False)
    force_symmorphic: bool = attrib(converter=bool, default=False)
    use_all_frac: bool = attrib(converter=bool, default=False)
    occupations: str = attrib(converter=str, default='smearing')
    one_atom_occupations: bool = attrib(converter=bool, default=False)
    starting_spin_angle: bool = attrib(converter=bool, default=False)
    degauss: float = attrib(converter=float, default=0.0)
    smearing: str = attrib(converter=str, default='gaussian')
    nspin: int = attrib(converter=int, default=1)
    noncolin: bool = attrib(converter=bool, default=False)
    ecfixed: float = attrib(converter=float, default=0.0)
    qcutz: float = attrib(converter=float, default=0.0)
    q2sigma: float = attrib(converter=float, default=0.1)
    input_dft: str = attrib(converter=str, default='PBE0')
    exx_fraction: float = attrib(converter=float, default=0.25)
    screening_parameter: float = attrib(converter=float, default=0.106)
    exxdiv_treatment: str = attrib(converter=str, default='gygi-baldereschi')
    x_gamma_extrapolation: bool = attrib(converter=bool, default=True)
    ecutvcut: float = attrib(converter=float, default=0.0)
    nqx1: int = attrib(converter=int, default=1)
    nqx2: int = attrib(converter=int, default=1)
    nqx3: int = attrib(converter=int, default=1)
    lda_plus_u: bool = attrib(converter=bool, default=False)
    lda_plus_u_kind: int = attrib(converter=int, default=0)
    Hubbard_U: List[float] = attrib(converter=lambda x: list(map(float, x)), factory=list)
    Hubbard_J0: List[float] = attrib(converter=lambda x: list(map(float, x)), factory=list)
    Hubbard_alpha: List[float] = attrib(converter=lambda x: list(map(float, x)), factory=list)
    Hubbard_beta: List[float] = attrib(converter=lambda x: list(map(float, x)), factory=list)
    Hubbard_J: List[float] = attrib(converter=lambda x: list(map(float, x)), factory=list)
    starting_ns_eigenvalue: float = attrib(converter=float, default=-1.0)
    U_projection_type: str = attrib(converter=str, default='atomic')
    edir: int = attrib(converter=int, default=1)
    emaxpos: float = attrib(converter=float, default=0.5)
    eopreg: float = attrib(converter=float, default=0.1)
    eamp: float = attrib(converter=float, default=0.001)
    angle1: List[float] = attrib(converter=lambda x: list(map(float, x)), factory=list)
    angle2: List[float] = attrib(converter=lambda x: list(map(float, x)), factory=list)
    constrained_magnetization: str = attrib(converter=str, default='none')
    fixed_magnetization: List[float] = attrib(converter=lambda x: list(map(float, x)), factory=list)
    # lambda: float = attrib(converter=float, default=1.0)
    report: int = attrib(converter=int, default=100)
    lspinorb: bool = attrib(converter=bool, default=False)
    assume_isolated: str = attrib(converter=str, default='none')
    esm_bc: str = attrib(converter=str, default='pbc')
    esm_w: float = attrib(converter=float, default=0.0)
    esm_efield: float = attrib(converter=float, default=0.0)
    esm_nfit: int = attrib(converter=int, default=4)
    fcp_mu: float = attrib(converter=float, default=0.0)
    vdw_corr: str = attrib(converter=str, default='none')
    london: bool = attrib(converter=bool, default=False)
    london_s6: float = attrib(converter=float, default=0.75)
    london_c6: List[float] = attrib(converter=lambda x: list(map(float, x)), factory=list)
    london_rvdw: List[float] = attrib(converter=lambda x: list(map(float, x)), factory=list)
    london_rcut: int = attrib(converter=int, default=200)
    ts_vdw_econv_thr: float = attrib(converter=float, default=1e-06)
    ts_vdw_isolated: bool = attrib(converter=bool, default=False)
    xdm: bool = attrib(converter=bool, default=False)
    xdm_a1: float = attrib(converter=float, default=0.6836)
    xdm_a2: float = attrib(converter=float, default=1.5045)
    space_group: int = attrib(converter=int, default=0)
    uniqueb: bool = attrib(converter=bool, default=False)
    origin_choice: int = attrib(converter=int, default=1)
    rhombohedral: bool = attrib(converter=bool, default=True)
    zgate: float = attrib(converter=float, default=0.5)
    relaxz: bool = attrib(converter=bool, default=False)
    block: bool = attrib(converter=bool, default=False)
    block_1: float = attrib(converter=float, default=0.45)
    block_2: float = attrib(converter=float, default=0.55)
    block_height: float = attrib(converter=float, default=0.1)


@attrs
class ElectronsNamelist(Namelist):
    name: str = attrib("ELECTRONS")
    electron_maxstep: int = attrib(converter=int, default=100)
    scf_must_converge: bool = attrib(converter=bool, default=True)
    conv_thr: float = attrib(converter=float, default=1e-06)
    adaptive_thr: bool = attrib(converter=bool, default=False)
    conv_thr_init: float = attrib(converter=float, default=0.001)
    conv_thr_multi: float = attrib(converter=float, default=0.1)
    mixing_mode: str = attrib(converter=str, default='plain')
    mixing_beta: float = attrib(converter=float, default=0.7)
    mixing_ndim: int = attrib(converter=int, default=8)
    mixing_fixed_ns: int = attrib(converter=int, default=0)
    diagonalization: str = attrib(converter=str, default='david')
    ortho_para: int = attrib(converter=int, default=0)
    diago_thr_init: float = attrib(converter=float, default=1e-06)
    diago_cg_maxiter: int = attrib(converter=int, default=400)
    diago_david_ndim: int = attrib(converter=int, default=4)
    diago_full_acc: bool = attrib(converter=bool, default=False)
    efield: float = attrib(converter=float, default=0.0)
    efield_cart: tuple = attrib(converter=tuple, default=(0.0, 0.0, 0.0))
    efield_phase: str = attrib(converter=str, default='none')
    startingpot: str = attrib(converter=str, default='atomic')
    startingwfc: str = attrib(converter=str, default='atomic+random')
    tqr: bool = attrib(converter=bool, default=False)


@attrs
class IonsNamelist(Namelist):
    name: str = attrib("IONS")
    ion_dynamics: str = attrib(converter=str, default='bfgs')
    ion_positions: str = attrib(converter=str, default='default')
    pot_extrapolation: str = attrib(converter=str, default='atomic')
    wfc_extrapolation: str = attrib(converter=str, default='none')
    remove_rigid_rot: bool = attrib(converter=bool, default=False)
    ion_temperature: str = attrib(converter=str, default='not_controlled')
    tempw: float = attrib(converter=float, default=300.0)
    tolp: float = attrib(converter=float, default=100.0)
    delta_t: float = attrib(converter=float, default=1.0)
    nraise: int = attrib(converter=int, default=1)
    refold_pos: bool = attrib(converter=bool, default=False)
    upscale: float = attrib(converter=float, default=100.0)
    bfgs_ndim: int = attrib(converter=int, default=1)
    trust_radius_max: float = attrib(converter=float, default=0.8)
    trust_radius_min: float = attrib(converter=float, default=0.001)
    trust_radius_ini: float = attrib(converter=float, default=0.5)
    w_1: float = attrib(converter=float, default=0.01)
    w_2: float = attrib(converter=float, default=0.5)


@attrs
class CellNamelist(Namelist):
    name: str = attrib("CELL")
    cell_dynamics: str = attrib(converter=str, default='none')
    press: float = attrib(converter=float, default=0.0)
    wmass: float = attrib(converter=float, default=0.001)
    cell_factor: float = attrib(converter=float, default=2.0)
    press_conv_thr: float = attrib(converter=float, default=0.5)
    cell_dofree: str = attrib(converter=str, default='all')
