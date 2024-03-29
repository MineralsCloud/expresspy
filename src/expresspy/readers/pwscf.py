#!/usr/bin/env python
# -*- coding: utf-8 -*-


import io
import operator
import re
import warnings
from collections import namedtuple, OrderedDict
from typing import List, Tuple, MutableMapping, Optional

import f90nml
import numpy as np

from ..cards import MonkhorstPackGrid, AtomicSpecies, AtomicPosition, CellParametersCard, \
    KPointsCard, GammaPoint, AtomicSpeciesCard, AtomicPositionCard
from ..inputs import PWscfInput


def namelist_identifiers() -> Tuple[str, ...]:
    return '&CONTROL', '&SYSTEM', '&ELECTRONS', '&IONS', '&CELL'


def card_identifiers() -> Tuple[str, ...]:
    return 'ATOMIC_SPECIES', 'ATOMIC_POSITIONS', 'K_POINTS', 'CELL_PARAMETERS', 'OCCUPATIONS', 'CONSTRAINTS', 'ATOMIC_FORCES'


class RangeIndices(namedtuple('RangeIndices', ['begin', 'end'])):
    def __str__(self) -> str:
        return "'begin: {0}, end: {1}'".format(self.begin, self.end)


class PWscfInputReader(object):
    """
    This class reads a standard Quantum ESPRESSO PWscf input file or string in, and read it.
    """

    def __init__(self, inp: Optional[str] = None):
        self.newline = "[\r\n,]"  # TODO: This will fail when ',' is inside a value of a parameter.
        self.namelist_sep = r"/\s*[\r\n]"
        self.input = inp

    @property
    def text_content(self):
        with open(self.input) as f:
            return io.StringIO(f.read()).getvalue()

    def get_namelist_identifier_positions(self, pos: int = 0, include_heading: bool = True,
                                          include_ending: bool = False) -> MutableMapping[str, RangeIndices]:
        """
        For example, a typical returned result will look like

        .. code-block:: python

            {'&CONTROL': 'begin: 1, end: 440', '&SYSTEM': 'begin: 443, end: 629', '&ELECTRONS': 'begin: 632, end: 684'}

        :param pos:
        :param include_heading:
        :param include_ending:
        :return:
        """
        match_records = dict()
        identifiers: Tuple[str, ...] = namelist_identifiers()
        s: str = self.text_content
        for pattern in identifiers:
            # ``re.compile`` will produce a regular expression object, on which we can use its ``search`` method.
            m0 = re.compile(pattern, flags=re.IGNORECASE).search(s, match_records.get(pattern, pos))
            if not m0:
                continue
            m1 = re.compile(self.namelist_sep, flags=re.IGNORECASE).search(s, m0.end())
            match_records[pattern] = RangeIndices(begin={True: m0.start(), False: m0.end()}[include_heading],
                                                  end={True: m1.end(), False: m1.start()}[include_ending])
        # The following one-line code first sort the ``dict`` *positions* by its values, i.e., a ``RangeIndices`` tuple.
        # Then we get a list of tuples, with first entries to be the identifiers and second indices to be
        # the indices.
        # For example, if ``x = {'a': 2, 'b': 3, 'c': 1}``, then
        # >>> sorted(x.items(), key=operator.itemgetter(1))
        # [('c', 1), ('a', 2), ('b', 3)]
        # Tuples are compared lexicographically using comparison of corresponding elements, thus compared with their
        # *begin* entry.
        m: List[Tuple[str, RangeIndices]] = sorted(match_records.items(), key=operator.itemgetter(1))
        return OrderedDict(m)

    def get_card_identifier_positions(self, pos: Optional[int] = None) -> MutableMapping[str, RangeIndices]:
        """


        :param pos:
        :return:
        """
        match_records = dict()
        identifiers: Tuple[str, ...] = card_identifiers()
        s: str = self.text_content
        if not pos:
            pos = list(self.get_namelist_identifier_positions().values())[-1].end
        for pattern in identifiers:
            m = re.compile(pattern, flags=re.IGNORECASE).search(s, match_records.get(pattern, pos))
            if not m:
                continue
            match_records[pattern] = m.start()
        sorted_records = sorted(match_records.items(), key=operator.itemgetter(1))
        keys = list(map(operator.itemgetter(0), sorted_records))
        start_indices = list(map(operator.itemgetter(1), sorted_records))
        end_indices = list(map(lambda x: x - 1, start_indices[1:] + [len(s)]))
        return OrderedDict(zip(keys, (RangeIndices(begin=b, end=e) for b, e in zip(start_indices, end_indices))))

    @property
    def namelists_found(self) -> Tuple[str, ...]:
        """
        Check whether an input contains all the namelists necessary for Quantum ESPRESSO to do computations. If the
        input is validated, all namelists found in the input will be returned.

        :return: All namelists found in the input.
        """
        keys = self.get_namelist_identifier_positions().keys()
        if not isempty({"&CONTROL", "&SYSTEM", "&ELECTRONS"}.difference(set(keys))):
            raise RuntimeError("Namelists must contain 'CONTROL', 'SYSTEM' and 'ELECTRONS'!")
        for i, k in enumerate(keys):
            if not k == namelist_identifiers()[i]:
                # Namelists must appear in the order given below.
                raise RuntimeError("Namelists must be in order 'CONTROL', 'SYSTEM', 'ELECTRONS', 'IONS', 'CELL'!")
        else:
            return tuple(map(lambda s: s[1:], keys))

    @property
    def cards_found(self) -> Tuple[str, ...]:
        """
        Check whether an input contains all the cards necessary for Quantum ESPRESSO to do computations. If the
        input is validated, all cards found in the input will be returned.

        :return: All cards found in the input.
        """
        keys = self.get_card_identifier_positions().keys()
        if not isempty({"ATOMIC_SPECIES", "ATOMIC_POSITIONS", "K_POINTS"}.difference(set(keys))):
            raise RuntimeError("'ATOMIC_SPECIES', 'ATOMIC_POSITIONS', and 'K_POINTS' are necessary cards!")
        else:
            return tuple(keys)

    def get_card(self, identifier) -> Optional[List[str]]:
        if identifier in self.cards_found:
            begin, end = self.get_card_identifier_positions()[identifier]
            return re.split(self.newline, self.text_content[begin:end])
        else:
            warnings.warn("Identifier '{0}' is not found in input!".format(identifier), stacklevel=2)

    def read_namelists(self):
        """
        A generic method to read a namelist.
        Note you cannot write more than one parameter in each line!
        :return: a dictionary that stores the inputted information of the intended card
        """
        return OrderedDict(f90nml.read(self.input))

    def read_atomicspecies(self):
        s: Optional[List[str]] = self.get_card("ATOMIC_SPECIES")
        if not s:  # If the returned result is ``None``.
            warnings.warn("'ATOMIC_SPECIES' not found in input!", stacklevel=2)
        else:
            atomic_species = []
            for line in s[1:]:
                # Skip the title line, any empty line, or a line of comment.
                if not line.strip() or line.strip().startswith('!'):
                    continue
                match = re.match(r"(\S+)\s*(-?\d*\.?\d*)\s*(\S+)\s*", line.strip())
                if match is None:
                    warnings.warn("No match found in the line {0}!".format(line))
                else:
                    name, mass, pseudopotential = match.groups()
                    atomic_species.append(AtomicSpecies(name, mass, pseudopotential))
            return AtomicSpeciesCard(data=atomic_species)

    def read_atomicpositions(self):
        s: Optional[List[str]] = self.get_card("ATOMIC_POSITIONS")
        if not s:  # If the returned result is ``None``.
            warnings.warn("'ATOMIC_POSITIONS' not found in input!", stacklevel=2)
        else:
            atomic_positions = []
            title_line = s[0]
            match = re.match("ATOMIC_POSITIONS\s*(?:[({])?\s*(\w*)\s*(?:[)}])?", title_line, flags=re.IGNORECASE)
            if match is None:
                raise RuntimeError("No match found in the line '{0}'! Something went wrong!".format(title_line))
            option = match.group(1)
            if option == '':
                warnings.warn("No option is found, default option 'alat' will be set! "
                              "Not specifying units is DEPRECATED and will no longer be allowed in the future",
                              category=DeprecationWarning)
                option = 'alat'
            for line in s[1:]:
                # If this line is an empty line or a line of comment.
                if line.strip() == '' or line.strip().startswith('!'):
                    continue
                if line.strip() == '/':
                    raise RuntimeError('Do not start any line in cards with a "/" character!')
                if re.match("{.*}", line):
                    match = re.match(
                        "(\w+)\s*(-?\d+\.\d+)\s*(-?\d+\.\d+)\s*(-?\d+\.\d+)\s*{\s*([01])?\s*([01])?\s*([01])?\s*}",
                        line.strip())
                    name, x, y, z, if_pos1, if_pos2, if_pos3 = match.groups()
                    atomic_positions.append(AtomicPosition(atom=name, position=[x, y, z]))
                else:
                    match = re.match("(\w+)\s*(-?\d+\.\d+)\s*(-?\d+\.\d+)\s*(-?\d+\.\d+)", line.strip())
                    if match is None:
                        warnings.warn("No match found in the line {0}!".format(line))
                    else:
                        name, x, y, z = match.groups()
                        atomic_positions.append(AtomicPosition(atom=name, position=[x, y, z]))
            return AtomicPositionCard(option=option, data=atomic_positions)

    # TODO: finish this method
    def read_kpoints(self):
        """
        Find 'K_POINTS' line in the file, and read the k-mesh.
        We allow options and comments on the same line as 'K_POINTS':
        >>> test_strs = ['K_POINTS { crystal }','K_POINTS {crystal}','K_POINTS  crystal','K_POINTScrystal', \
        'K_POINTScrystal! This is a comment.','K_POINTS ! This is a comment.','K_POINTS']
        >>> [re.match("K_POINTS\s*{?\s*(\w*)\s*}?", s, re.IGNORECASE).group(1) for s in test_strs]
        ['crystal', 'crystal', 'crystal', 'crystal', 'crystal', '', '']
        :return: a named tuple defined above
        """
        s: Optional[List[str]] = self.get_card("K_POINTS")
        title_line = s[0]
        match = re.match("K_POINTS\s*(?:[({])?\s*(\w*)\s*(?:[)}])?", title_line, flags=re.IGNORECASE)
        if match is None:
            raise RuntimeError("Match not found! Check your option!")
        option = match.group(1)  # The first parenthesized subgroup will be `option`.
        if option == '':
            raise RuntimeError("Option is not given! you must give one!")
        elif option == 'gamma':
            return KPointsCard(option=option, points=GammaPoint())
        elif option == 'automatic':
            for line in s[1:]:
                if line.strip() == '' or line.strip().startswith('!'):
                    continue
                if line.strip() == '/':
                    raise RuntimeError('Do not start any line in cards with a "/" character!')
                line = line.split()
                grid, offsets = list(map(int, line[0:3])), list(map(int, line[3:7]))
                return KPointsCard(option=option, points=MonkhorstPackGrid(grid=grid, offsets=offsets))
        elif option in {'tpiba', 'crystal', 'tpiba_b', 'crystal_b', 'tpiba_c', 'crystal_c'}:
            return NotImplemented
        else:
            raise ValueError("Unknown option '{0}' given!".format(option))

    def read_cellparameters(self):
        """
        Read 3 lines that follows 'CELL_PARAMETERS' string, so there must be no empty line between 'CELL_PARAMETERS' and
        the real cell parameters!
        :return: a numpy array that stores the cell parameters
        """
        if not self.get_card('CELL_PARAMETERS'):  # If returned result is ``None``.
            warnings.warn("'CELL_PARAMETERS' not found in input!", stacklevel=2)
        else:
            cell_params = []
            title_line = self.get_card('CELL_PARAMETERS')[0]
            match = re.match("CELL_PARAMETERS\s*{?\s*(\w*)\s*}?", title_line, re.IGNORECASE)
            if match is None:
                # The first line should be 'CELL_PARAMETERS blahblahblah', if it is not, either the regular expression
                # wrong or something worse happened.
                raise RuntimeError("No match found! Check you 'CELL_PARAMETERS' line!")
            option = match.group(1)
            if option == '':
                warnings.warn('Not specifying unit is DEPRECATED and will no longer be allowed in the future!',
                              category=DeprecationWarning)
                option = 'bohr'
            for line in self.get_card('CELL_PARAMETERS')[1:]:
                if line.strip() == '/':
                    raise RuntimeError('Do not start any line in cards with a "/" character!')
                if re.match("(-?\d*\.\d*)\s*(-?\d*\.\d*)\s*(-?\d*\.\d*)\s*", line.strip()):
                    v1, v2, v3 = re.match("(-?\d*\.\d*)\s*(-?\d*\.\d*)\s*(-?\d*\.\d*)\s*", line.strip()).groups()
                    cell_params.append([v1, v2, v3])
            return CellParametersCard.from_array(option, np.array(cell_params))

    def create_input_object(self):
        namelists: MutableMapping = self.read_namelists()
        return PWscfInput(**namelists, atomicspecies=self.read_atomicspecies(),
                          atomicpositions=self.read_atomicpositions(), kpoints=self.read_kpoints(),
                          cellparameters=self.read_cellparameters())


def isempty(iterable):
    if len(iterable) == 0:
        return True
    return False
