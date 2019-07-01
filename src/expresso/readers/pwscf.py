#!/usr/bin/env python
# -*- coding: utf-8 -*-


import operator
import re
import warnings
from collections import namedtuple, OrderedDict
from typing import *


class RangeIndices(namedtuple('RangeIndices', ['begin', 'end'])):
    def __str__(self) -> str:
        return "'begin: {0}, end: {1}'".format(self.begin, self.end)


class PWscfInputLexer:
    """
    This class reads a standard Quantum ESPRESSO PWscf input file or string in, and lex it.
    """

    def __init__(self, inp: Optional[str] = None, **kwargs):
        self.newline = "[\r\n,]"  # TODO: This will fail when ',' is inside a value of a parameter.
        self.namelist_sep = "/\s*[\r\n]"
        self.__text_stream = open(inp, 'r')

    @property
    def namelist_identifiers(self) -> Tuple[str, ...]:
        return '&CONTROL', '&SYSTEM', '&ELECTRONS', '&IONS', '&CELL'

    @property
    def card_identifiers(self) -> Tuple[str, ...]:
        return 'ATOMIC_SPECIES', 'ATOMIC_POSITIONS', 'K_POINTS', 'CELL_PARAMETERS', 'OCCUPATIONS', 'CONSTRAINTS', 'ATOMIC_FORCES'

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
        identifiers: Tuple[str, ...] = self.namelist_identifiers
        s: str = self.__text_stream.content
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
        identifiers: Tuple[str, ...] = self.card_identifiers
        s: str = self.__text_stream.content
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
    def namelists_found(self) -> Optional[Set[str]]:
        """
        Check whether an input contains all the namelists necessary for Quantum ESPRESSO to do computations. If the
        input is validated, all namelists found in the input will be returned.

        :return: All namelists found in the input.
        """
        keys = list(self.get_namelist_identifier_positions().keys())
        for i, k in enumerate(keys):
            if not k == self.namelist_identifiers[i]:
                # Namelists must appear in the order given below.
                raise RuntimeError("Namelists Must be in order 'CONTROL', 'SYSTEM', 'ELECTRONS', 'IONS', 'CELL'!")
        else:
            return set(map(lambda s: s[1:], keys))

    @property
    def cards_found(self) -> Optional[Set[str]]:
        """
        Check whether an input contains all the cards necessary for Quantum ESPRESSO to do computations. If the
        input is validated, all cards found in the input will be returned.

        :return: All cards found in the input.
        """
        keys = set(self.get_card_identifier_positions().keys())
        if keys < {'ATOMIC_SPECIES', 'ATOMIC_POSITIONS', 'K_POINTS'}:
            # 'ATOMIC_SPECIES', 'ATOMIC_POSITIONS', and 'K_POINTS' are necessary cards.
            warnings.warn('Not enough necessary cards given!')
        else:
            return keys
