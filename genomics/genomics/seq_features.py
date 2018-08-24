#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Sequence features/words library for RUST, peak and similar analyses
"""

import itertools
import numpy as np
from genomics import codons, structure
from specialized import cached_objects


__author__ = "Bulak Arpat"
__copyright__ = "Copyright 2018, Bulak Arpat"
__license__ = "GPLv3"
__version__ = "0.0.1"
__maintainer__ = "Bulak Arpat"
__email__ = "Bulak.Arpat@unil.ch"
__status__ = "Development"


# Genetic code
NTS = codons.NTS
CODONTABLE = codons.CODON_TABLE
STOP_CODONS = codons.STOP_CODONS
# Side chain charges by codons
CHARGE_TABLE = codons.CHARGE_TABLE
CHARGE = codons.CHARGE
# Codon usage
CODON_TAI = codons.CODON_TAI['mouse']
CODON_USAGE = codons.CODON_USAGE['mouse']
CODON_STAI = codons.CODON_STAI['mouse-liver']
# Secondary structure
SECONDARY_STRUCT = structure.SECONDARY_STRUCTURE
TOPOLOGY = structure.TOPOLOGY


class _BaseException(Exception):
    """
    A simple internal Exception class for all custom exceptions in this module
    """
    def __init__(self, m):
        message = "[seq_features] {}".format(m)
        super(_BaseException, self).__init__(message)


class TranslationException(_BaseException):
    """
    A custom Exception class for issues raised during translation.
    All other functions setting up features should call FeatureException
    """


class FeatureException(_BaseException):
    """
    A custom Exception class for issues raised while setting up features
    """


def trans(seq, word_dict, word_size, step_size=3):
    """
    Converts (translates) a sequence into sequence features (words)
    """
    translated = []
    try:
        word_mod = word_size % step_size
    except TypeError:
        raise TranslationException("word_size and step_size has to be 'int'")
    if word_mod == 0:
        word_mod = step_size
    step_range = list(range(step_size + 1 - word_mod))
    try:
        seq_len = len(seq)
    except TypeError:
        raise TranslationException("seq has to be 'str' or 'str'-like obj")
    extra = max(word_size - step_size, 0) + seq_len % step_size
    try:
        for i in range(0, seq_len-extra, step_size):
            translated.append(tuple(
                word_dict[seq[i+word_start:i+word_start+word_size]] for
                word_start in step_range))
    except KeyError as err:
        raise TranslationException(
            "seq must have only letters defined in word_dict: {}".format(err))
    except TypeError as err:
        if str(err)[:15] == "unhashable type":
            raise TranslationException(
                "seq must be sliceable and slices must be hashable")
        raise TranslationException(
            "word_dict must be a 'dict' features indexed by letters of seq")
    return translated


def pass_tru():
    """
    Will return always a word dictinary with 0-sized word
    """
    word_size = 0
    words = []
    word_dict = {}
    return word_size, words, word_dict


def nucleotide(size):
    """
    Nucleotide based features, size specifies the size of each word (nt)
    """
    word_size = size
    words = sorted(
        ["".join(word) for word in itertools.product(NTS, repeat=word_size)])
    word_dict = {c: i for i, c in enumerate(words)}
    return word_size, words, word_dict


def aminoacid(size):
    """
    Aminoacid based features, size specifies the size of each word (aa)
    """
    word_size = 3 * size
    words = sorted(["".join(word) for word in
                    itertools.product(set(CODONTABLE.values()), repeat=size)])
    word_dict = {"".join(codons): words.index("".join(aas)) for codons, aas in
                 [zip(*codes) for codes in itertools.product(CODONTABLE.items(),
                                                             repeat=size)]}
    return word_size, words, word_dict


def topology():
    """
    Topology based feature, feature size is fixed at 1 aa (3 nt on nuc acid seq)
    """
    word_size = 3
    words = sorted(TOPOLOGY.values())
    word_dict = {code: words.index(topo) for code, topo in TOPOLOGY.items()}
    return word_size, words, word_dict


def secondary_structure(size):
    """
    Secondary structure based features, size specifies the size of each word (aa)
    """
    word_size = 3 * size
    words = sorted(["".join(word) for word in itertools.product(
        set(SECONDARY_STRUCT.values()), repeat=size)])
    word_dict = {"".join(codons): words.index("".join(sec)) for codons, sec in
                 [zip(*codes) for codes in itertools.product(
                    SECONDARY_STRUCT.items(), repeat=size)]}
    return word_size, words, word_dict


def gc_content(size):
    """
    GC content based features, size specifies the size of each word (aa)
    """
    word_size = 3 * size
    gc_count = {code: str(code.count('G') + code.count('C')) for code in
                ("".join(codons) for codons in
                 itertools.product(CODONTABLE.keys(), repeat=size))}
    words = sorted(set(gc_count.values()))
    word_dict = {code: words.index(gc_count[code]) for
                 code in gc_count}
    return word_size, words, word_dict


def codon_usage_mean(size, usage_index, bins=7):
    """
    Codon usage based features

    Args:
        size: :obj:`int` specifies the size of each word (aa)
        usage_index: :obj:`dict` of codon -> usage metric {codon usage, codon
                     tRNA adaptation index (tAI), species-specific tAI (stAI)}
        bins: :obj:`int` specificies the number of bins that will be used to
              categorize the words. High values are not tested.
    """
    word_size = 3 * size
    usage = {"".join(words): sum([usage_index[codon] for codon in words]) / size
             for words in itertools.product(usage_index.keys(), repeat=size)}
    vals = np.array(sorted(usage.values()))
    lims = np.percentile(vals, np.arange(0, 100, 100 / bins))
    words_rank = np.searchsorted(lims, vals, side="right") - 1
    words = [str(i + 1) for i in sorted(set(words_rank))]
    word_dict = {word: words_rank[np.where(vals == val)[0][0]] for
                 word, val in usage.items()}
    return word_size, words, word_dict


def charge(size, mean=True):
    """
    Side-chain charge based features. Intended for short words (< 5 aa)

    Args:
        size: :obj:`int` specifies the size of each word (aa)
        mean: :obj:`bool` specifies if output should be the mean charge of the
              stretch (True), or exact charge composition (e.g. PUUNNUU) (False)
    """
    if size > 4:
        raise FeatureException(
            "'charge()' can only handle size < 5. Use, charge_bins()")
    def _conv(charges, mean=True):
        if mean:
            return sum(CHARGE[charge] for charge in charges) / len(charges)
        return "".join(charges)
    word_size = 3 * size
    charges = {"".join(codons): _conv(charges, mean) for codons, charges in
               [zip(*codes) for codes in itertools.product(CHARGE_TABLE.items(),
                                                           repeat=size)]}
    words = [str(i) for i in sorted(set(charges.values()))]
    word_dict = {"".join(codes): words.index(
                 str(charges["".join([CODONTABLE[code] for code in codes])]))
                 for codes in itertools.product(set(CODONTABLE.keys()), repeat=size)}
    return word_size, words, word_dict


def charge_bins(size, bins):
    """
    Side-chain charge based features. Uses the mean charge and categorizes into
    bins. Intended for long streches - memory efficient but slow

    Args:
        size: :obj:`int` specifies the size of each word (aa)
        bins: :obj:`int` specificies the number of bins that will be used to
              categorize the words. If it is set to 0 or equal to the size of
              all possible output values then individual values will be used
              for category labels, if it is smaller, then inclusive interval
              notations of quantiles (n=BINS) will be used.
    """
    word_size = 3 * size
    vals = np.array(sorted(set([sum(charges) / size for charges in
        itertools.combinations_with_replacement(CHARGE.values(), size)])))
    if bins == 0:
        bins = len(vals)
    lims = np.percentile(vals, np.arange(0, 100, 100 / bins))
    words_rank = np.searchsorted(lims, vals, side="right") - 1
    if bins == len(vals):
        words = [str(i) for i in vals]
    else:
        words = ["[" + ", ".join(
            [str(i) for i in vals[np.where(words_rank == val)[0][np.array([0, -1])]]])
                 + "]" for val in sorted(set(words_rank))]
    def _charge_indexer(code):
        _charge = [CHARGE[CHARGE_TABLE[code[i]]] for i in range(0, size)]
        return sum(_charge) / len(_charge)
    charges = cached_objects.CappedDictWithProducer(_charge_indexer, cap=8000000)
    def _word_indexer(code):
        return words_rank[np.where(
            vals == charges["".join([CODONTABLE[code[i:i+3]] for i in
                                     range(0, 3 * size, 3)])])[0][0]]
    word_dict = cached_objects.CappedDictWithProducer(_word_indexer, cap=8000000)
    return word_size, words, word_dict


def charge_aa(size1, size2):
    """
    Side-chain charge + following aa based features. This is an experimental
    feature. Words should be composed of the mean charge of size1 aa + identity
    of size2 aa.
    size1 shouldn't be too high - no binning
    size2 shouldn't be larger than 2 - too many combinations
    Example: mean charge of 5 amino acids (size1=5) followed by 2 amino acids
    A word would be then '0.5KD'
    """
    word_size = (size1 + size2) * 3  # 12 for charge + 3 for 1 aa
    charges = {"".join(words): sum([CHARGE[CHARGE_TABLE[aa]] for aa in words]) / size1 for
               words in itertools.product(CHARGE_TABLE.keys(), repeat=size1)}
    aas = {"".join(codons): "".join(aas) for codons, aas in
           [zip(*codes) for codes in itertools.product(CODONTABLE.items(), repeat=size2)]}
    words = sorted(set([charge + aa for charge, aa in itertools.product([str(i) for i
             in set(charges.values())], set(aas.values()))]))
    def _word_indexer(code):
        return words.index(
            str(charges["".join([CODONTABLE[code[i:i+3]] for i in range(0, 3 * size1, 3)])])
            + "".join([CODONTABLE[code[i:i+3]] for i in range(3 * size1, 3 * (size1 + size2), 3)]))
    word_dict = cached_objects.CappedDictWithProducer(_word_indexer, cap=8000000)
    return word_size, words, word_dict



def setup_words(function):
    """
    Sets up the features (words) given a function. It always returns a tuple of
    word_size, words, word_dict explained below:

    Args:
        function: :obj:`str` describes which method to use to set the words
    Returns:
        word_size: :obj:`int` size of words (in units of original seq)
        words: :obj:`list` of words/features used in new translation
        word_dict :obj:`dict` of original seq string: index of word in words
    """
    choices = {
        "pass": pass_tru,
        "nucleotide": lambda: nucleotide(1),
        "dinucleotide": lambda: nucleotide(2),
        "codon": lambda: nucleotide(3),
        "6mer": lambda: nucleotide(6),
        "9mer": lambda: nucleotide(9),
        "aminoacid": lambda: aminoacid(1),
        "dipeptide": lambda: aminoacid(2),
        "tripeptide": lambda: aminoacid(3),
        "topology": topology,
        "secondary_struct": lambda: secondary_structure(1),
        "secondary_struct_2p": lambda: secondary_structure(2),
        "gc1": lambda: gc_content(1),
        "gc2": lambda: gc_content(2),
        "gc3": lambda: gc_content(3),
        "codon_usage": lambda: codon_usage_mean(1, usage_index=CODON_USAGE),
        "codon_usage2p": lambda: codon_usage_mean(2, usage_index=CODON_USAGE),
        "codon_stai": lambda: codon_usage_mean(1, usage_index=CODON_STAI),
        "codon_stai_2p": lambda: codon_usage_mean(2, usage_index=CODON_STAI),
        "codon_tai": lambda: codon_usage_mean(1, usage_index=CODON_TAI),
        "codon_tai_2p": lambda: codon_usage_mean(2, usage_index=CODON_TAI),
        "charge1a": lambda: charge(1, mean=False),
        "charge2p": lambda: charge(2, mean=False),
        "charge2p_mean": lambda: charge(2, mean=True),
        "charge4p_mean": lambda: charge(4, mean=True),
        "charge5p_mean": lambda: charge_bins(5, bins=10),
        "charge4p_aa": lambda: charge_aa(4, 1)
    }
    if function not in choices:
        raise Exception("Unknown function")
    return choices[function]()


def filter_words(words):
    """
    This functions filters out 'semi-smartly' any words that has a stop
    codon or '*' aa. From structure related words, it filters out words
    containing Turn.
    """
    def _tri_slice(word):
        pos = 0
        sliced = []
        while pos < len(word):
            sliced.append(word[pos:pos+3])
            pos += 3
        return sliced
    if words[0][0] == '*':
        filtered_words = [word for word in words if '*' not in word]
    else:
        filtered_words = []
        for word in words:
            if set(_tri_slice(word)) & STOP_CODONS == set():
                filtered_words.append(word)
    if words[0][:5] == 'Helix':
        filtered_words = [word for word in words if 'Turn' not in word]
                          # and 'Unstructured' not in word)]
    return filtered_words
