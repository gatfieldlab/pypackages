# -*- coding: utf-8 -*-

"""
Utilities for collections items
"""


import collections


__author__ = "Bulak Arpat"
__copyright__ = "Copyright 2017, Bulak Arpat"
__license__ = "GPLv3"
__version__ = "0.1.0"
__maintainer__ = "Bulak Arpat"
__email__ = "Bulak.Arpat@unil.ch"
__status__ = "Development"


def fix_depth_recursive_defaultdict(depth=3, t=int):
    """
    Generates a recursive defaultdict object with a given depth. The
    innermost values of the object will of type 't'
    """
    if depth==1:
        return collections.defaultdict(t)
    else:
        return collections.defaultdict(
            lambda: fix_depth_recursive_defaultdict(depth-1, t))


recursive_defaultdict = lambda: collections.defaultdict(recursive_defaultdict)
recursive_defaultdict.__doc__ = "Generates an endless recursive defaultdict object"

def int_factory():
    """
    A default dict of integers. Useful to use within other defaultdicts
    """
    return collections.defaultdict(int)


def dict_factory():
    """
    A default dict of dicts. Useful to use within other defaultdicts
    """
    return collections.defaultdict(dict)
