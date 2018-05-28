# -*- coding: utf-8 -*-

"""
Various utility functions
"""

from contextlib import contextmanager
import time
import sys

__author__ = "Bulak Arpat"
__copyright__ = "Copyright 2017, Bulak Arpat"
__license__ = "GPLv3"
__version__ = "0.1.0"
__maintainer__ = "Bulak Arpat"
__email__ = "Bulak.Arpat@unil.ch"
__status__ = "Development"


#
# Context manager for timing
#
@contextmanager
def measureTime(title, handle=sys.stderr, unit='second'):
    t1 = time.clock()
    yield
    t2 = time.clock()
    delta = t2 - t1
    if unit == 'minute':
        delta = delta / 60
    elif unit == 'hour':
        delta = delta / 3600
    else:
        unit = 'second'
    handle.write('%s: %0.2f %s(s) elapsed\n' % (title, delta, unit))

#
# Utility function to convert various true/false strings into boolean variables
#
def check_true(s):
    s = s.lower()
    return s in ['true', '1', 't', 'y', 'yes', 'ok']
