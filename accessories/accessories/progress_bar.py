# -*- coding: utf-8 -*-

"""
A simple progress bar for python CLI apps
"""

import sys

__author__ = "Bulak Arpat"
__copyright__ = "Copyright 2017, Bulak Arpat"
__license__ = "GPLv3"
__version__ = "0.1.0"
__maintainer__ = "Bulak Arpat"
__email__ = "Bulak.Arpat@unil.ch"
__status__ = "Development"


class ProgressBar:
    def __init__(self, size, bar=40, barchar='='):
        self.size = size
        self.barlength = bar
        self.barchar = barchar
    def start(self):
        self.update(0)
    def updater(self, c, new_line=False):
        self.current = min(c, self.size)
        p = 1.0 * self.current / self.size
        pg = int(p*self.barlength)
        tail = ">"
        eol = ""
        if pg == self.barlength:
            tail = ""
        if new_line:
            tail = ""
            eol = "\n"
        bar = ("[" + self.barchar*pg + tail +
               " "*(self.barlength-len(tail)-pg) + "]")
        sys.stderr.write("\r{:4.0%}{} {:,.0f} {}".format(p, bar, self.current,
                                                         eol))
        sys.stderr.flush()
    def update(self, c):
        self.updater(c)
    def finish(self, c):
        self.updater(c, new_line=True)

def main():
    c = 183932840
    a = 0.0
    update_step = 0
    bar_length = 40
    p = ProgressBar(c, 50, 'o')
    p.start()
    while True:
        a += 1024
        update_step += 1024
        if a > c:
            a = c
            break
        if update_step > 15000:
            p.update(a)
            update_step = 0
    p.finish(a)
