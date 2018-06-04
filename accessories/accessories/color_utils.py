#!/usr/bin/env python3
# -*- coding: utf-8 -*-


"""
Accessory functions for easy and smoothed gradient palettes
"""

from colour import Color

__author__ = "Bulak Arpat"
__copyright__ = "Copyright 2017, Bulak Arpat"
__license__ = "GPLv3"
__version__ = "0.1.1"
__maintainer__ = "Bulak Arpat"
__email__ = "Bulak.Arpat@unil.ch"
__status__ = "Development"


PALETTE_TYPES = ('quantitative', 'qualitative', 'sequential',
                 'diverging', 'unknown')

# TODO: i) Refactor code, ii) Provide Palette versions for gradients

class Palette(object):
    """
    A palette class to contain color objects from colour package
    """
    def __init__(self, name, family, colors=None, pal_type=None):
        if isinstance(name, str):
            self.name = name
        else:
            raise TypeError("'name' of the palette has to be a string type")
        if isinstance(family, str):
            self.family = family
        else:
            raise TypeError("'family' of palette has to be a string type")
        self.colors = []
        if isinstance(colors, (list, tuple)):
            self.colors = [colorize(color) for color in colors]
        if pal_type in PALETTE_TYPES:
            self.pal_type = pal_type
        else:
            self.pal_type = 'unknown'

    def __getattr__(self, label):
        if 'get_' + label in self.__class__.__dict__:
            return getattr(self, 'get_' + label)()
        else:
            raise AttributeError("'%s' not found" % label)

    # getters
    def get_ncols(self):
        return len(self.colors)
    def get_hsl(self):
        return [color.hsl for color in self.colors]
    def get_hex(self):
        return [color.hex for color in self.colors]
    def get_hex_l(self):
        return [color.hex_l for color in self.colors]
    def get_rgb(self):
        return [color.rgb for color in self.colors]

    # pretty
    def __str__(self):
        return "{}".format(self.name)
    def __repr__(self):
        return "<Palette {} ('{}', {} colors)>".format(
            self.name, self.pal_type, self.ncols)


class ColorGradient():
    """
    A class to conveniently return color gradients
    """
    valid_gradients = ('linear1', 'linear2', 'bezier')
    valid_returns = ('hex', 'rgb')
    def __init__(self, color_s, gradient_type='linear1', return_type='hex'):
        self.gradient_type = self._check_gradient_type(gradient_type)
        self.return_type = self._check_return_type(return_type)
        if not color_s:
            color_s = ['White', 'Black']
        if isinstance(color_s, tuple):
            color_s = list(color_s)
        if not isinstance(color_s, list):
            color_s = [color_s]
        color_s = [colorize(c) for c in color_s]
        self.edge_colors = color_s
        self.init_color = color_s[0]
        try:
            self.last_color = color_s[-1]
        except IndexError:
            self.last_color = Color('black')
    def _check_gradient_type(self, gradient_type):
        if gradient_type and gradient_type in self.valid_gradients:
            return gradient_type
        elif hasattr(self, 'gradient_type'):
            return self.gradient_type
        return self.valid_gradients[0]
    def _check_return_type(self, return_type):
        if return_type and return_type in self.valid_returns:
            return return_type
        elif hasattr(self, 'return_type'):
            return self.return_type
        return self.valid_returns[0]
    def gradient(self, num_col=5, gradient_type=None, return_type=None):
        """
        Returns a gradient
        """
        gradient_type = self._check_gradient_type(gradient_type)
        return_type = self._check_return_type(return_type)
        if gradient_type == 'linear1':
            cur_gradient = linear_gradient1(self.init_color,
                                            self.last_color, num_col)
        elif gradient_type == 'linear2':
            cur_gradient = linear_gradient2(self.init_color,
                                            self.last_color, num_col)
        elif gradient_type == 'bezier':
            cur_gradient = bezier_gradient(self.edge_colors, num_col)
        if return_type == 'hex':
            gradient = [c.hex_l for c in cur_gradient]
        elif return_type == 'rgb':
            gradient = [c.rgb for c in cur_gradient]
        return gradient
    def gradient_palette(self, name, num_col=5, gradient_type=None,):
        """
        Returns a gradient palette
        """
        colors = self.gradient(num_col, gradient_type, return_type="hex")
        return Palette(name, gradient_type + "_gradient", colors=colors)


def colorize(color):
    if isinstance(color, Color):
        return color
    if isinstance(color, int):
        return Color(red=color)
    if isinstance(color, str):
        return Color(color)
    if isinstance(color, (tuple, list)):
        if len(color) >= 3:
            return Color(rgb=color[:3])
    raise TypeError("'{}' type can't be colorized.".format(type(color).__name__))


def linear_gradient1(init_color, finish_color=Color('black'), num_col=10):
    """
    Returns a gradient list of (n) colors between two Color objects.
    Based on Ben Southgate's function on 
    https://bsou.io/posts/color-gradients-with-python
    """
    # Starting and ending colors in RGB form
    start = init_color.rgb
    finish = finish_color.rgb
    # Initilize a list of the output colors with the starting color
    rgb_list = [start]
    # Calcuate a color at each evenly spaced value of t from 1 to n
    for num in range(1, num_col):
        # Interpolate RGB vector for color at the current value of t
        curr_vector = [start[j] + num / (num_col - 1) * (finish[j] - start[j])
                       for j in range(3)]
        # Add it to our list of output colors
        rgb_list.append(curr_vector)
    return [Color(rgb=rgb) for rgb in rgb_list]


def linear_gradient2(init_color, finish_color=Color('black'), num_col=10):
    """
    Returns a gradient list of (n) colors between two Color objects.
    Based on 'colour' module's internal range_to function
    """
    return init_color.range_to(finish_color, num_col)


FACT_CACHE = {}
def fact(num):
    """
    Memoized factorial function
    """
    try:
        return FACT_CACHE[num]
    except KeyError:
        if num == 1 or num == 0:
            result = 1
        else:
            result = num * fact(num-1)
        FACT_CACHE[num] = result
        return result


def bernstein(frac, num, inx):
    """
    Bernstein coefficient
    """
    binom = fact(num) / (fact(inx) * fact(num - inx))
    return binom * (1 - frac)**(num - inx) * (frac**inx)


def bezier_gradient(colors, num_out=100):
    """
    Returns a "bezier gradient" dictionary using a given list of
    colors as control points. Based on Ben Southgate's function on
    https://bsou.io/posts/color-gradients-with-python
    """
    # RGB vectors for each color, use as control points
    rgb_list = [color.rgb for color in colors]
    num = len(rgb_list) - 1

    def bezier_interp(frac):
        """
        Defines an interpolation function for this specific curve
        """
        # List of all summands
        summands = [list(map(lambda x: bernstein(frac, num, i) * x, c))
                    for i, c in enumerate(rgb_list)]
        # Output color
        out = [0.0, 0.0, 0.0]
        # Add components of each summand together
        for vector in summands:
            for inx in range(3):
                out[inx] += vector[inx]
        return out
    gradient = [bezier_interp(t/(num_out - 1)) for t in range(num_out)]
    return [Color(rgb=rgb) for rgb in gradient]
