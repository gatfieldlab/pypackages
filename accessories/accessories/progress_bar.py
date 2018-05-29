# -*- coding: utf-8 -*-

"""
This module provides a simple progress bar for python CLI apps

Typical usage would be:

    from progress_bar import ProgressBar

    workload = 5000000
    my_bar = ProgressBar(workload, bar_len=50, bar_char='-')
    my_bar.start()
    for progress in range(0, workload, 1024):
        my_bar.update(progress)
    my_bar.finish(workload)
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

    """
    Creates, updates and prints a CLI based progress bar

    Args:
        workload (int): This is the value of the total workload for which the
            progress will be shown in the bar. It is unitless and its value
            depends on the implementation of the bar. For example, it could be
            the file size in kb (2938), full percentage (100), number of days
            in a year (365) etc. What is important is that it has to be an
            integer - floats are rounded and therefore won't work - and, the
            following updates has to be in the same unit to make sense.
        bar_size (int, optional): Number of characaters (bar_char) that will
            make up the prgress bar
        bar_char (:obj:`str`): Character that will be used to draw the bar.
            IMPORTANT: If longer than one, the real size of the bar will be
            proportionally longer.
    Returns:
        None
    """

    def __init__(self, workload, bar_len=40, bar_char='='):
        self.size = workload
        self.barlength = bar_len
        self.barchar = bar_char
        self.current = 0

    def start(self):
        """
        Starts the progress bar at 0. For starting at other values try:
            my_bar.update(some_value)
        """
        self.update(0)

    def updater(self, cur_progress, new_line=False):
        """
        For regular cases, use `update` function. `updater` is lower level
        and provides the `new_line` argument that could be useful when the
        output is piped into a log file.

        Args:
            cur_progress (int): Value of current progress.
            new_line (:obj:`bool`, optional): For normal updating this needs to
                be `False`. However, when the progress bar output will be
                forwarded to a log file, multi-line updating is desired,
                and it should be `True`

        Returns:
            None
        """
        self.current = min(cur_progress, self.size)
        progress = 1.0 * self.current / self.size
        progress_len = int(progress * self.barlength)
        tail = ">"
        eol = ""
        if progress_len == self.barlength:
            tail = ""
        if new_line:
            tail = ""
            eol = "\n"
        bar_str = ("[" + self.barchar * progress_len + tail
                   + " " * (self.barlength - len(tail) - progress_len) + "]")
        sys.stderr.write("\r{:4.0%}{} {:,.0f} {}".format(
            progress, bar_str, self.current, eol))
        sys.stderr.flush()

    def update(self, cur_progress):
        """
        Updates the progress bar to the value of the `cur_progress` and
        prints it.

        Args:
            cur_progress: :obj:`int` value of current progress
        Returns:
            None
        """
        self.updater(cur_progress)

    def finish(self, final_progress):
        """
        Finishes the progress bar with the final value of the progress and
        makes a new line. The final_progress does not to be equal to workload.
        This is intentional; in cases where a download breaks before the
        whole file is downloaded, the final progress will indicate this.

        Args:
            final_progress: :obj:`int` final value of progress
        Returns:
            None
        """
        self.updater(final_progress, new_line=True)


def example():
    """
    This function illustrates the use the of the ProgressBar class
    """
    # Our total work is worth 1839328401 apples
    workload = 1839328401
    # We start with 0 apples and progress in steps of 1024 apples
    cur_progress = 0
    progress_step = 1024
    # We would like to update our bar after batches of 15000 apples
    update_threshold = 15000
    update_load = 0
    # We init our bar, which will be composed of 50 'o' characters
    my_bar = ProgressBar(workload, bar_len=50, bar_char='o')
    # Start is nice to draw an empty bar
    my_bar.start()
    while True:
        cur_progress += progress_step
        update_load += progress_step
        if cur_progress > workload:
            break
        if update_load > update_threshold:
            my_bar.update(cur_progress)
            update_load = 0
    # Finishing provides a full bar and a new line. Don't skip it.
    # Even if the progress is larger than workload, max value displayed will
    # be the workload! We don't go beyond 100%
    my_bar.finish(cur_progress)
