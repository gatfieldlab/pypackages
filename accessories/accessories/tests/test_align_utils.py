"""
Unit test for the align_utils module
"""

from unittest import TestCase
from accessories import align_utils

class TestStrMod(TestCase):
    """

    Sets up the unit tests for `str_mod` function

    """
    def setUp(self):
        self.string = "ABCDEF"
        self.insert1 = "ghijklmn"
        self.insert2 = "gh"
    def test_insert_at_start(self):
        """
        Unit test for insertion at start with defaults
        """
        res = align_utils.str_mod(self.string, 0, self.insert1)
        self.assertEqual(res, "ghijklmnABCDEF")
    def test_insert_at_end(self):
        """
        Unit test for insertion at the end with defaults
        """
        res = align_utils.str_mod(self.string, 6, self.insert1)
        self.assertEqual(res, "ABCDEFghijklmn")
    def test_insert_at_middle(self):
        """
        Unit test for insertion into middle with defaults
        """
        res = align_utils.str_mod(self.string, 3, self.insert1)
        self.assertEqual(res, "ABCghijklmnDEF")
    def test_replace_at_start(self):
        """
        Unit test for insertion at start with replacement
        """
        res = align_utils.str_mod(self.string, 0, self.insert1, replace=True)
        self.assertEqual(res, "ghijklmn")
    def test_replace_at_start2(self):
        """
        Unit test for short insertion at start with replacement
        """
        res = align_utils.str_mod(self.string, 0, self.insert2, replace=True)
        self.assertEqual(res, "ghCDEF")
    def test_replace_at_middle(self):
        """
        Unit test for insertion into middle with replacement
        """
        res = align_utils.str_mod(self.string, 3, self.insert1, replace=True)
        self.assertEqual(res, "ABCghijklmn")
    def test_replace_at_middle2(self):
        """
        Unit test for short insertion into middle with replacement
        """
        res = align_utils.str_mod(self.string, 3, self.insert2, replace=True)
        self.assertEqual(res, "ABCghF")
    def test_replace_at_end(self):
        """
        Unit test for insertion at the end with replacement
        """
        res = align_utils.str_mod(self.string, 6, self.insert1, replace=True)
        self.assertEqual(res, "ABCDEFghijklmn")
    def test_replace_at_end2(self):
        """
        Unit test for short insertion at the end with replacement
        """
        res = align_utils.str_mod(self.string, 6, self.insert2, replace=True)
        self.assertEqual(res, "ABCDEFgh")
    def test_insert_beyond(self):
        """
        Unit test for insertion at an index beyond the string with defaults
        """
        res = align_utils.str_mod(self.string, 12, self.insert1, replace=True)
        self.assertEqual(res, "ABCDEFghijklmn")

class TestCigarMD(TestCase):
    """

    Base class for unit tests for CIGAR and MD parsing functions

    """
    def setUp(self):
        self.valid_triples = [
            ("AAACTCAGATCCATCGGTTTTCCATCTTGGTCAGGGTATCTACCATAGAAGAA",
             "3S6M100N8M1D12M10N12M12S", "14^T24"),
            ("CCTTATCCACCTTCCGCTTTACAGCCTCAATGGCGGGAGCATCTGTTGAG",
             "4S38M12N8M", "4C40T"),
            ("TACTTTTGAGAGTTACAACTTGATATATTTTATCCTGAATATGGTTTCGCAATTGAAGTGCAAGAAGAACAGCATGGAAAAATATATAGAATT",
             "11S63M1I8M3I7M", "24A14G2C3A6G9C14"),
            ("", "", ""),
            ("", "", ""),
            ("", "", "")]