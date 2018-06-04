# -*- coding: utf-8 -*-

"""
Provides several utility functions to help processing sequence alignmnets

In articular, it facilitates parsing, converting and printing alignments from
structures obtained from a SAM file, i.e CIGAR and MD-tags.

From SAM specification:

Op  Description
‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾
M   alignment match (can be a sequence match or mismatch)
I   insertion to the reference
D   deletion from the reference
N   skipped region from the reference
S   soft clipping (clipped sequences present in SEQ)
H   hard clipping (clipped sequences NOT present in SEQ)
P   padding (silent deletion from padded reference)
=   sequence match
X   sequence mismatch

• H can only be present as the first and/or last operation.
• S may only have H operations between them and the ends of the string.
• For mRNA-to-genome alignment, an N operation represents an intron.
  For other types of alignments, the interpretation of N is not defined.
• Sum of the lengths of the M/I/S/=/X operations shall equal the length of SEQ.
"""


import re


__author__ = "Bulak Arpat"
__copyright__ = "Copyright 2017, Bulak Arpat"
__license__ = "GPLv3"
__version__ = "0.1.1"
__maintainer__ = "Bulak Arpat"
__email__ = "Bulak.Arpat@unil.ch"
__status__ = "Development"


# Compiled regex
MD_REGEX = re.compile(r'([0-9]+|[A-Z]+|\^[A-Z]+)')
CIGAR_REGEX = re.compile(r'[0-9]+[MIDNSHP=X]')

# Function to replace/insert str in str
def str_mod(string, pos, insert, replace=False):
    """
    This function facilitates insertion/replacement of a string into another
    one at a certain position

    Args:
        string (:obj:`str`): String INTO which the insert will be put
        pos (int): Position (0-based) of string where the insert will be put.
            NOTE, negative values will produce unintended results and are not
            checked!
        insert (:obj:`str`): String that will be inserted
        replace (:obj:`bool`, optinal): Whether or not the insertion should
            replace the characters in the original string
    Returns:
        :obj:`str` The new string with the insertion
    """
    next_pos = pos + len(insert) * replace
    return string[:pos] + insert + string[next_pos:]

# Function to parse a CIGAR
def parse_cigar(cigar):
    """
    Parses a CIGAR string into a list

    Args:
        cigar (:obj:`str`): The CIGAR string. IMPORTANT: Any parts of the CIGAR
            that do not confirm the SAM specifications will be silently ignored
    Returns:
        List[:obj:`str`] where each item matches r'[0-9]+[MIDNSHP=X]'
    """
    return CIGAR_REGEX.findall(cigar)

# Function to convert CIGAR to a list of operations
def cigar_to_list(cigar):
    """
    Parses a CIGAR string into a list of tuples, where each tuple describes a
    single CIGAR operation

    Args:
        cigar (:obj:`str`): The CIGAR string. IMPORTANT: Any parts of the CIGAR
            that do not confirm the SAM specifications will be silently ignored
    Returns:
        List[Tuple(:obj:`str`, int)] where each item describes a single
            operation: ('operant', size in nt)
    """
    cigar_ops = []
    for operant in parse_cigar(cigar):
        cigar_ops.append((operant[-1], int(operant[:-1])))
    return cigar_ops

# Function to calculate the length of SEQ based on cigar
def seq_len_from_cigar(cigar):
    """
    Calculates the sequence length from the CIGAR string

    Args:
        cigar (:obj:`str`): The CIGAR string. IMPORTANT: Any parts of the CIGAR
            that do not confirm the SAM specifications will be silently ignored
    Returns:
        int length of sequence in nt
    """
    seq_len = 0
    for operant, op_len in cigar_to_list(cigar):
        if operant in ('M', 'I', 'S', '=', 'X'):
            seq_len += op_len
    return seq_len

def aln_len_from_cigar(cigar):
    """
    Calculates the alignment length from the CIGAR string

    Args:
        cigar (:obj:`str`): The CIGAR string. IMPORTANT: Any parts of the CIGAR
            that do not confirm the SAM specifications will be silently ignored
    Returns:
        int length of alingment in bp
    """
    aln_len = 0
    for operant, op_len in cigar_to_list(cigar):
        if operant in ('M', 'D', 'N', 'P'):
            aln_len += op_len
    return aln_len

# Function to parse an MD tag (not fully tested)
def parse_md(md_tag):
    """
    Parses a MD string into a list

    Args:
        md_tag (:obj:`str`): The MD string. IMPORTANT: Any parts of the MD
            that do not confirm the SAM specifications will be silently ignored
    Returns:
        List[:obj:`str`] where each item matches r'([0-9]+|[A-Z]+|\^[A-Z]+)'
    """
    return MD_REGEX.findall(md_tag)

# Function to confirm a seq, CIGAR, MD tag triple confirm each other
def assert_aln(read_seq, cigar, md_tag):
    """
    Args:
        read_seq (:obj:`str`): The (query) sequence used in alignment
        cigar (:obj:`str`): The CIGAR string. IMPORTANT: Any parts of the CIGAR
            that do not confirm the SAM specifications will be silently ignored
        md_tag (:obj:`str`): The MD string. IMPORTANT: Any parts of the MD
            that do not confirm the SAM specifications will be silently ignored
    Returns:
        Tuple(None, None, None) if assertion fails
        Tuple(parsed_cigar, parsed_md, mod_cigar) if assertion succeeds
    """
    is_ok = True
    parsed_md = parse_md(md_tag)
    parsed_cigar = cigar_to_list(cigar)
    md_len = 0
    md_dels = []
    for operant in parsed_md:
        if operant.isdigit():
            md_len += int(operant)
        elif operant[0] == '^':
            md_dels.append((md_len, len(operant)-1))
        else:
            md_len += len(operant)
    cigar_len = 0
    cigar_len2 = 0
    mod_cigar = []
    for operant, op_len in parsed_cigar:
        if operant == 'D':
            cigar_del = (cigar_len2, op_len)
            if cigar_del not in md_dels:
                is_ok = False
            else:
                md_dels.remove(cigar_del)
        elif operant in ('I', 'S'):
            cigar_len += op_len
        elif operant not in ('H', 'N'):
            cigar_len += op_len
            cigar_len2 += op_len
        if operant == 'I':
            mod_cigar.append((cigar_len2, operant, op_len))
    # print("Seq len: {}\nCigar Len: {}\nMD len: {}\nCigar Len2: {}\nis_ok: {}".format(
    #     len(read_seq), cigar_len, md_len, cigar_len2, is_ok))
    if md_dels or not (cigar_len2 == md_len and cigar_len == len(read_seq)):
        is_ok = False
    if is_ok:
        return (parsed_cigar, parsed_md, mod_cigar)
    return (None, None, None)

def edit_dist_range(cigar, md_tag):
    """
    Calculates the edit distance by position from a pair of CIGAR and MD

    Args:
        cigar (:obj:`str`): The CIGAR string. IMPORTANT: Any parts of the CIGAR
            that do not confirm the SAM specifications will be silently ignored
        md_tag (:obj:`str`): The MD string. IMPORTANT: Any parts of the MD
            that do not confirm the SAM specifications will be silently ignored
    Returns:
        List(int) where each item denotes the edit distance at its index. The
            edit distance is calculated cumulatively from the 5' end where it
            starts with 0. With each edit, it increases by 1.
    """
    md_by_cigar = []
    md_by_pos = []
    parsed_cigar = cigar_to_list(cigar)
    parsed_md = parse_md(md_tag)
    cur_edit = 0
    for operant in parsed_md:
        if operant.isdigit():
            md_by_pos += [cur_edit] * int(operant)
        elif operant[0] == '^':
            md_by_pos += range(cur_edit+1, cur_edit+len(operant))
            cur_edit = md_by_pos[-1]
        else:
            md_by_pos += range(cur_edit+1, cur_edit+len(operant)+1)
            cur_edit = md_by_pos[-1]
    cur_edit = 0
    for operant, op_len in parsed_cigar:
        if operant in ('S', 'H', 'P'):
            pass
        elif operant in ('M', 'D'):
            md_by_cigar += [cur_edit] * op_len
        elif operant == 'I':
            cur_edit += op_len
        elif operant == 'N':
            md_by_pos = (md_by_pos[:len(md_by_cigar)]
                         + [md_by_pos[len(md_by_cigar)-1]] * op_len
                         + md_by_pos[len(md_by_cigar):])
            md_by_cigar += [cur_edit] * op_len
    if not len(md_by_cigar) == len(md_by_pos):
        raise Exception(
            "Length of edit distance arrays by cigar and md do not match")
    edit_range = [md_by_cigar[i] + md_by_pos[i] for i
                  in range(len(md_by_cigar))]
    return edit_range


# Function to convert a CIGAR and MD tag to a gapped alignment
def sam2gapped(read_seq, cigar, md_tag, skip_mode="skip"):
    """
    This function takes the read sequence, CIGAR and MD tag
    as inputs and produces a gapped alignment of the read sequence
    against the reference sequence without the need of reference
    sequence as an input:

    ref     AGTGCCTTGGGTGTTCA-----ATCCCCATGCAACAACC
    aln     || ||||    ||||||     | |||||||||||||||
    read    AGAGCCTC---TGTTCACATAGACCCCCATGCAACAACC

    """
    # Assert and parse the CIGAR and MD tag
    cigar_ops, md_ops, mod_cigar = assert_aln(read_seq, cigar, md_tag)
    if not cigar_ops:
        raise Exception('Problem: {},{},{}'.format(read_seq, cigar, md_tag))

    # Initialize the ref seq and the pairwise aligment
    ref_seq = read_seq
    pairwise = '|' * len(read_seq)

    # First loop over the MD tag
    # Skip soft clipped bps (S) at the 5'
    soft_pos = 0
    for i in (0, min(1, len(cigar_ops) - 1)):
        if cigar_ops[i][0] not in ('H', 'S'):
            break
        elif cigar_ops[i][0] == 'S':
            soft_pos += cigar_ops[i][1]
    cur_pos = soft_pos
    for md_op in md_ops:
        while mod_cigar and mod_cigar[0][0] <= cur_pos - soft_pos:
            insertion = mod_cigar.pop(0)
            cur_pos += insertion[2]
        # a match
        if md_op.isdigit():
            op_len = int(md_op)
        # a deletion
        elif md_op[0] == '^':
            op_len = len(md_op) - 1
            read_seq = str_mod(read_seq, cur_pos, '-' * op_len)
            ref_seq = str_mod(ref_seq, cur_pos, md_op[1:])
            pairwise = str_mod(pairwise, cur_pos, ' ' * op_len)
        # a mismatch
        else:
            op_len = len(md_op)
            ref_seq = str_mod(ref_seq, cur_pos, md_op, replace=True)
            pairwise = str_mod(pairwise, cur_pos, '.' * op_len, replace=True)
        cur_pos += op_len

    # Then loop over CIGAR
    cur_pos = 0
    for cigar_op, op_len in cigar_ops:
        if cigar_op in ('I', 'S'):
            ref_seq = str_mod(ref_seq, cur_pos, '-'*op_len, replace=True)
            pairwise = str_mod(pairwise, cur_pos, ' '*op_len, replace=True)
        if cigar_op == 'N':
            if skip_mode == 'full':
                insert = 'N' * op_len
            elif skip_mode == 'skip':
                insert = ''
                op_len = 0
            elif skip_mode == 'short':
                insert = '[' + str(op_len) + 'N]'
                op_len = len(insert)
            else:
                raise Exception("skip_mode not known")
            ref_seq = str_mod(ref_seq, cur_pos, insert)
            pairwise = str_mod(pairwise, cur_pos, ' ' * op_len)
            read_seq = str_mod(read_seq, cur_pos, '-' * op_len)
        if not cigar_op in ('H', 'P'):
            cur_pos += op_len
    return '\n'.join(['ref     ' + ref_seq,
                      'aln     ' + pairwise,
                      'read    ' + read_seq])
