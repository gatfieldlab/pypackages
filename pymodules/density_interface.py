#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""A lightwight module for interfacing a CLI app with density profiles
providing consistent argument parser, certain checks and data structures
"""

import sys
import re
import argparse
import calc_codons_bgz as cc
from extended import argparse_types


__author__ = "Bulak Arpat"
__copyright__ = "Copyright 2016, 2017, Bulak Arpat"
__license__ = "GPLv3"
__version__ = "0.8.0"
__maintainer__ = "Bulak Arpat"
__email__ = "Bulak.Arpat@unil.ch"
__status__ = "Development"


class InterfaceException(Exception):
    """Exceptions for density_interface"""
    def __init__(self, m):
        message = "[density_interface] {}".format(m)
        super(InterfaceException, self).__init__(message)


# disome project specific settings
# ideally these kind of variables should be
# read via a local database
ZTS = [0, 2, 12]
READ_TYPES = ('Mono', 'Di', 'Tr')
READ_LEN = {'Mono': {'min': 26, 'max': 35},
            'Di': {'min': 45, 'max': 70},
            'Tr': {'min': 21, 'max': 60}}
READ_RANGE = {k: range(v['min'], v['max']+1) for k, v in READ_LEN.items()}
ASITE_OFF = {
    'Mono': {'std': ('5p', {i: 15+int(i > 30) for i in READ_RANGE['Mono']}),
             '3p': ('3p', {i: 0 for i in READ_RANGE['Mono']})},
    'Di': {'selective': ('5p', {59: 45, 60: 45, 62: 48, 63: 48}),
           'permissive': ('5p', {57: 45, 58: 45, 59: 45, 60: 45,
                                 61: 46, 62: 48, 63: 48, 64: 49}),
           'very_permissive': ('5p', {
               i: int(i >= READ_LEN['Di']['min'] and i <= 57)*(i-15) +
               int(i == 58)*44 +
               int(i > 58 and i < 61)*45 +
               int(i == 61)*46 +
               int(i > 61 and i <= 63)*48 +
               int(i == 64)*49 +
               int(i > 64 and i <= 67)*51 for i in
               READ_RANGE['Di'] if i not in [54, 55, 68, 69, 70]}),
           '3p': ('3p', {i: 0 for i in READ_RANGE['Di']})},
    'Tr': {'std': ('5p', {i: i // 2 for i in READ_RANGE['Tr']}),
           '3p': ('3p', {i: 0 for i in READ_RANGE['Tr']})}}
ASITE_KEYS = tuple(set([k for t in READ_TYPES for k in ASITE_OFF[t]]))
DEFAULTS = {'Mono': 'std', 'Di': 'selective', 'Tr': 'std'}
DB_SAMPLE, DB_REP, DB_FOLDER, DB_FILE, DB_EXT, DB_FLAG, DB_NORM = range(1, 8)
RP_TR = r'({})'.format('|'.join(READ_TYPES))
ZT = r'{}(\d+)'.format(RP_TR)


aggregate_codons = cc.aggregate_codons


def get_argparser(prog_name, descr, exclude=None,
                  parser_class=argparse.ArgumentParser, max_help_pos=27):
    """Sets up an argparse based CLI argument, help, usage utility and parser_class
    the arguments. Certain arguments can be left out via exlude"""

    def custom_formatter():
        """A custom formatter lambda function to set the max_help_pos of the
        argparse.ArgumentDefaultsHelpFormatter with a user set value"""
        return lambda prog: argparse.ArgumentDefaultsHelpFormatter(
            prog, max_help_position=max_help_pos)

    if not exclude:
        exclude = []
    elif isinstance(exclude, str):
        exclude = [exclude]
    elif not isinstance(exclude, (list, tuple)):
        raise InterfaceException(
            "get_argparser's exclude argument has to be either None, str, or a"
            " a list/tuple (of strings)")
    parser = parser_class(prog=prog_name, description=descr,
                          formatter_class=custom_formatter())
    parser.add_argument(
        '-s', '--offset', help='offset to A-site. Possible options depends on'
        ' the type of read. Accepted values include positive integers and one'
        ' of {}'.format(('default',)+ASITE_KEYS), type=str, default='default')
    parser.add_argument(
        '-l', '--lengths', type=argparse_types.mix_int_list, default='28:32',
        help='length(s) of footprints included, for example 28 is 28 only; '
        '29:31 is [29-31]; 28,30:32 is [28,30-32]')
    parser.add_argument(
        '-d', '--density-db-file', type=str, required=True, help='specially '
        'formatted text file indexing the density files for all samples')
    parser.add_argument(
        '-i', '--tr-list', type=str, help='a file containing the list of '
        'transcripts to include in the computation. Note the final list of '
        'transcripts will be an intersection of this list and what is '
        'available in the CDS file')
    parser.add_argument(
        '-f', '--flag', action='store_true', help='consider flags when parsing'
        ' CDS file')
    if 'normalize' not in exclude:
        parser.add_argument(
            '-n', '--normalize', action='store_true', help='normalize the '
            'counts?')
    if 'summarize' not in exclude:
        parser.add_argument(
            '-u', '--summarize', type=str, default='sum', help='summarize the '
            'counts over files, and if so how?',
            choices=('',) + cc.SUMMARY_OPERATIONS)
    region_parser = parser.add_mutually_exclusive_group()
    region_parser.add_argument(
        '-m', '--margins', metavar='L,R', help="Margins: L and R nts from 5' "
        "and 3' termini, respectively, will be excluded from any calculations",
        type=argparse_types.capped_tuple(mins=(0, 0), maxs=(240, 240)),
        default=(0, 0))
    region_parser.add_argument(
        '-g', '--region', type=str, help='A string variable defining the '
        'region of the transcript that will be used for computations. Examples'
        ' include: "cds": CDS region, "cds:3utr": CDS+3UTR, "*cds:cds*": from '
        'start to end of CDS, same as "cds", "*cds+20:3utr*-50": from CDS-star'
        't+20 nt to 3UTR-end-50 nt, "50:150": from nt-50 to nt-150 of the TR')
    parser.add_argument(
        '-c', '--cds', type=str, help="A CDS file of type 'prepared_cds'. This"
        " file is necessary if CDS info can not be read from the FASTA file, "
        "and the info read here will override any CDS information extracted "
        "from the FASTA file")
    parser.add_argument(
        '-z', '--zt', type=argparse_types.int_list, default=[ZTS[0]],
        help='comma separated ZT values.')
    if 'type' not in exclude:
        parser.add_argument(
            '-t', '--type', type=str, choices=READ_TYPES,
            help='Type of read to process.', default='Mono')
    parser.add_argument(
        '-r', '--rep', type=argparse_types.int_list, default='1', help='comma'
        ' separated REP values. Possible values are 1 1,2 2 1,3 etc depending '
        'on replicate numbers.')
    return parser


def get_offset(offset_str, read_type=None, read_lengths=None):
    """Checks and returns the correct A-site offsetting values for the given
    read type and offset argument. Restricts to given read lengths if provided
    by read_lengths"""
    def _fixed_offset(i_val):
        return {read_len: i_val for read_len in READ_RANGE[read_type]}
    # overwrite the read_type if one is supplied through offset_str
    try:
        read_type, offset_str = offset_str.split(':', 1)
    except ValueError:
        if read_type is None:
            raise InterfaceException(
                'When read type is not set, offset should be given in the form'
                ' read_type:offset')
    if offset_str == 'default':
        offset = ASITE_OFF[read_type][DEFAULTS[read_type]]
    elif offset_str in ASITE_OFF[read_type]:
        offset = ASITE_OFF[read_type][offset_str]
    elif offset_str.isdigit():
        offset = ('5p', _fixed_offset(int(offset_str)))
    else:
        raise InterfaceException(
            "Requested offset '{}' is not defined for read type '{}'\n".format(
                offset_str, read_type))
    if read_lengths:
        offset = (offset[0], {read_len: offset[1][read_len] for read_len in
                              read_lengths if read_len in offset[1]})
        if len(offset[1]) == 0:
            raise InterfaceException(
                "Read length restriction results in empty offset structure")
    return offset


def get_db_files(db_file, zts, read_type, reps):
    """Reads and parses the given db flat file and extracts file paths
    corresponding to the supplied read type, ZT and replicate specs

    Args:
        db_file: path to the DB file, should be openable by open()
        zts: an iterable :obj: of integers for time-points
        read_type: :obj:`str` for read type
        reps: an iterable :obj: of integers for replicate numbers
    Returns:
        density_files: :obj:`list` of path string to density files
        norm_factors: :obj:`list` of float normalization factor of each file
    """
    zt_regex = re.compile(ZT)
    rp_regex = re.compile(RP_TR)
    density_files = []
    norm_factors = []
    try:
        with open(db_file) as db_handle:
            for line in db_handle:
                parsed = line.strip().split('\t')
                if parsed[DB_FLAG] == 'Y':
                    cur_sample = parsed[DB_SAMPLE]
                    rp = rp_regex.search(cur_sample).group(1)
                    if not rp == read_type:
                        continue
                    zt = int(zt_regex.search(cur_sample).group(2))
                    if zt not in zts:
                        continue
                    rep = int(parsed[DB_REP])
                    if rep not in reps:
                        continue
                    cur_file = "{}{}{}".format(parsed[DB_FOLDER],
                                               parsed[DB_FILE],
                                               parsed[DB_EXT])
                    density_files.append(cur_file)
                    norm_factors.append([float(parsed[DB_NORM])])
    except IOError:
        raise InterfaceException(
            'Could not read the database from {}\n'.format(db_file))
    except (AttributeError, TypeError, IndexError) as e:
        raise InterfaceException(
            'Could not parse the database: {}\nError: {}\n'.format(db_file, e)
        )
    return density_files, norm_factors


def process_cds(cds_file, use_flag=False):
    """Parses a specially formatted file containing CDS information. The file
    has to of 'prepared_cds' format"""
    tr_db = {}
    p_index = 6, 7, 5
    t_index = 2,
    try:
        with open(cds_file) as cds_handle:
            for line in cds_handle:
                parsed = line.strip().split('\t')
                tr_id = '|'.join([parsed[i] for i in t_index])
                if parsed[3] == 'composite' or (use_flag and parsed[8] != '*'):
                    continue
                tr_db[tr_id] = tuple([int(parsed[i]) for i in p_index])
    except IOError:
        raise InterfaceException(
            "Could not read the CDS info: {}\n".format(cds_file))
    except IndexError:
        raise InterfaceException(
            "Could not parse the CDS info: {}\n".format(cds_file))
    return tr_db


def get_trs(tr_file):
    """Parses (any white space) and extracts TR-ids from a flat file"""
    if not tr_file:
        return None
    tr_list = []
    with open(tr_file) as tr_f:
        for line in tr_f:
            tr_list.extend(line.strip().split())
    return tr_list


class DensityInterface(object):
    """A class for interfacing density file database and accessing density
    data directly via an iterator which uses calc_codons"""

    def __init__(self, density_db_file, region, zt, rtype, rep,
                 cds_file, asite_off, tr_list_file=None, flag=False):
        self.meta = {'density_db': density_db_file, 'region': region,
                     'zt': zt, 'type': rtype, 'rep': rep,
                     'cds_file': cds_file, 'offset': asite_off,
                     'tr_list': tr_list_file, 'flag': flag}
        self.region = region
        self.offset = asite_off
        self.files, self.norm_factors = get_db_files(density_db_file, zt,
                                                     rtype, rep)
        self.tr_db = process_cds(cds_file, flag)
        tr_list = get_trs(tr_list_file)
        if tr_list is not None:
            self.transcripts = list(set(self.tr_db.keys()) & set(tr_list))
        else:
            self.transcripts = list(self.tr_db.keys())
        self.success = 0
        self.fail = 0
        self.source = cc.DensityParser(self.tr_db, self.transcripts,
                                       self.files, self.region)

    def iterate(self, normalize=True, operation='sum'):
        """Iterates over transcripts in the selected density files. To extract
        density information, calc_codons is used. Data can be normalized and
        summarized (mean, median, sum, ...?)"""
        for (gid, tr_id, limits, extracts,
             errors) in self.source.get_densitybin(asite_off=self.offset):
            at_least_one_pass = False
            for error in errors:
                if not error:
                    at_least_one_pass = True
            if at_least_one_pass:
                seq_id = "|".join([gid, tr_id])
                norm_factors = normalize and self.norm_factors
                normalized = cc.normalize_densitybins(extracts, norm_factors)
                if operation:
                    summarized = cc.summarize_densitybins(normalized,
                                                          operation)
                else:
                    summarized = normalized
                self.success += 1
                yield (seq_id, summarized, limits)
            else:
                self.fail += 1
                yield (tr_id, None, None)

        def get_meta(self):
            """Updates and returns the metadata of the class isinstance"""
            return self.meta.update({'files': self.files, '#fail': self.fail,
                                     '#tr': len(self.transcripts),
                                     '#success': self.success})


def main():
    """Reads CLI arguments and simply outputs the processed densities per
    transcript"""
    my_parser = get_argparser('density_interface.py', __doc__)
    args = my_parser.parse_args()
    if not args.region:
        region = "*cds+{}:cds*-{}".format(*args.margins)
    else:
        region = args.region
    offset = get_offset(args.type, args.offset)
    my_den = DensityInterface(args.density_db_file, region, args.zt,
                              args.type, args.rep, args.cds, offset,
                              args.tr_list, args.flag)
    for seq_id, density in my_den.iterate(args.normalize, args.summarize):
        if density is not None:
            print(seq_id, density)


if __name__ == "__main__":
    sys.exit(main())
