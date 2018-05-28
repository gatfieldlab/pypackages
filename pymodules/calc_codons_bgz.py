#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 17 17:23:09 2013

@author: barpat
"""

import sys, re, operator
import argparse
import numpy as np
import get_density_bgz as gd
from collections import defaultdict
import codecs


__author__ = "Bulak Arpat"
__copyright__ = "Copyright 2014-2017, Bulak Arpat"
__license__ = "GPLv3"
__version__ = "1.2.0"
__maintainer__ = "Bulak Arpat"
__email__ = "Bulak.Arpat@unil.ch"
__status__ = "Development"


np.seterr(invalid='raise')

REGIONS = ('5utr', 'cds', '3utr', 'transcript')
SUMMARY_OPERATIONS = ('mean', 'sum', 'median')
CDS_START = 0
CDS_END = 1
TR_LEN = 2
MAX_READ_LEN = 100
#ASITE_OFF = ('5p', {i: 15+int(i>30) for i in range(1,MAX_READ_LEN+1)})
ASITE_OFF = ('5p', {i: 15 for i in range(1,MAX_READ_LEN+1)})
#ASITE_OFF = ('3p', {i: 15 for i in range(1,MAX_READ_LEN+1)})
OK_FLAGS = ['OKEY', 'SNGL']
REGION_PATTERN = r'^((?P<region>{region})|(?P<number>\d+)|(?P<star>(?P<pre>\*)?(?=({region}))({region})(?(pre)|\*?))(?P<op>[\+\-])(?P<offset>\d+))$'.format(region='|'.join(REGIONS))
REGION_MATCH = re.compile(REGION_PATTERN)
OPS = {'+': operator.add, '-': operator.sub}

class densityError(Exception):
    def __init__(self, m):
        self.message = m

class regionError(Exception):
    def __init__(self, m):
        self.message = m

class Region(object):
    '''
    [ar]{,1}\d+
    cds = [cds_start - cds_end]
    cds:3utr = [cds_start - 3utr_end]
    *cds = [cds_start - cds_start]
    *cds:cds* = [cds_start - cds_end]
    *cds-50:cds = [cds_start-50 - cds_end]
    *cds-50:cds*+50 = [cds_start-50 - cds_end+50]
    '''
    def __init__(self, region):
        self.parse_words = []
        parsed_regions = region.split(":")
        for parsed_region in parsed_regions:
            match = REGION_MATCH.search(parsed_region)
            if match:
                if match.group('region'):
                    self.parse_words.append([match.group('region')])
                elif match.group('number'):
                    self.parse_words.append([match.group('number')])
                else:
                    self.parse_words.append([match.group('star'),
                                             match.group('op'),
                                             match.group('offset')])
            else:
                raise regionError("'{}' could not be match to region regex".format(parsed_region))

    def get_limits(self, tr_info):
        limits = []
        for sentence in self.parse_words:
            limits.append([])
            for word in sentence:
                if word.isdigit():
                    limits[-1].append(int(word))
                elif word in ['+','-']:
                    limits[-1].append(word)
                else:
                    left, right = word[0] == '*', word[-1] == '*'
                    if left:
                        word = word[1:]
                    if right:
                        word = word[:-1]
                    if word == '5utr':
                        word_limits = (0, tr_info[CDS_START])
                    elif word == 'cds':
                        word_limits = (tr_info[CDS_START], tr_info[CDS_END])
                    elif word == '3utr':
                        word_limits = (tr_info[CDS_END], tr_info[TR_LEN])
                    elif word == 'transcript':
                        word_limits = (0, tr_info[TR_LEN])
                    if left:
                        limits[-1].append(word_limits[0])
                    if right:
                        limits[-1].append(word_limits[1])
                    if not (left or right):
                        limits[-1].extend(word_limits)
        result = []
        for limit in limits:
            if len(limit) == 1 and not (type(limit[0]) is int):
                raise regionException('Fatal error: {}'.format(limits))
            while '-' in limit or '+' in limit:
                o_i, op = next((i,o) for i,o in enumerate(limit) if o == '-' or o == '+')
                try:
                    res = OPS[op](limit[o_i-1], limit[o_i+1])
                except KeyError:
                    raise regionException('Fatal error: {}'.format(limits))
                limit = limit[:max(0, o_i-1)] + [res] + limit[min(len(limit), o_i+2):]
            result.append(limit)
        lmin = min([min(l) for l in result])
        lmax = max([max(l) for l in result])
        return (lmin, lmax)

class DensityParser(object):
    '''
    A class
    '''

    def __init__(self, tr_db, tr_ids, density_files, region):
        self.tr_db = tr_db
        self.tr_ids = tr_ids
        self.density_files = []
        for density_file in density_files:
            try:
                self.density_files.append(gd.get_index_tr(density_file, self.tr_ids))
            except gd.fileError as e:
                raise densityError(e.message)

        self.region = Region(region)

    def get_densitybin(self, asite_off=ASITE_OFF, pileup=False, filter_size=0):
        if asite_off[0] not in ['5p', '3p']:
            raise densityError("Supplied A-site offset rule is not correct: {}\n".format(asite_off))
        while True:
            extracts = []
            is_eof = []
            for densities in self.density_files:
                try:
                    extracts.append(next(densities))
                except gd.fileError as e:
                    raise densityError(e.message)
                except StopIteration:
                    is_eof.append(True)
                else:
                    is_eof.append(False)
            if all(is_eof):
                raise StopIteration
            elif any(is_eof):
                raise densityError("Some density files reached EOF before others.\n")
            else:
                tr_ids, extracts, errors = list(zip(*extracts))
                is_error = [error is not None for error in errors]
                trs_same = len(set(tr_ids)) == 1
                if not trs_same:
                    raise densityError("Some density files return TR-IDs out-of-sync.\n")
                if is_error:
                    pass # for now
                cur_err = []
                result = []
                tr_id = tr_ids[0]
                gid = []
                try:
                    tr_info = self.tr_db[tr_id]
                except KeyError:
                    yield (None, tr_id, None, None, [densityError("<{}> could not be found in DB\n".format(tr_id))])
                else:
                    slice_limits = self.region.get_limits(tr_info)
                    if slice_limits[0] < 0 or slice_limits[1] > tr_info[TR_LEN]:
                        yield (None, tr_id, None, None, [densityError("Slice {} can't be applied to <{}>\n".format(slice_limits, tr_id))])
                        continue
                    for i, extract in enumerate(extracts):
                        density = []
                        cur_gene = ''
                        if isinstance(errors[i], gd.indexError):
                            cur_err.append(densityError(errors[i].message))
                        elif extract == '':
                            cur_err.append(densityError("<{}> returned empty list from density file\n".format(tr_id)))
                        else:
                            ids = extract.split('\t',1)[0]
                            cur_gene, cur_tr = ids.split('|')
                            if not cur_tr == tr_id:
                                cur_err.append(densityError("TR-ID from density extract does not match the returned TR-ID\n"))
                            else:
                                for line in extract.strip().split('\n'):
                                    parsed = line.split('\t')
                                    try:
                                        read_pos = int(parsed[1])
                                        read_len = int(parsed[2])
                                    except IndexError:
                                        cur_err.append(densityError("Can't parse the density data returned by get_density from '{}'!\n".format(density_file)))
                                        break
                                    except ValueError:
                                        cur_err.append(densityError("Can't interpret density data returned by get_density!\n"))
                                        break
                                    if filter_size and read_len != filter_size:
                                        continue
                                    if pileup:
                                        density.extend([i for i in range(read_pos, read_pos + read_len) if i < tr_info[TR_LEN]])
                                    else:
                                        try:
                                            offset = asite_off[1][read_len]
                                        except KeyError:
                                            continue
                                            # sys.stderr.write("Found a read with a size outside of specified range of 1 to {}\nWill continue setting size to {}!\n".format(MAX_READ_LEN, MAX_READ_LEN))
                                            # offset = asite_off[1][MAX_READ_LEN]
                                        if asite_off[0] == '5p':
                                            asite = read_pos + offset
                                        elif asite_off[0] == '3p':
                                            asite = read_pos + read_len - offset
                                        density.append(asite)
                            cur_err.append(None)
                        if density:
                            binned_density=np.bincount(density, minlength=tr_info[TR_LEN])
                        else:
                            binned_density=np.zeros(tr_info[TR_LEN])
                        result.append(binned_density[slice_limits[0]:slice_limits[1]])
                        if cur_gene:
                            gid.append(cur_gene)
                    assert len(result) == len(cur_err)
                    gid = set(gid)
                    if not gid:
                        gid.add('')
                    assert len(gid) == 1, "{} : {} : {}".format(tr_id, gid, cur_gene)
                    yield (gid.pop(), tr_id, slice_limits, result, cur_err)

def summarize_densitybins(densitybins, operation='sum', **kwargs):
    if operation in SUMMARY_OPERATIONS:
        try:
            op = getattr(np, operation)
        except AttributeError:
            sys.stderr.write("{} is not a numpy attribute??".operation)
            sys.exit(1)
        try:
            res = op(np.array(densitybins), axis=0, **kwargs)
        except TypeError as e:
            sys.stderr.write("Numpy error during '{}' operation: {}".format(operation, str(e)))
            sys.exit(1)
        else:
            return res

def normalize_densitybins(densitybins, norm_factors):
    if norm_factors:
        return np.array(densitybins) / np.array(norm_factors)
    else:
        return np.array(densitybins)

def aggregate_codons(densitybin, operation='sum', **kwargs):
    if operation in SUMMARY_OPERATIONS:
        try:
            op = getattr(np, operation)
        except AttributeError:
            sys.stderr.write("{} is not a numpy attribute??".operation)
            sys.exit(1)
        try:
            num_codons = np.ceil(np.divide(len(densitybin),3.0))
            codons = np.array_split(np.array(densitybin), num_codons)
        except TypeError as e:
            sys.stderr.write("Numpy error during '{}' operation: {}".format(operation, str(e)))
            sys.exit(1)
        aggregated = []
        for codon in codons:
            aggregated.append(op(codon, axis=0, **kwargs))
        return np.array(aggregated)

def unescaped_str(arg_str):
    return codecs.decode(str(arg_str), 'unicode_escape')

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('files', nargs='+', metavar='density_file', type=str)
    parser.add_argument('-t', '--tr_ids', type=str, default='tr_ids.txt')
    parser.add_argument('-r', '--region', type=str)
    parser.add_argument('-c', '--cds', type=str, default='prepared_cds.txt')
    parser.add_argument('-l', '--max_read_len', type=int, default=MAX_READ_LEN,
                        help='specify max length of reads, default={}'.format(MAX_READ_LEN))
    parser.add_argument('-d', '--delimiter', type=unescaped_str, default='\t',
                        help='delimiter used in output table')
    parser.add_argument('-n', '--no_skip', action='store_true',
                        help='Include failed TR-IDs in the output with an NA')
    report = parser.add_mutually_exclusive_group(required=False)
    report.add_argument('-p', '--pileup', help='report pileup counts per nt',
                        action='store_true')
    report.add_argument('-q', '--quantify', action='store_true',
                        help='report only total read count in the region')
    report.add_argument('-f', '--fraction', action='store_true',
                        help='report fraction of reads to all reads in region')
    args = parser.parse_args()
    args.delimiter = args.delimiter
    tr_ids = []
    try:
        with open(args.tr_ids) as tr_file:
            for line in tr_file:
                tr_ids.append(line.strip())
    except IOError:
        sys.stderr.write('Could not read the TR-IDs from {}\n'.format(args.tr_ids))
        sys.exit(1)

    tr_db = {}
    try:
        with open(args.cds) as cds_file:
            for line in cds_file:
                parsed = line.strip().split('\t')
                if parsed[2] in tr_ids:
                    (tr_len, cds_start, cds_end) = (int(parsed[i]) for i in (5,6,7))
                    tr_db[parsed[2]] = (cds_start, cds_end, tr_len)
    except IOError:
        sys.stderr.write('Could not read the CDS models from {}\n'.format(args.cds))
        sys.exit(1)

    density_extractor = DensityParser(tr_db, tr_ids, args.files, args.region)
    for gid, tr_id, slice_limits, extracts, errors in density_extractor.get_densitybin(pileup=args.pileup):
        if len(errors) == 1 and isinstance(errors[0], densityError):
            sys.stderr.write(errors[0].message)
            if args.no_skip:
                print(tr_id, 'NA', sep=args.delimiter)
        else:
            summarized = summarize_densitybins(extracts)
            total_count = np.sum(summarized)
            if args.fraction:
                fractions = summarized / max(total_count, 1)
                for i, pos in enumerate(summarized):
                    print(tr_id, i+1, pos, fractions[i], sep=args.delimiter)
            elif args.quantify:
                print(tr_id, total_count, sep=args.delimiter)
            else:
                for i, pos in enumerate(summarized):
                    print(tr_id, i+1, pos, sep=args.delimiter)


if __name__ == "__main__":
    main()
#         densitybins = get_densitybin(tr_ids, density_file)
#         while True:
#             try:
#                 densitybin = next(densitybins)
#                 if isinstance(densitybin, densityError):
#                     raise densitybin
#             except StopIteration:
#                 break
#             except densityError as e:
#                 sys.stderr.write(e.message)
#             else:
#                 cur_tr = densitybin[0][1]
#                 bp = 0
#                 read_sum = np.sum(densitybin[1])
#                 if read_sum == 0:
#                     read_sum = 1
#                 for count in densitybin[1]:
#                     bp += 1
#                     codon = (bp - 1) // 3 + 1
#                     frame = (bp - 1) % 3 + 1
#                     try:
#                         print(cur_tr, bp, codon, frame, count, count/read_sum, sep='\t')
#                     except FloatingPointError:
#                         sys.stderr.write("{}\n".format('\t'.join(['FlotaingPointError', cur_tr, str(codon), str(frame), str(count), str(read_sum)])))
#                 #print(tr_db[cur_tr]) # debug
