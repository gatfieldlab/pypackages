# -*- coding: utf-8 -*-

"""
Based on gtf_utils, was created on Sat Oct 27 16:43:25 2012
Created on Fri May 29 15:26:25 2015
Now writing pretty much from scratch

@author: bulak
"""

import warnings
import operator
import HTSeq
import sys
from collections import defaultdict
from operator import attrgetter

class GtfUtilsError(Exception):

    """Base class for gtf_utils exceptions."""

class GtfUtilsWarning(Warning):

    """Base class for gtf_utils warnings."""


class SpecialGenomicFeature(HTSeq.GenomicFeature):
    direction_repr = {'+': '->', '-': '<-', '.': '<>'}
    def __repr__(self):
        return "<{}: '{}' at {}: {} {} {} (strand '{}')>".format(
        self.__class__.__name__, self.name,
        self.iv.chrom, self.iv.start_d, self.direction_repr[self.iv.strand], self.iv.end_d, self.iv.strand)

#
# Exon
#
class Exon(SpecialGenomicFeature):
    accepted_sources = ("BED", "GTF")
    def __init__(self, genomic_feature, source_type):
        SpecialGenomicFeature.__init__(self, str(genomic_feature.iv), "exon", genomic_feature.iv)
        self.type = 'exon'
        if source_type == 'GTF':
            self.attr = {'gene_id': genomic_feature.attr['gene_id']}
        elif source_type == 'BED':
            self.attr = {}

class Intron(SpecialGenomicFeature):
    def __init__(self, genomic_feature, transcript_id = None):
        self.type = 'intron'
        self.attr = getattr(genomic_feature, 'attr',
                                {'transcript_id': transcript_id or ""})
        if transcript_id:
            self.attr['transcript_id'] = transcript_id
        transcript_id = self.attr['transcript_id']
        if not transcript_id:
            raise GtfUtilsError("An Intron object has to have a transcript_id")
        SpecialGenomicFeature.__init__(self, genomic_feature.name, self.type, genomic_feature.iv)

class Transcript(SpecialGenomicFeature):
    def __init__(self, genomic_feature, transcript_id = None, gene_id = None):
        self.attr = getattr(genomic_feature, 'attr', {'gene_id': gene_id or "", 'transcript_id': transcript_id or ""})
        if transcript_id:
            self.attr['transcript_id'] = transcript_id
        if gene_id:
            self.attr['gene_id'] = gene_id
        transcript_id = self.attr['transcript_id']
        gene_id = self.attr['gene_id']
        if transcript_id == gene_id:
            raise GtfUtilsError("A Transcript object has to have a transcript_id that is different than its gene_id")
        if not gene_id:
            raise GtfUtilsError("A Transcript object has to have a gene_id")
        SpecialGenomicFeature.__init__(self, transcript_id, 'transcript', genomic_feature.iv)
        self.type = 'transcript'
        self.exons = []
        self.introns = []
        self.gene_id = gene_id
    def __repr__(self):
        return "<{}: '{}' from '{}' at {}: {} -> {} (strand '%s')>".format(
        self.__class__.__name__, self.name, self.gene_id,
        self.iv.chrom, self.iv.start_d, self.iv.end_d, self.iv.strand)
    def add_exon(self, exon):
        if isinstance(exon, Exon):
            self.exons.append(exon)
        return len(self.exons)
    def sorted_exons(self):
        return sorted(self.exons, key=attrgetter('iv.start'))

class Gene(SpecialGenomicFeature):
    f_types = ('transcripts', 'exons', 'flat_exons', 'introns')
    export_ftypes = ('flat_exons', 'introns')
    def __init__(self, genomic_feature, gene_id = None):
        self.attr = getattr(genomic_feature, 'attr', {'gene_id': gene_id or ""})
        if gene_id:
            self.attr['gene_id'] = gene_id
        gene_id = self.attr['gene_id']
        if not gene_id:
            raise GtfUtilsError("A Gene object has to have a gene_id")
        SpecialGenomicFeature.__init__(self, gene_id, 'gene', genomic_feature.iv)
        self.type = 'gene'
        self.transcripts = {}
        self.exons = {}
        self.flat_exons = []
        self.introns = []
    def add_transcript(self, transcript):
        if isinstance(transcript, Transcript):
            self.transcripts[transcript.name] = transcript
    def add_exon(self, exon, transcript_id=None):
        if isinstance(exon, Exon):
            if not exon.name in self.exons:
                self.exons[exon.name] = exon
            else:
                exon = self.exons[exon.name]
            if not transcript_id:
                warnings.warn("Added exon without transcript binding", GtfUtilsWarning)
            else:
                transcript = self.get_transcript(transcript_id)
                if transcript:
                    transcript.add_exon(exon)
                else:
                    raise GtfUtilsError("Transcript ID '{}' does not exists in '{}': can not associate exon with this transcript".format(transcript_id, self.name))
    def get_transcript(self, transcript_id):
        if transcript_id in self.transcripts:
            return self.transcripts[transcript_id]
        else:
            return None
    def transcript_count(self):
        return len(self.transcripts)
    def exon_count(self):
        return len(self.exons)
    def sorted_exons(self):
        return sorted(self.exons, key=attrgetter('iv.start'))
    def print_transcripts(self):
        print self.name
        for transcript in self.transcripts.values():
            print transcript.name
            print " : ".join([str(exon.iv) for exon in transcript.sorted_exons()])
    def extend_termini(self, relative_to_gene=False):
        # the idea is to extend the 5' and 3' ends of those exons
        # that are not surpassed by an exon-exon junction
        # [===]--- ok
        #   [=]--- ok
        #   [====]-ok
        #      [=]-NOT ok
        #
        # Idea: we can check the size of difference and apply a filter there
        exons_5 = sorted(self.exons.values(), key=attrgetter('iv.end'))
        exons_3 = sorted(self.exons.values(), key=attrgetter('iv.start'))
        p5_terminus = exons_3[0].iv.end
        p3_terminus = exons_5[-1].iv.start
        if relative_to_gene:
            left_border = self.iv.start
            right_border = self.iv.end
        else:
            left_border = exons_3[0].iv.start
            right_border = exons_5[-1].iv.end
        for exon in self.exons.values():
            if exon.iv.start <= p5_terminus:
                exon.iv.start = left_border
            if exon.iv.end >= p3_terminus:
                exon.iv.end = right_border
    def reduce_transcripts(self):
        if self.transcript_count() > 1:
            cur_trs = self.transcripts.values()
            nonoverlapping_set = []
            while cur_trs:
                cur_tr = cur_trs.pop()
                cur_exons = tuple(sorted([(exon.iv.start, exon.iv.end) for exon in cur_tr.exons]))
                if not nonoverlapping_set:
                    nonoverlapping_set.append((cur_exons,[cur_tr.name]))
                    #start = (cur_exons,[cur_tr.name])
                else:
                    added = False
                    for index, (exons, tr_ids) in enumerate(nonoverlapping_set):
                        a = b = None
                        is_a_superset = False
                        #if exons[0][0] == cur_exons[0][0] and exons[-1][1] == cur_exons[-1][1]:
                        if len(exons) >= len(cur_exons):
                            a = exons
                            b = cur_exons
                        else:
                            a = cur_exons
                            b = exons
                        if a:
                            is_a_superset = True
                            b_in_a = []
                            for exon in b:
                                try:
                                    b_in_a.append(a.index(exon))
                                except ValueError:
                                    is_a_superset = False
                                    break
                            if is_a_superset:
                                bsum = (b_in_a[-1]*(b_in_a[-1]+1) - b_in_a[0]*(b_in_a[0]-1))/2
                                if not sum(b_in_a) == bsum:
                                    is_a_superset = False
                        if is_a_superset:
                            tr_ids.append(cur_tr.name)
                            nonoverlapping_set[index] = (a, tr_ids)
                            added = True
                            break
                    if not added:
                        nonoverlapping_set.append((cur_exons,[cur_tr.name]))
            if len(nonoverlapping_set) == 1:
                gene_type = 'MULTIPLE_INCLUSIVE'
            else:
                gene_type = 'MULTIPLE_NONINCLUSIVE'
        else:
            gene_type = 'SINGLE_TRANSCRIPT'
        print gene_type
        print nonoverlapping_set
    def flatten_exons(self, call_introns=True):
        exon_starts, exon_ends = zip(*[(exon.iv.start, exon.iv.end) for exon in self.exons.values()])
        exon_starts = sorted(list(exon_starts))
        exon_ends = sorted(list(exon_ends))
        cur_junctions = []
        split_exons = []
        while exon_ends:
            added = False
            if exon_starts:
                if exon_starts[0] < exon_ends[0]:
                    cur_junctions.append(exon_starts.pop(0))
                    added = True
            if not added:
                cur_junctions.append(exon_ends.pop(0))
            if len(exon_starts) == len(exon_ends):
                cur_junctions = sorted(set(cur_junctions))
                for n, junc in enumerate(cur_junctions[:-1]):
                    split_exons.append((junc, cur_junctions[n+1]))
                cur_junctions = []
        for n, (start, end) in enumerate(split_exons):
            cur_exonic_part_number = "{:03d}".format(n+1)
            genomic_feature = HTSeq.GenomicFeature("{}#{}".format(self.name, cur_exonic_part_number),
                                                    "exonic_part",
                                                    HTSeq.GenomicInterval(self.iv.chrom, start, end, self.iv.strand))
            self.flat_exons.append(Exon(genomic_feature, 'BED'))
            self.flat_exons[-1].attr = {"exon_number": cur_exonic_part_number, "gene_id": self.name,
                                        "transcript_id": self.name+'_flatten',
                                        #"transcripts": '+'.join(sorted(overlapping_attr['transcript'])),
                                        "gene_name": self.attr['gene_name'], "gene_biotype": self.attr['gene_biotype']}
            self.flat_exons[-1].source = "gtf_utils2.py"
        if call_introns:
            cur_intronic_part = 0
            for n, (start, end) in enumerate(split_exons[:-1]):
                if split_exons[n+1][0] > end:
                    cur_intronic_part += 1
                    cur_intronic_part_number = "{:03d}".format(cur_intronic_part)
                    genomic_feature = HTSeq.GenomicFeature("{}#{}".format(self.name, cur_intronic_part_number),
                                                           "intronic_part",
                                                           HTSeq.GenomicInterval(self.iv.chrom, end, split_exons[n+1][0],
                                                           self.iv.strand))
                    genomic_feature.attr = {"intron_number": cur_intronic_part_number, "gene_id": self.name,
                                            "transcript_id": self.name+'_flatten', "gene_name": self.attr['gene_name'],
                                            "gene_biotype": self.attr['gene_biotype']}
                    self.introns.append(Intron(genomic_feature, self.name+'_flatten'))
                    self.introns[-1].source = "gtf_utils2.py"
    def export_features(self, export_ftypes=None):
        export_ftypes = export_ftypes or self.export_ftypes
        fexport = {}
        for f_type in export_ftypes:
            fexport[f_type] = []
            features = getattr(self, f_type)
            if isinstance(features, dict):
                features = features.values()
            for feature in features:
                fexport[f_type].append(dict(feature.attr.items() +
                                            {'name': feature.name, 'chrom': feature.iv.chrom,
                                             'start': feature.iv.start, 'end': feature.iv.end,
                                             'strand': feature.iv.strand,
                                             'source': getattr(feature, 'source', "N/A")}.items()
                                             ))
        return fexport
    def export_data(self, export_ftypes=None):
        return {'id': self.name, 'chrom':self.iv.chrom, 'start': self.iv.start, 'end': self.iv.end,
                'strand': self.iv.strand, 'gff_line': self.get_gff_line(),
                'features': self.export_features(export_ftypes), 'f_types': self.f_types}
    def get_features_gff(self, f_type='flat_exons'):
        features = getattr(self, f_type, [])
        if isinstance(features, dict):
            features = features.values()
        gff_line = ""
        for feature in features:
            gff_line += feature.get_gff_line()
        return gff_line

class Genome(object):
    def __init__(self, name):
        self.name = name
        self.genomic_array = HTSeq.GenomicArrayOfSets('auto')
        self.genes = {}
    def add_gene(self, gene):
        if isinstance(gene, Gene):
            self.genes[gene.name] = gene
            self.genomic_array[gene.iv] += gene.name
    def get_chroms(self):
        return self.genomic_array.chrom_vectors
    def get_size(self):
        return len(self.genes)
    def get_gene(self, gene_id):
        if gene_id in self.genes:
            return self.genes[gene_id]
        else:
            return None
    def get_array(self, iv):
        if isinstance(iv, HTSeq.GenomicInterval):
            return self.genomic_array[iv].steps()
        else:
            warnings.warn("Genome.get_array only accepts a HTSeq.GenomicInterval object as argument", GtfUtilsWarning)
            return None
    def add_exon(self, feature, source ):
        if isinstance(feature, HTSeq.GenomicFeature):
            try:
                attr = feature.attr
            except AttributeError:
                raise GtfUtilsError("Can not add an exon to a genome without GTF attributes")
            try:
                gene_id = attr['gene_id']
            except KeyError:
                raise GtfUtilsError("Can not add an exon to a genome without geneID")
            cur_gene = self.get_gene(gene_id)
            if not cur_gene:
                raise GtfUtilsError("Exon with GeneID '{}' can not be attached to any genes in genome".format(gene_id))
            transcript_id = attr.get('transcript_id')
            cur_gene.add_exon(Exon(feature, source), transcript_id)
    def add_transcript(self, transcript):
        if isinstance(transcript, Transcript):
            try:
                gene_id = transcript.attr['gene_id']
            except KeyError:
                raise GtfUtilsError("Can not add a transcript to a genome without geneID")
            try:
                cur_gene = self.genes[gene_id]
            except KeyError:
                raise GtfUtilsError("Transcript with TranscriptID '{}' can not be attaced to any genes".format(transcript.name))
            cur_gene.add_transcript(transcript)



#
# GFF s locus
#
class Locus(object):

    def __init__(self, id):
        self.id = id
        self.iv = None
        self.gff_line = None
        self.biotype = None
        self.name = None
        self.transcripts = set()
        self.features = {'CDS':[], 'gene':[], 'exon':[], 'intron':[], 'reduced_exon':[], 'intronic_part':[], 'exonic_aggregate':[], 'exonic_part':[]}
        self.f_types = self.features.keys()
        self.export_ftypes = ['exonic_part', 'intronic_part']

    def extend_termini(self):
        # the idea is to extend the 5' and 3' ends of those exons
        # that are not surpassed by an exon-exon junction
        # [===]--- ok
        #   [=]--- ok
        #   [====]-ok
        #      [=]-NOT ok
        exons_5 = sorted(self.features['exon'], key=lambda x: x.iv.end)
        exons_3 = sorted(self.features['exon'], key=lambda x: x.iv.start)
        p5_terminus = exons_3[0].iv.end
        p3_terminus = exons_5[-1].iv.start
        for exon in self.features['exon']:
            if exon.iv.start <= p5_terminus:
                exon.iv.start = self.iv.start
            if exon.iv.end >= p3_terminus:
                exon.iv.end = self.iv.end
    def set_iv(self, base='exon'):
        if not self.iv:
            start = float("inf")
            end = 0
            for exon in self.features[base]:
                if exon.iv.start < start:
                    start = exon.iv.start
                if exon.iv.end > end:
                    end = exon.iv.end
            self.iv = HTSeq.GenomicInterval(exon.iv.chrom, start, end, exon.iv.strand)
        else:
            warnings.warn("Locus IV already set. Will not set it again.", GtfUtilsWarning)
    def merge_transcripts(self):
        pass
    def get_gene(self):
        self.features['gene'].append(self.get_feature(self.iv.start, self.iv.end, 1, ftype='gene', source=None))
    def get_gff_line(self):
        return self.gff_line
    def add_feature(self, feature):
        if feature.type in self.f_types:
            self.features[feature.type].append(feature)
    def export_features(self):
        fexport = {}
        for f_type in self.export_ftypes:
            fexport[f_type] = []
            for feature in self.features[f_type]:
                fexport[f_type].append(dict({'name': feature.name, 'chrom': feature.iv.chrom,
                                             'start': feature.iv.start, 'end': feature.iv.end,
                                             'strand': feature.iv.strand,
                                             'source': feature.source}.items() + feature.attr.items()
                                            ))
        return fexport
    def export_data(self):
        return {'id': self.id, 'start': self.iv.start, 'end': self.iv.end,
                'strand': self.iv.strand, 'gff_line': self.gff_line,
                'features': self.export_features(),
                'f_types': self.f_types}
    def flatten_exons(self, keep_original=True):
        split_exons = set([])
        exon_quicksect = count_utils.create_quicksect({self.id: self}, ['exon'])
        for exon in self.features['exon']:
            cur_iv = exon.iv
            overlapping_exons = count_utils.get_features_within_interval(cur_iv.chrom, cur_iv.start, cur_iv.end-1, exon_quicksect)
            overlapping_exons.sort(key=operator.itemgetter(2))
            left_border = cur_iv.start
            right_border = cur_iv.end
            split_sites = set([left_border, right_border])
            for overlapping_exon in overlapping_exons:
                if int(overlapping_exon[2]) > left_border:
                    split_sites.add(int(overlapping_exon[2]))
                if int(overlapping_exon[3]) < right_border:
                    split_sites.add(int(overlapping_exon[3])+1)
            split_sites = list(split_sites)
            split_sites.sort()
            #print split_sites
            for i in range(len(split_sites)-1):
                #print split_sites[i], split_sites[i+1]
                # this is in python coordinates [0,4) = 1-4
                split_exons.add((exon.iv.chrom, split_sites[i], split_sites[i+1]))
        split_exons = sorted(split_exons)
        if not keep_original:
            self.features['exon'] = []
            ftype = 'exon'
        else:
            ftype = 'exonic_part'
        cur_exonic_part_number = 0
        for (chrom, start, end) in split_exons:
            overlapping_exons = count_utils.get_features_within_interval(chrom, start, end-1, exon_quicksect)
            overlapping_attr = {'transcript': set([]), 'source': set([])}
            for exon in overlapping_exons:
                overlapping_attr['transcript'].add(exon[6]['transcript_id'])
                overlapping_attr['source'].add(exon[4])
            cur_exonic_part_number += 1
            cur_exonic_part = HTSeq.GenomicFeature(self.id+"#"+'%03d'%cur_exonic_part_number, ftype,
                                             HTSeq.GenomicInterval(chrom, start, end, self.iv.strand))
            cur_exonic_part.attr = {"{}_number".format(ftype): '%03d'%cur_exonic_part_number, "gene_id": self.id,
                                    "transcript_id": self.id+'_flatten',
                                    "transcripts": '+'.join(sorted(overlapping_attr['transcript'])),
                                    "gene_name": self.name, "gene_biotype": self.biotype}
            cur_exonic_part.source = "+".join(sorted(overlapping_attr['source']))
            #print cur_exonic_part.get_gff_line()
            self.add_feature(cur_exonic_part)
    def get_feature(self, start, end, num, ftype='intronic_part', source='intronic_part'):
        if not source:
            source = ftype
        f = HTSeq.GenomicFeature(self.id, ftype,
                                 HTSeq.GenomicInterval(self.iv.chrom, start, end, self.iv.strand))
        f.__setattr__('attr', {'gene_id':self.id, ftype+'_number': '%03d'%num,
                               'transcript_id': self.id+'_aggtr',
                               'gene_name':self.name, 'gene_biotype':self.biotype})
        f.__setattr__('source', source)
        return f
    def amalgamate(self, source='amalgamate_exon'):
        one_exon = self.get_feature(self.iv.start, self.iv.end, 1, ftype='exon', source=source)
        self.features['exon'] = [one_exon]
    def call_introns(self):
        if len(self.features['exonic_part']) == 0:
            self.flatten_exons()
        if len(self.features['exonic_part']) < 2:
            return
        intron_num = 0
        locus_start = self.iv.start
        locus_end =   self.iv.end
        for exon in self.features['exonic_part']:
            if exon.iv.start == locus_start:
                cur_intronic_start = exon.iv.end
            else:
                cur_intronic_end = exon.iv.start
                if cur_intronic_end - cur_intronic_start > 0:
                    intron_num += 1
                    self.features['intronic_part'].append(self.get_feature(cur_intronic_start, cur_intronic_end, intron_num))
                if not exon.iv.end == locus_end:
                    cur_intronic_start = exon.iv.end
    def call_introns2(self, add_exons=False, source=None):
        """ This function does not need flat/split exons """
        exons = [(exon.iv.start, exon.iv.end) for exon in self.features['exon']]
        exons.sort()
        start = 0
        end = 1
        fused_exons = []
        while exons:
            if len(exons) == 1 or exons[0][end] < exons[1][start]:
                fused_exons.append(exons.pop(0))
            else:
                if exons[0][end] < exons[1][end]:
                    exons = [(exons[0][start], exons[1][end])] + exons[2:]
                else:
                    exons.pop(1)
        intron_num = 0
        if add_exons:
            self.features['exonic_aggregate'].append(self.get_feature(fused_exons[0][start], fused_exons[0][end], 1, 'exonic_aggregate', source))
        if len(fused_exons) < 2:
            return
        for i in range(1, len(fused_exons)):
            if add_exons:
                self.features['exonic_aggregate'].append(self.get_feature(fused_exons[i][start], fused_exons[i][end], i+1, 'exonic_aggregate', source))
            intron_num += 1
            self.features['intronic_part'].append(self.get_feature(fused_exons[i-1][end], fused_exons[i][start], intron_num, source=source))
