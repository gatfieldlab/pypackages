# -*- coding: utf-8 -*-
"""
Created on Sat Oct 27 16:43:25 2012

@author: bulak
"""

import operator
import HTSeq
from accessories import count_utils
#
# GFF single locus
#
class gff_single_locus:

    accepted_types = {'all': set(['protein_coding', 'lincRNA']),
                      'lincRNA': set(['processed_transcript', 'lincRNA', 'antisense'])}
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
        
    def reduce_exons(self):
        #exon_map = set([])
        if not self.biotype in self.accepted_types:
            cbiotype = 'all'
        else:
            cbiotype = self.biotype
        exon_types = set([e.source for e in self.features['exon']])
        if len(exon_types) > 1:
            filtered_types = exon_types & self.accepted_types[cbiotype]
        else:
            filtered_types = exon_types
        if not filtered_types:
            print 'I have got nothing', self.id
        #print exon_types, "=>", filtered_types
        #for exon in self.features['exon']:
            #pass
            #if exon.source in accepted_types:
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
    def get_gene(self):
        self.features['gene'].append(self.get_feature(self.iv.start, self.iv.end, 1, ftype='gene', source=None))
    def merge_small_exons(self, cutoff=100):
        count = 0
        for exonic_part in range(1,len(self.features['exonic_part'])-1):
            common_left_and_right_transcripts = set()
            cur_exons = self.features['exonic_part'][exonic_part-1:exonic_part+2]
            if cur_exons[1].iv.length > cutoff:
                continue
            for e in (range(exonic_part)+range(exonic_part+1,len(self.features['exonic_part']))):
                common_left_and_right_transcripts |= set(self.features['exonic_part'][e].attr['transcripts'].split('+'))
            is_overlapping_with_intron = common_left_and_right_transcripts - set(cur_exons[1].attr['transcripts'].split('+'))
            consecutive = False
            #is_overlapping_with_intron = (set(cur_exons[0].attr['transcripts'].split('+')) & set(cur_exons[2].attr['transcripts'].split('+'))) - set(cur_exons[1].attr['transcripts'].split('+'))
            #print is_overlapping_with_intron
            # joinable 
            # consider some transcripts that extend introns for two or more exons in the other transcripts
            # does the above method work?
            # how to join - join small ones first?
            # perhaps loop first, identify joinable ones (non intron overlapping)
            # then join until > 100?? also how to keep track of this?
            # so that feature list is updated??
            # exon_parts need to be renumbered...
            if cur_exons[1].iv.start == cur_exons[0].iv.end:
                consecutive = True
#                print "I am left consecutive", cur_exons[0].iv, cur_exons[1].iv
            if cur_exons[1].iv.end == cur_exons[2].iv.start:
                consecutive = True
#                print "I am right consecutive", cur_exons[1].iv, cur_exons[2].iv
            if bool(len(is_overlapping_with_intron)) and consecutive:
                count += 1
        return count
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


            
            


