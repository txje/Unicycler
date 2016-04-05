from __future__ import print_function
from __future__ import division

import subprocess
import sys
import os
import re
from assembly_graph import Segment

class SegmentWithAlignments(object):
    '''
    This class holds a segment and its long read alignments.
    '''
    def __init__(self, segment, long_reads_fastq, graphmap_path, working_dir, allowed_end_clips):
        '''
        The constructor actually conduct the GraphMap alignment of the long reads to the segment.
        '''
        self.segment = segment
        self.allowed_end_clips = allowed_end_clips
        temp_segment_fasta = os.path.join(working_dir, 'NODE_' + str(segment.number) + '.fasta')
        temp_segment_sam = os.path.join(working_dir, 'NODE_' + str(segment.number) + '.sam')
        segment.save_to_fasta(temp_segment_fasta)
        command = [graphmap_path, '-r', temp_segment_fasta, '-d', long_reads_fastq, '-o',
                   temp_segment_sam, '-Z']
        process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        _, _ = process.communicate()
        self.alignments = load_alignments(temp_segment_sam, segment)
        self.alignments = self.filter_bad_alignments(self.alignments)
        os.remove(temp_segment_fasta)
        os.remove(temp_segment_fasta + '.gmidx')
        os.remove(temp_segment_fasta + '.gmidxsec')
        # os.remove(temp_segment_sam)

    def filter_bad_alignments(self, alignments):
        '''
        Removes alignments with unaligned parts. Specifically, it looks for cases where the start
        of the alignment has soft clipping but the alignment did not begin at the start of the
        segment, or where the end of the alignment has soft clipping but the alignment did not end
        at the end of the segment.
        '''
        filtered_alignments = []
        for alignment in alignments:
            bad_alignment = False
            segment_start = alignment.segment_start
            segment_unaligned_end = alignment.segment_unaligned_end
            if alignment.reverse_complement:
                segment_start, segment_unaligned_end = segment_unaligned_end, segment_start
            if alignment.get_start_soft_clips() > self.allowed_end_clips and \
               segment_start > self.allowed_end_clips:
                bad_alignment = True
            elif alignment.get_end_soft_clips() > self.allowed_end_clips and \
               segment_unaligned_end > self.allowed_end_clips:
                bad_alignment = True
            if not bad_alignment:
                filtered_alignments.append(alignment)
        return filtered_alignments

    def get_contained_alignments(self):
        '''
        Returns all of the alignments where the read is entirely contained within the contig. I.e.
        there is no significant clipping on either end.
        '''
        return [x for x in self.alignments \
                if x.get_start_soft_clips() < self.allowed_end_clips and \
                x.get_end_soft_clips() < self.allowed_end_clips]
        




class LongRead(object):
    '''
    This class holds a long read, e.g. from PacBio or Oxford Nanopore.
    '''
    def __init__(self, name, sequence, qualities):
        self.name = name
        self.sequence = sequence
        self.qualities = qualities
        self.alignments = []



class Alignment(object):
    '''
    This class describes an alignment between a long read and a contig.
    '''
    def __init__(self, sam_line, segment):
        sam_parts = sam_line.split('\t')
        self.read_name = sam_parts[0]
        self.flag = int(sam_parts[1])
        self.reverse_complement = bool(self.flag & 0x10)
        self.segment_name = sam_parts[2]
        self.mapping_quality = int(sam_parts[4])
        self.cigar = sam_parts[5]
        self.cigar_parts = re.findall(r'\d+\w', self.cigar)
        self.read_sequence = sam_parts[9]
        self.read_length = len(self.read_sequence)
        self.read_quality = sam_parts[10]
        self.flags = sam_parts[11:]
        self.edit_distance = None
        self.alignment_score = None
        self.e_value = None
        self.segment = segment
        for flag in self.flags:
            if flag.startswith('NM:i:'):
                self.edit_distance = int(flag[5:])
            if flag.startswith('AS:i:'):
                self.alignment_score = int(flag[5:])
            if flag.startswith('ZE:f:'):
                self.e_value = float(flag[5:])
        self.segment_start = int(sam_parts[3])
        self.segment_end = self.segment_start - 1
        for cigar_part in self.cigar_parts:
            self.segment_end += get_reference_shift_from_cigar_part(cigar_part)
        self.segment_unaligned_end = segment.get_length() - self.segment_end
        self.read_start_fraction = self.get_start_soft_clips() / self.read_length
        self.read_end_fraction = (self.read_length - self.get_end_soft_clips()) / self.read_length
        self.segment_start_fraction = (self.segment_start - 1) / segment.get_length()
        self.segment_end_fraction = self.segment_end / segment.get_length()

    def __repr__(self):
        return self.read_name + '_' + self.segment_name

    def get_start_soft_clips(self):
        '''
        Returns the number of soft-clipped bases at the start of the alignment.
        '''
        match = re.search(r'^\d+S', self.cigar)
        if not match:
            return 0
        else:
            return int(match.group(0)[:-1])

    def get_end_soft_clips(self):
        '''
        Returns the number of soft-clipped bases at the start of the alignment.
        '''
        match = re.search(r'\d+S$', self.cigar)
        if not match:
            return 0
        else:
            return int(match.group(0)[:-1])


def load_long_reads(fastq_filename):
    '''
    This function loads in long reads from a FASTQ file and returns a dictionary where key = read
    name and value = LongRead object.
    '''
    reads = {}
    fastq = open(fastq_filename, 'r')
    for line in fastq:
        name = line.strip()[1:]
        sequence = next(fastq).strip()
        _ = next(fastq)
        qualities = next(fastq).strip()
        reads[name] = LongRead(name, sequence, qualities)
    fastq.close()
    return reads


def load_alignments(sam_filename, segment):
    '''
    This function returns a list of Alignment objects from the given SAM file.
    '''
    alignments = []
    sam_file = open(sam_filename, 'r')
    for line in sam_file:
        if line.startswith('@'):
            continue
        line = line.strip()
        if line.split('\t')[2] != '*':
            alignments.append(Alignment(line, segment))
    return alignments

def get_reference_shift_from_cigar_part(cigar_part):
    '''
    This function returns how much a given cigar moves on a reference.
    Examples:
      * '5M' returns 5
      * '5S' returns 0
      * '5D' returns 5
      * '5I' returns 0
    '''
    if cigar_part[-1] == 'M':
        return int(cigar_part[:-1])
    if cigar_part[-1] == 'D':
        return int(cigar_part[:-1])
    if cigar_part[-1] == 'S':
        return 0
    if cigar_part[-1] == 'I':
        return 0


