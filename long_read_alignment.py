from __future__ import print_function
from __future__ import division

import subprocess
import sys
import os
import re
import datetime
from assembly_graph import Segment
from assembly_graph import reverse_complement



class LongRead(object):
    '''
    This class holds a long read, e.g. from PacBio or Oxford Nanopore.
    '''
    def __init__(self, name, sequence, qualities):
        self.name = name
        self.sequence = sequence
        self.qualities = qualities
        self.alignments = []

    def add_alignment(self, new_alignment):
        self.alignments.append(new_alignment)
        self.alignments = sorted(self.alignments, key=lambda x: x.read_start_pos)



class Alignment(object):
    '''
    This class describes an alignment between a long read and a contig.
    '''
    def __init__(self, sam_line, segment_lengths):

        # Load all important parts from the SAM line.
        sam_parts = sam_line.split('\t')
        self.read_name = sam_parts[0].split('/')[0]
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
        for flag in self.flags:
            if flag.startswith('NM:i:'):
                self.edit_distance = int(flag[5:])
            if flag.startswith('AS:i:'):
                self.alignment_score = int(flag[5:])
            if flag.startswith('ZE:f:'):
                self.e_value = float(flag[5:])

        # Determine the position of the alignment in the segment.
        # Positions are stored with a Python style 0-based index.
        self.segment_length = segment_lengths[self.segment_name]
        self.segment_start_pos = int(sam_parts[3]) - 1
        self.segment_end_pos = self.segment_start_pos
        for cigar_part in self.cigar_parts:
            self.segment_end_pos += get_reference_shift_from_cigar_part(cigar_part)

        # If the match is to the reverse complement, flip the sequence and alignment.
        if self.reverse_complement:
            self.read_sequence = reverse_complement(self.read_sequence)
            self.read_quality = self.read_quality[::-1]
            self.cigar_parts = self.cigar_parts[::-1]
            self.cigar = ''.join(self.cigar_parts)
            old_start = self.segment_start_pos
            self.segment_start_pos = self.segment_length - self.segment_end_pos
            self.segment_end_pos = self.segment_length - old_start

        self.segment_end_gap = self.segment_length - self.segment_end_pos

        # Determine the position of the alignment in the read.
        self.read_start_pos = self.get_start_soft_clips()
        self.read_end_pos = self.read_length - self.get_end_soft_clips()
        self.read_end_gap = self.get_end_soft_clips()

    def __repr__(self):
        if self.reverse_complement:
            strand = '-'
        else:
            strand = '+'
        return self.read_name + ' (' + str(self.read_start_pos) + '-' + str(self.read_end_pos) + \
               '), ' + self.segment_name + ' (' + str(self.segment_start_pos) + '-' + \
               str(self.segment_end_pos) + ', strand: ' + strand + ')'

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

def run_blasr_alignment_all_segments(graph, long_reads_fasta, sam_file, blasr_path, working_dir):
    '''
    Runs BLASR to produce a SAM file of alignments.
    '''
    graph_fasta = os.path.join(working_dir, 'graph.fasta')
    graph.save_to_fasta(graph_fasta)
    run_blasr(graph_fasta, long_reads_fasta, sam_file, blasr_path, working_dir)
    os.remove(graph_fasta)

def run_graphmap_alignment_all_segments(graph, long_reads_fastq, sam_file, graphmap_path, working_dir):
    '''
    Runs GraphMap to produce a SAM file of alignments.
    '''
    graph_fasta = os.path.join(working_dir, 'graph.fasta')
    graph.save_to_fasta(graph_fasta)
    run_graphmap(graph_fasta, long_reads_fastq, sam_file, graphmap_path, working_dir)
    os.remove(graph_fasta)

def run_graphmap_owler(graph, long_reads_fastq, sam_file, graphmap_path, working_dir):
    graph_fasta = os.path.join(working_dir, 'graph.fasta')
    graph.save_to_fasta(graph_fasta)
    command = [graphmap_path, '-r', graph_fasta, '-d', long_reads_fastq, '-o',
               sam_file, '-w', 'owler', '-L', 'paf', '-Z']
    process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    _, _ = process.communicate()
    os.remove(graph_fasta + '.gmidxowl')
    os.remove(graph_fasta)

def run_graphmap_alignment_one_segment_at_a_time(graph, long_reads_fastq, sam_file,
                                                 graphmap_path, working_dir):
    final_sam = open(sam_file, 'w')

    for segment in graph.segments.values():
        segment_short_name = 'NODE_' + str(segment.number)
        segment_fasta = os.path.join(working_dir, segment_short_name + '.fasta')
        segment.save_to_fasta(segment_fasta)
        segment_sam_filename = os.path.join(working_dir, segment_short_name + '.sam')
        run_graphmap(segment_fasta, long_reads_fastq, segment_sam_filename, graphmap_path,
                     working_dir)

        # Copy the segment's SAM alignments to the final SAM file.
        segment_sam = open(segment_sam_filename, 'r')
        for line in segment_sam:
            if not line.startswith('@') and line.split('\t', 3)[2] != '*':
                final_sam.write(line)
        segment_sam.close()

        # Clean up
        os.remove(segment_fasta)
        os.remove(segment_sam_filename)

    final_sam.close()


def run_graphmap(fasta, long_reads_fastq, sam_file, graphmap_path, working_dir):

    # First build the index. This may not be necessary in the future, but a development version of
    # GraphMap I was using could crash if you didn't do this first.
    command = [graphmap_path, '-I', '-r', fasta]
    process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    _, _ = process.communicate()

    # Now actually run GraphMap.
    command = [graphmap_path, '-r', fasta, '-d', long_reads_fastq, '-o',
               sam_file, '-Z', '-F', '1.0', '-t', '8']
    process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    _, _ = process.communicate()

    # Clean up.
    os.remove(fasta + '.gmidx')
    os.remove(fasta + '.gmidxsec')




def run_blasr(fasta, long_reads_fasta, sam_file, blasr_path, working_dir):

    my_env = os.environ.copy()
    my_env['DYLD_LIBRARY_PATH'] = '/Users/Ryan/Applications/blasr_install/blasr/libcpp/alignment:/Users/Ryan/Applications/blasr_install/blasr/libcpp/hdf:/Users/Ryan/Applications/blasr_install/blasr/libcpp/pbdata:/Users/Ryan/Applications/blasr_install/hdf5/lib'
    command = [blasr_path, long_reads_fasta, fasta, '-out', sam_file, '-sam', '-nproc', '8', '-affineAlign']
    process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, env=my_env)
    _, _ = process.communicate()

def load_alignments(sam_filename, segment_lengths):
    '''
    This function returns a list of Alignment objects from the given SAM file.
    '''
    alignments = []
    sam_file = open(sam_filename, 'r')
    for line in sam_file:
        if not line.startswith('@') and line.split('\t', 3)[2] != '*':
            alignments.append(Alignment(line, segment_lengths))
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

def filter_by_end_gaps(alignments, allowed_gaps):
    '''
    Removes alignments with unaligned parts. Specifically, it looks for cases where the start
    of the alignment has soft clipping but the alignment did not begin at the start of the
    segment, or where the end of the alignment has soft clipping but the alignment did not end
    at the end of the segment.
    '''
    good_alignments = []
    bad_alignments = []

    for alignment in alignments:
        s_start_gap = alignment.segment_start_pos
        s_end_gap = alignment.segment_end_gap
        r_start_gap = alignment.read_start_pos
        r_end_gap = alignment.read_end_gap

        if (s_start_gap > allowed_gaps and r_start_gap > allowed_gaps) or \
           (s_end_gap > allowed_gaps and r_end_gap > allowed_gaps):
            bad_alignments.append(alignment)
        else:
            good_alignments.append(alignment)

    return good_alignments, bad_alignments

