from __future__ import print_function
from __future__ import division

import subprocess
import sys
import os
import re
import datetime
from assembly_graph import Segment
from assembly_graph import reverse_complement
from assembly_graph import AssemblyGraph



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
    def __init__(self, sam_line, references):

        # Load all important parts from the SAM line.
        sam_parts = sam_line.split('\t')
        self.read_name = sam_parts[0].split('/')[0]
        self.flag = int(sam_parts[1])
        self.reverse_complement = bool(self.flag & 0x10)
        self.reference_name = sam_parts[2]
        self.mapping_quality = int(sam_parts[4])
        self.cigar = sam_parts[5]
        self.cigar_parts = re.findall(r'\d+\w', self.cigar)
        self.full_read_sequence = sam_parts[9]
        self.read_length = len(self.full_read_sequence)
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

        # Determine the position of the alignment in the reference.
        self.reference_length = len(references[self.reference_name])
        self.reference_start_pos = int(sam_parts[3]) - 1
        self.reference_end_pos = self.reference_start_pos
        for cigar_part in self.cigar_parts:
            self.reference_end_pos += get_reference_shift_from_cigar_part(cigar_part)
        self.reference_end_gap = self.reference_length - self.reference_end_pos

        # Determine the position of the alignment in the read.
        self.read_start_pos = self.get_start_soft_clips()
        self.read_end_pos = self.read_length - self.get_end_soft_clips()
        self.read_end_gap = self.get_end_soft_clips()

        # Extend the alignment so it is fully semi-global, reaching the end of the sequence.
        self.extend_alignment()

        # Count matches, mismatches, insertions and deletions.
        # Insertions and deletions are counted per base. E.g. 5M3I4M has 3 insertions, not 1.
        self.matches = 0
        self.mismatches = 0
        self.insertions = 0
        self.deletions = 0
        self.percent_identity = 0.0
        self.tally_up_alignment(references)

    def __repr__(self):
        if self.reverse_complement:
            strand = '-'
        else:
            strand = '+'
        return self.read_name + ' (' + str(self.read_start_pos) + '-' + str(self.read_end_pos) + \
               '), ' + self.reference_name + ' (' + str(self.reference_start_pos) + '-' + \
               str(self.reference_end_pos) + ', strand: ' + strand + '), ' + \
               '%.2f' % self.percent_identity + '%'

    def get_alignment_length_read(self):
        '''
        Returns the length of the aligned read sequence.
        '''
        return self.read_end_pos - self.read_start_pos

    def get_alignment_length_reference(self):
        '''
        Returns the length of the aligned reference sequence.
        '''
        return self.reference_end_pos - self.reference_start_pos

    def extend_alignment(self):
        '''
        This function extends the alignment as much as possible in both directions so the alignment
        only terminates when it reaches the end of either the read or the reference.
        It does not actually perform the alignment - it just counts each alignment as a match. This
        means that very long extensions will probably result in terrible alignments, but that's
        okay because we'll filter alignments by quality later.
        '''
        missing_bases_at_start = min(self.read_start_pos, self.reference_start_pos)
        missing_bases_at_end = min(self.read_end_gap, self.reference_end_gap)

        if missing_bases_at_start:
            # Adjust the start of the reference.
            self.reference_start_pos -= missing_bases_at_start

            # Adjust the start of the read and fix up the CIGAR to match.
            self.read_start_pos -= missing_bases_at_start
            self.cigar_parts.pop(0)
            if self.cigar_parts[0][-1] == 'M':
                new_match_length = missing_bases_at_start + int(self.cigar_parts[0][:-1])
                self.cigar_parts.pop(0)
                new_cigar_part = str(new_match_length) + 'M'
            else:
                new_cigar_part = str(missing_bases_at_start) + 'M'
            self.cigar_parts.insert(0, new_cigar_part)
            if self.read_start_pos > 0:
                self.cigar_parts.insert(0, str(self.read_start_pos) + 'S')
            self.cigar = ''.join(self.cigar_parts)

        if missing_bases_at_end:
            # Adjust the end of the reference.
            self.reference_end_pos += missing_bases_at_end
            self.reference_end_gap -= missing_bases_at_end

            # Adjust the end of the read and fix up the CIGAR to match.
            self.read_end_pos += missing_bases_at_end
            self.read_end_gap -= missing_bases_at_end
            self.cigar_parts.pop()
            if self.cigar_parts[-1][-1] == 'M':
                new_match_length = missing_bases_at_end + int(self.cigar_parts[-1][:-1])
                self.cigar_parts.pop()
                new_cigar_part = str(new_match_length) + 'M'
            else:
                new_cigar_part = str(missing_bases_at_end) + 'M'
            self.cigar_parts.append(new_cigar_part)
            if self.read_end_gap > 0:
                self.cigar_parts.append(str(self.read_end_gap) + 'S')
            self.cigar = ''.join(self.cigar_parts)

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

    def tally_up_alignment(self, references):
        '''
        Counts the matches, mismatches, indels and deletions. Also calculates the percent identity,
        which it does like BLAST: matches / alignment positions.
        '''
        # Get the aligned parts of the read and reference sequences.
        read_seq = self.full_read_sequence[self.read_start_pos:self.read_end_pos]
        ref_seq = references[self.reference_name][self.reference_start_pos:self.reference_end_pos]

        # Remove the soft clipping parts of the CIGAR string.
        cigar_parts = self.cigar_parts[:]
        if cigar_parts[0][-1] == 'S':
            cigar_parts.pop(0)
        if cigar_parts[-1][-1] == 'S':
            cigar_parts.pop()

        # Step through the alignment, counting as we go.
        read_i = 0
        ref_i = 0
        align_i = 0
        for cigar_part in cigar_parts:
            cigar_count = int(cigar_part[:-1])
            cigar_type = cigar_part[-1]
            if cigar_type == 'I':
                self.insertions += cigar_count
                read_i += cigar_count
            elif cigar_type == 'D':
                self.deletions += cigar_count
                ref_i += cigar_count
            else: # match/mismatch
                for _ in range(cigar_count):
                    if read_seq[read_i] == ref_seq[ref_i]:
                        self.matches += 1
                    else:
                        self.mismatches += 1
                    read_i += 1
                    ref_i += 1
            align_i += cigar_count
        self.percent_identity = 100.0 * self.matches / align_i

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

def load_alignments(sam_filename, references):
    '''
    This function returns a list of Alignment objects from the given SAM file.
    '''
    alignments = []
    sam_file = open(sam_filename, 'r')
    for line in sam_file:
        if not line.startswith('@') and line.split('\t', 3)[2] != '*':
            alignments.append(Alignment(line, references))
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

