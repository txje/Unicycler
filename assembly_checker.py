#!/usr/bin/env python
'''
Long read assembly checker

Author: Ryan Wick
email: rrwick@gmail.com
'''
from __future__ import print_function
from __future__ import division

import sys
import os
import re
import random
import argparse
import math
from multiprocessing import cpu_count

from semi_global_long_read_aligner import AlignmentScoringScheme, Read, Reference, load_references, \
                                          load_long_reads, quit_with_error, get_nice_header, \
                                          get_random_sequence_alignment_error_rates, \
                                          reverse_complement, int_to_str, float_to_str, \
                                          get_ref_shift_from_cigar_part, print_progress_line

'''
VERBOSITY controls how much the script prints to the screen.
'''
VERBOSITY = 0

def main():
    '''
    Script execution starts here.
    '''
    args = get_arguments()
    references = load_references(args.ref, VERBOSITY)
    reference_dict = {x.name: x for x in references}
    read_dict, _ = load_long_reads(args.reads, VERBOSITY)
    scoring_scheme = get_scoring_scheme_from_sam(args.sam)
    random_seq_error_rate = get_random_sequence_error_rate(scoring_scheme)
    alignments = load_sam_alignments(args.sam, read_dict, reference_dict, scoring_scheme)

    count_depth_and_errors_per_base(references, reference_dict, alignments)
    count_depth_and_errors_per_window(references, args.window)

    if VERBOSITY > 0:
        produce_console_output(references)
    if args.window_tables:
        produce_window_tables(references, args.window_tables)
    if args.base_tables:
        produce_base_tables(references, args.base_tables)
    if args.html:
        produce_html_files(references, args.html)

    sys.exit(0)

def get_arguments():
    '''
    Specifies the command line arguments required by the script.
    '''
    parser = argparse.ArgumentParser(description='Long read assembly checker',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('--sam', type=str, required=True, default=argparse.SUPPRESS,
                        help='Input SAM file of alignments')
    parser.add_argument('--ref', type=str, required=True, default=argparse.SUPPRESS,
                        help='FASTA file containing one or more reference sequences')
    parser.add_argument('--reads', type=str, required=True, default=argparse.SUPPRESS,
                        help='FASTQ file of long reads')
    parser.add_argument('--window', type=int, required=False, default=1000,
                        help='Window size for error summaries')
    parser.add_argument('--window_tables', type=str, required=False, default=argparse.SUPPRESS,
                        help='Path and/or prefix for table files summarising reference errors for '
                             'reference windows (default: do not save window tables)')
    parser.add_argument('--base_tables', type=str, required=False, default=argparse.SUPPRESS,
                        help='Path and/or prefix for table files summarising reference errors at '
                             'each base (default: do not save base tables)')
    parser.add_argument('--html', type=str, required=False, default=argparse.SUPPRESS,
                        help='Path and/or prefix for html files with plots (default: do not save '
                             'html files)')
    parser.add_argument('--threads', type=int, required=False, default=argparse.SUPPRESS,
                        help='Number of CPU threads used to align (default: the number of '
                             'available CPUs)')
    parser.add_argument('--verbosity', type=int, required=False, default=1,
                        help='Level of stdout information (0 to 2)')

    args = parser.parse_args()

    global VERBOSITY
    VERBOSITY = args.verbosity

    # If some arguments weren't set, set them to None/False. We don't use None/False as a default
    # in add_argument because it makes the help text look weird.
    try:
        args.window_tables
    except AttributeError:
        args.window_tables = None
    try:
        args.base_tables
    except AttributeError:
        args.base_tables = None
    try:
        args.html
    except AttributeError:
        args.html = None
    try:
        args.threads
    except AttributeError:
        args.threads = cpu_count()
        if VERBOSITY > 2:
            print('\nThread count set to', args.threads)

    return args

def get_scoring_scheme_from_sam(sam_filename):
    '''
    Looks for the 'SC' tag in the SAM file to get the alignment scoring scheme.
    '''
    sam_file = open(sam_filename, 'r')
    for line in sam_file:
        line = line.strip()

        # If we've reached the end of the header and still not found the scoring scheme, just
        # return a simple generic one.
        if not line.startswith('@'):
            return AlignmentScoringScheme('1,-1,-1,-1')

        line_parts = line.split('\t')
        for part in line_parts:
            if part.startswith('SC:'):
                scoring_scheme_string = part[3:]
                if scoring_scheme_string.count(',') == 3:
                    return AlignmentScoringScheme(scoring_scheme_string)

    return AlignmentScoringScheme('1,-1,-1,-1')

def get_random_sequence_error_rate(scoring_scheme):
    '''
    Returns the expected number of errors per reference base for an alignment of random sequences
    using the given scoring scheme.
    '''
    # I've precalculated the error rate for some typical scoring schemes.
    scoring_scheme_str = str(scoring_scheme)
    if scoring_scheme_str == '1,0,0,0':
        return 0.587
    elif scoring_scheme_str == '0,-1,-1,-1':
        return 0.526
    elif scoring_scheme_str == '1,-1,-1,-1':
        return 0.533
    elif scoring_scheme_str == '5,-4,-8,-6':
        return 0.527
    elif scoring_scheme_str == '5,-6,-10,0':
        return 1.012
    elif scoring_scheme_str == '2,-5,-2,-1':
        return 0.713
    elif scoring_scheme_str == '1,-3,-5,-2':
        return 0.544
    elif scoring_scheme_str == '5,-11,-2,-4':
        return 0.707
    elif scoring_scheme_str == '3,-6,-5,-2':
        return 0.641
    elif scoring_scheme_str == '2,-3,-5,-2':
        return 0.546
    elif scoring_scheme_str == '1,-2,0,0':
        return 0.707
    elif scoring_scheme_str == '0,-6,-5,-3':
        return 0.575
    elif scoring_scheme_str == '2,-6,-5,-3':
        return 0.578
    elif scoring_scheme_str == '1,-4,-6,-1':
        return 0.812

    # If the scoring scheme doesn't match a previously known one, we will use the C++ code to get
    # an error rate estimate.
    else:
        error_rate_str = get_random_sequence_alignment_error_rates(1000, 100, scoring_scheme)
        return float(error_rate_str.split('\n')[1].split('\t')[6])

def load_sam_alignments(sam_filename, read_dict, reference_dict, scoring_scheme):
    '''
    This function returns a list of Alignment objects from the given SAM file.
    '''
    if VERBOSITY > 0:
        print('Loading alignments')
        print('------------------')
        num_alignments = sum(1 for line in open(sam_filename) if not line.startswith('@'))
        print_progress_line(0, num_alignments)

    sam_alignments = []
    sam_file = open(sam_filename, 'r')
    for line in sam_file:
        line = line.strip()
        if line and not line.startswith('@') and line.split('\t', 3)[2] != '*':
            sam_alignments.append(Alignment(line, read_dict, reference_dict, scoring_scheme))
            if VERBOSITY > 0:
                print_progress_line(len(sam_alignments), num_alignments)

    if VERBOSITY > 0:
        if len(sam_alignments) < num_alignments:
            print_progress_line(len(sam_alignments), len(sam_alignments))
        print('\n')
    return sam_alignments

def count_depth_and_errors_per_base(references, reference_dict, alignments):
    '''
    Counts up the depth and errors for each base of each reference and stores the counts in the
    Reference objects.
    '''
    for ref in references:
        ref_length = ref.get_length()
        ref.depths = [0] * ref_length
        ref.error_counts = [0] * ref_length
        ref.error_rates = [0.0] * ref_length
        ref.alignment_count = 0

    for alignment in alignments:
        ref = reference_dict[alignment.ref.name]
        ref.alignment_count += 1
        for i in range(alignment.ref_start_pos, alignment.ref_end_pos):
            ref.depths[i] += 1
        for i in alignment.ref_mismatch_positions:
            ref.error_counts[i] += 1
        for i in alignment.ref_insertion_positions:
            ref.error_counts[i] += 1
        for i in alignment.ref_deletion_positions:
            ref.error_counts[i] += 1

    for ref in references:
        ref_length = ref.get_length()
        for i in range(ref_length):
            if ref.depths[i] > 0:
                ref.error_rates[i] = ref.error_counts[i] / ref.depths[i]

def count_depth_and_errors_per_window(references, window_size):
    '''
    Counts up the depth and errors for each window of each reference and stores the counts in the
    Reference objects.
    '''
    for ref in references:
        ref_length = ref.get_length()
        window_count = max(1, int(round(ref_length / window_size)))
        ref.window_size = ref_length / window_count

        ref.window_starts = []
        ref.window_depths = []
        ref.window_error_rates = []
        ref.min_window_depth = None
        ref.min_window_error_rate = None
        ref.max_window_depth = 0.0
        ref.max_window_error_rate = 0.0

        for i in range(window_count):
            window_start = int(round(ref.window_size * i))
            window_end = int(round(ref.window_size * (i + 1)))
            ref.window_starts.append(window_start)
            this_window_size = window_end - window_start

            total_window_depth = 0
            total_window_error_rate = 0
            for j in range(window_start, window_end):
                total_window_depth += ref.depths[j]
                total_window_error_rate += ref.error_rates[j]

            window_depth = total_window_depth / this_window_size
            window_error_rate = total_window_error_rate / this_window_size

            ref.window_depths.append(window_depth)
            ref.window_error_rates.append(window_error_rate)

            if ref.min_window_depth is None:
                ref.min_window_depth = window_depth
            else:
                ref.min_window_depth = min(window_depth, ref.min_window_depth)

            if ref.min_window_error_rate is None:
                ref.min_window_error_rate = window_error_rate
            else:
                ref.min_window_error_rate = min(window_error_rate, ref.min_window_error_rate)

            ref.max_window_depth = max(window_depth, ref.max_window_depth)
            ref.max_window_error_rate = max(window_error_rate, ref.max_window_error_rate)

def produce_console_output(references):
    '''
    Write a summary of the results to std out.
    '''
    for ref in references:
        print(ref.name)
        print('-' * len(ref.name))
        ref_length = ref.get_length()
        max_v = max(100, ref_length)

        print('Length:         ', int_to_str(ref_length, max_v) + ' bp')
        print('Alignments:     ', int_to_str(ref.alignment_count, max_v))
        print('Min depth:      ', float_to_str(ref.min_window_depth, 2, max_v))
        print('Max depth:      ', float_to_str(ref.max_window_depth, 2, max_v))
        print('Min error rate: ', float_to_str(ref.min_window_error_rate * 100.0, 2, max_v) + '%')
        print('Max error rate: ', float_to_str(ref.max_window_error_rate * 100.0, 2, max_v) + '%')

        print()




def produce_window_tables(references, window_tables_filename):
    '''
    Write tables of depth and error rates per reference window.
    '''





def produce_base_tables(references, base_tables_filename):
    '''
    Write tables of depth and error counts per reference base.
    '''





def produce_html_files(references, html_filename):
    '''
    Write html files containing plots of results.
    '''






# def summarise_errors(references, long_reads, table_prefix):
#     '''
#     Writes a table file summarising the alignment errors in terms of reference sequence position.
#     Works in a brute force manner - could be made more efficient later if necessary.
#     '''
#     # If we are not making table files or printing a summary, then quit now because there's nothing
#     # else to do.
#     if not table_prefix and VERBOSITY == 0:
#         return

#     # # We are be willing to throw out the worst alignments for any particular location. This value
#     # # specifies how many of our alignments we'll keep.
#     # # E.g. if 0.75, we'll throw out up to 25% of the alignments for each reference location.
#     # frac_alignments_to_keep = percent_alignments_to_keep / 100.0

#     # Ensure the table prefix is nicely formatted and any necessary directories are made.
#     if table_prefix:
#         if os.path.isdir(table_prefix) and not table_prefix.endswith('/'):
#             table_prefix += '/'
#         if table_prefix.endswith('/') and not os.path.isdir(table_prefix):
#             os.makedirs(table_prefix)
#         if not table_prefix.endswith('/'):
#             directory = os.path.dirname(table_prefix)
#             if directory and not os.path.isdir(directory):
#                 os.makedirs(directory)

#     if VERBOSITY > 0:
#         # So the numbers align nicely, we look for and remember the largest number to be displayed.
#         max_v = max(100,
#                     sum([len(x.alignments) for x in long_reads.itervalues()]),
#                     max([len(x.sequence) for x in references]))
#         print('Alignment summaries per reference')
#         print('---------------------------------')

#     for reference in references:
#         # Create a table file for the reference.
#         if table_prefix:
#             header_for_filename = clean_str_for_filename(reference.name)
#             if table_prefix.endswith('/'):
#                 table_filename = table_prefix + header_for_filename + '.txt'
#             else:
#                 table_filename = table_prefix + '_' + header_for_filename + '.txt'
#             table = open(table_filename, 'w')
#             table.write('\t'.join(['base',
#                                    'read depth',
#                                    'mismatches',
#                                    'deletions',
#                                    'insertions',
#                                    'mismatch rate',
#                                    'deletion rate',
#                                    'insertion rate',
#                                    'insertion sizes']) + '\n')

#         # Gather up the alignments for this reference and count up the depth for each reference
#         # position.
#         seq_len = len(reference.sequence)
#         depths = [0] * seq_len
#         alignments = []
#         for read in long_reads.itervalues():
#             for alignment in read.alignments:
#                 if alignment.ref.name == reference.name:
#                     alignments.append(alignment)
#                     for pos in xrange(alignment.ref_start_pos, alignment.ref_end_pos):
#                         depths[pos] += 1

#         # # Discard as many of the worst alignments as we can while keeping a sufficient read depth
#         # # at each position.
#         # sufficient_depths = [math.ceil(frac_alignments_to_keep * x) for x in total_depths]
#         # filtered_depths = total_depths[:]
#         # filtered_alignments = []
#         # alignments = sorted(alignments, key=lambda x: x.scaled_score)
#         # for alignment in alignments:
#         #     removal_okay = True
#         #     for pos in xrange(alignment.ref_start_pos, alignment.ref_end_pos):
#         #         if filtered_depths[pos] - 1 < sufficient_depths[pos]:
#         #             removal_okay = False
#         #             break
#         #     if removal_okay:
#         #         for pos in xrange(alignment.ref_start_pos, alignment.ref_end_pos):
#         #             filtered_depths[pos] -= 1
#         #     else:
#         #         filtered_alignments.append(alignment)

#         # Determine how much of the reference sequence has at least one aligned read.
#         aligned_len = sum([1 for x in depths if x > 0])
#         aligned_percent = 100.0 * aligned_len / seq_len

#         # Using the filtered alignments, add up mismatches, deletions and insertions and also get
#         # their rates (divided by the depth).
#         # Insertions are counted in two ways:
#         #   1) With multi-base insertions totalled up. E.g. 3I results in 3 counts all in the same
#         #      reference location.
#         #   2) With multi-base insertions counted only once. E.g. 3I results in 1 count at the
#         #      reference location.
#         mismatches = [0] * seq_len
#         deletions = [0] * seq_len
#         insertions = [0] * seq_len
#         insertion_sizes = {}

#         for alignment in alignments:
#             for pos in alignment.ref_mismatch_positions:
#                 mismatches[pos] += 1
#             for pos in alignment.ref_deletion_positions:
#                 deletions[pos] += 1
#             for pos, size in alignment.ref_insertion_positions_and_sizes:
#                 insertions[pos] += 1
#                 if pos not in insertion_sizes:
#                     insertion_sizes[pos] = []
#                 insertion_sizes[pos].append(size)
#         mismatch_rates = []
#         deletion_rates = []
#         insertion_rates = []
#         insertion_rates_multi_base = []
#         for i in xrange(seq_len):
#             mismatch_rate = 0.0
#             deletion_rate = 0.0
#             insertion_rate = 0.0
#             insertion_rate_multi_base = 0.0
#             if depths[i] > 0:
#                 mismatch_rate = mismatches[i] / depths[i]
#                 deletion_rate = deletions[i] / depths[i]
#                 insertion_rate = insertions[i] / depths[i]
#                 if i in insertion_sizes:
#                     insertion_rate_multi_base = sum(insertion_sizes[i]) / depths[i]
#             mismatch_rates.append(mismatch_rate)
#             deletion_rates.append(deletion_rate)
#             insertion_rates.append(insertion_rate)
#             insertion_rates_multi_base.append(insertion_rate_multi_base)

#         if table_prefix:
#             for i in xrange(seq_len):
#                 insertion_sizes_str = ''
#                 if i in insertion_sizes:
#                     insertion_sizes_str = ', '.join([str(x) for x in insertion_sizes[i]])
#                 table.write('\t'.join([str(i+1),
#                                        str(depths[i]),
#                                        str(mismatches[i]),
#                                        str(deletions[i]),
#                                        str(insertions[i]),
#                                        str(mismatch_rates[i]),
#                                        str(deletion_rates[i]),
#                                        str(insertion_rates[i]),
#                                        insertion_sizes_str]))
#                 table.write('\n')
#             table.close()

#         if VERBOSITY > 0:
#             contained_alignment_count = 0
#             overlapping_alignment_count = 0
#             for alignment in alignments:
#                 if alignment.is_whole_read():
#                     contained_alignment_count += 1
#                 else:
#                     overlapping_alignment_count += 1
#             print(reference.name)
#             print('  Length:       ', int_to_str(seq_len, max_v) + ' bp')
#             if alignments:
#                 mean_depth = sum(depths) / seq_len
#                 mean_mismatch_rate = 100.0 * sum(mismatch_rates) / aligned_len
#                 mean_insertion_rate = 100.0 * sum(insertion_rates_multi_base) / aligned_len
#                 mean_deletion_rate = 100.0 * sum(deletion_rates) / aligned_len
#                 if VERBOSITY > 0:
#                     print('  Alignments:         ', int_to_str(len(alignments), max_v))
#                 if VERBOSITY > 1:
#                     print('    Contained:        ', int_to_str(contained_alignment_count, max_v))
#                     print('    Overlapping:      ', int_to_str(overlapping_alignment_count, max_v))
#                 print('  Covered length:     ', int_to_str(aligned_len, max_v) + ' bp')
#                 print('  Covered fraction:   ', float_to_str(aligned_percent, 2, max_v) + '%')
#                 if VERBOSITY > 1:
#                     print('  Mean read depth:    ', float_to_str(mean_depth, 2, max_v))
#                     print('  Mismatch rate:      ',
#                           float_to_str(mean_mismatch_rate, 2, max_v) +'%')
#                     print('  Insertion rate:     ',
#                           float_to_str(mean_insertion_rate, 2, max_v) + '%')
#                     print('  Deletion rate:      ',
#                           float_to_str(mean_deletion_rate, 2, max_v) + '%')
#             else:
#                 print('  Alignments:         ', int_to_str(0, max_v))
#             print()

class Alignment(object):
    '''
    This class describes an alignment between a long read and a reference.
    '''
    def __init__(self, sam_line, read_dict, reference_dict, scoring_scheme):

        # Grab the important parts of the alignment from the SAM line.
        sam_parts = sam_line.split('\t')
        self.rev_comp = bool(int(sam_parts[1]) & 0x10)
        self.cigar = sam_parts[5]
        self.cigar_parts = re.findall(r'\d+\w', self.cigar)

        self.read = read_dict[sam_parts[0]]
        self.read_start_pos = self.get_start_soft_clips()
        self.read_end_pos = self.read.get_length() - self.get_end_soft_clips()
        self.read_end_gap = self.get_end_soft_clips()

        self.ref = reference_dict[get_nice_header(sam_parts[2])]
        self.ref_start_pos = int(sam_parts[3]) - 1
        self.ref_end_pos = self.ref_start_pos
        for cigar_part in self.cigar_parts:
            self.ref_end_pos += get_ref_shift_from_cigar_part(cigar_part)
        if self.ref_end_pos > len(self.ref.sequence):
            self.ref_end_pos = len(self.ref.sequence)
        self.ref_end_gap = len(self.ref.sequence) - self.ref_end_pos

        self.ref_mismatch_positions = []
        self.ref_deletion_positions = []
        self.ref_insertion_positions = []

        # Remove the soft clipping parts of the CIGAR string for tallying.
        cigar_parts = self.cigar_parts[:]
        if cigar_parts[0][-1] == 'S':
            cigar_parts.pop(0)
        if cigar_parts and cigar_parts[-1][-1] == 'S':
            cigar_parts.pop()
        if not cigar_parts:
            return

        if self.rev_comp:
            read_seq = reverse_complement(self.read.sequence)
        else:
            read_seq = self.read.sequence

        read_i = self.read_start_pos
        ref_i = self.ref_start_pos

        for cigar_part in cigar_parts:
            cigar_count = int(cigar_part[:-1])
            cigar_type = cigar_part[-1]
            if cigar_type == 'I':
                for _ in xrange(cigar_count):
                    self.ref_insertion_positions.append(ref_i)
                read_i += cigar_count
            elif cigar_type == 'D':
                for i in xrange(cigar_count):
                    self.ref_deletion_positions.append(ref_i + i)
                ref_i += cigar_count
            else: # match/mismatch
                for _ in xrange(cigar_count):
                    # If all is good with the CIGAR, then we should never end up with a sequence
                    # index out of the sequence range. But a CIGAR error (which has occurred in
                    # GraphMap) can cause this, so check here.
                    if read_i >= self.read.get_length() or ref_i >= self.ref.get_length():
                        break
                    if read_seq[read_i] != self.ref.sequence[ref_i]:
                        self.ref_mismatch_positions.append(ref_i)
                    read_i += 1
                    ref_i += 1

    def __repr__(self):
        read_start, read_end = self.read_start_end_positive_strand()
        return_str = self.read.name + ' (' + str(read_start) + '-' + str(read_end) + ', '
        if self.rev_comp:
            return_str += 'strand: -), '
        else:
            return_str += 'strand: +), '
        return_str += self.ref.name + ' (' + str(self.ref_start_pos) + '-' + \
                      str(self.ref_end_pos) + ')'
        error_count = len(self.ref_mismatch_positions) + len(self.ref_deletion_positions) + \
                      len(self.ref_insertion_positions)
        return_str += ', errors = ' + int_to_str(error_count)
        return return_str

    def get_start_soft_clips(self):
        '''
        Returns the number of soft-clipped bases at the start of the alignment.
        '''
        if self.cigar_parts[0][-1] == 'S':
            return int(self.cigar_parts[0][:-1])
        else:
            return 0

    def get_end_soft_clips(self):
        '''
        Returns the number of soft-clipped bases at the start of the alignment.
        '''
        if self.cigar_parts[-1][-1] == 'S':
            return int(self.cigar_parts[-1][:-1])
        else:
            return 0

    def read_start_end_positive_strand(self):
        '''
        This function returns the read start/end coordinates for the positive strand of the read.
        For alignments on the positive strand, this is just the normal start/end. But for
        alignments on the negative strand, the coordinates are flipped to the other side.
        '''
        if not self.rev_comp:
            return self.read_start_pos, self.read_end_pos
        else:
            start = self.read.get_length() - self.read_end_pos
            end = self.read.get_length() - self.read_start_pos
            return start, end

if __name__ == '__main__':
    main()
