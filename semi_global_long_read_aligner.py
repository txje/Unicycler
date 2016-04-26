#!/usr/bin/env python
'''
Semi-global long read aligner

This is a script to align error-prone long reads (e.g. PacBio or Nanopore) to one or more
references in a semi-global manner. It uses GraphMap to do the actual alignment and then filters
the resulting alignments and outputs summaries.

Semi-global alignment allows for unpenalised end gaps, but the alignment will continue until one of
the two sequences ends. This includes cases where the two sequences overlap and cases where one
sequence is contained within the other:

  AAAAA        AAAAAAAAAAA         AAAAAAAA     AAAAAAAA
  |||||          |||||||           |||||           |||||
BBBBBBBBB        BBBBBBB       BBBBBBBBB           BBBBBBBBB

This tool is intended for cases where the reads and reference are expected to match perfectly (or
at least as perfectly as error-prone long reads can match). An example of an appropriate case would
be if the reference sequences are assembled contigs of a bacterial strain and the long reads are
from the same strain.

Required inputs:
  1) FASTA file of one or more reference sequences
  2) FASTQ file of long reads

Required outputs:
  1) SAM file of alignments

Optional outputs:
  1) SAM file of raw (unfiltered) GraphMap alignments
  2) Table files of depth and errors per base of each reference

Author: Ryan Wick
email: rrwick@gmail.com
'''
from __future__ import print_function
from __future__ import division

import subprocess
import sys
import os
import re
import random
import argparse
import string
import math
from ctypes import CDLL, cast, c_char_p, c_int, c_double, c_void_p
from multiprocessing.dummy import Pool as ThreadPool

'''
VERBOSITY controls how much the script prints to the screen.
0 = nothing is printed
1 = a relatively simple output is printed
2 = a more thorough output is printed, including details on each Seqan alignment
3 = even more output is printed, including stuff from the C++ code
4 = tons of stuff is printed, including all k-mer positions in each Seqan alignment
'''
VERBOSITY = 0

C_LIB = CDLL(os.path.join(os.path.dirname(os.path.realpath(__file__)), 'seqan_align.so'))

'''
This function conducts a semi-global alignment surrounding the given line. Since the search is
limited to a narrow band it is much more efficient than an exhaustive search.
'''
C_LIB.bandedSemiGlobalAlignment.argtypes = [c_char_p, # Sequence 1
                                            c_char_p, # Sequence 2
                                            c_int,    # Sequence 1 length
                                            c_int,    # Sequence 2 length
                                            c_double, # Slope
                                            c_double, # Intercept
                                            c_int,    # K-mer size
                                            c_int,    # Band size
                                            c_int,    # Verbosity
                                            c_int,    # Match score
                                            c_int,    # Mismatch score
                                            c_int,    # Gap open score
                                            c_int,    # Gap extension score
                                            c_char_p] # K-mer location string
C_LIB.bandedSemiGlobalAlignment.restype = c_void_p    # String describing alignment

'''
This function looks for the most likely line representing the alignment. It is used to get a
line to be given to C_LIB.bandedSemiGlobalAlignment.
'''
C_LIB.findAlignmentLines.argtypes = [c_char_p, # Sequence 1
                                     c_char_p, # Sequence 2
                                     c_int,    # Sequence 1 length
                                     c_int,    # Sequence 2 length
                                     c_double, # Expected slope
                                     c_int]    # Verbosity
C_LIB.findAlignmentLines.restype = c_void_p    # String describing the found line(s)

'''
This function is used to conduct a short alignment for the sake of extending a GraphMap alignment.
'''
C_LIB.startExtensionAlignment.argtypes = [c_char_p, # Sequence 1
                                          c_char_p, # Sequence 2
                                          c_int,    # Sequence 1 length
                                          c_int,    # Sequence 2 length
                                          c_int,    # Verbosity
                                          c_int,    # Match score
                                          c_int,    # Mismatch score
                                          c_int,    # Gap open score
                                          c_int]    # Gap extension score
C_LIB.startExtensionAlignment.restype = c_void_p    # String describing alignment

'''
This function cleans up the heap memory for the C strings returned by the other C functions. It
must be called after them.
'''
C_LIB.free_c_string.argtypes = [c_void_p]
C_LIB.free_c_string.restype = None





def main():
    '''
    If this script is run on its own, execution starts here.
    '''
    args = get_arguments()
    check_file_exists(args.ref)
    check_file_exists(args.reads)
    temp_dir_exist_at_start = os.path.exists(args.temp_dir)
    if not temp_dir_exist_at_start:
        os.makedirs(args.temp_dir)

    global VERBOSITY
    VERBOSITY = args.verbosity
    if VERBOSITY > 0:
        print()

    references = load_references(args.ref)
    reads = load_long_reads(args.reads)

    reads = semi_global_align_long_reads(references, args.ref, reads, args.reads, args.sam,
                                         args.temp_dir, args.path, args.threads, args.partial_ref,
                                         AlignmentScoringScheme(args.scores))
    summarise_errors(references, reads, args.table, args.keep)

    if not temp_dir_exist_at_start:
        os.rmdir(args.temp_dir)
    sys.exit(0)

def get_arguments():
    '''
    Specifies the command line arguments required by the script.
    '''
    parser = argparse.ArgumentParser(description='Semi-global long read aligner',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--ref', type=str, required=True, default=argparse.SUPPRESS,
                        help='FASTA file containing one or more reference sequences')
    parser.add_argument('--partial_ref', action='store_true',
                        help='Set if some reads are not expected to align to the reference.')
    parser.add_argument('--reads', type=str, required=True, default=argparse.SUPPRESS,
                        help='FASTQ file of long reads')
    parser.add_argument('--sam', type=str, required=True, default=argparse.SUPPRESS,
                        help='SAM file of resulting alignments')
    parser.add_argument('--table', type=str, required=False,
                        help='Path and/or prefix for table files summarising reference errors')
    parser.add_argument('--temp_dir', type=str, required=False, default='align_temp',
                        help='Temp directory for working files')
    parser.add_argument('--path', type=str, required=False, default='graphmap',
                        help='Path to the GraphMap executable')
    parser.add_argument('--scores', type=str, required=False, default='3,-6,-5,-2',
                        help='Alignment scores: match, mismatch, gap open, gap extend')
    parser.add_argument('--threads', type=int, required=False, default=8,
                        help='Number of threads used by GraphMap')
    parser.add_argument('--verbosity', type=int, required=False, default=1,
                        help='Level of stdout information (0 to 4)')
    parser.add_argument('--keep', type=float, required=False, default=75.0,
                        help='Percentage of alignments to keep')
    return parser.parse_args()

def semi_global_align_long_reads(references, ref_fasta, reads, reads_fastq, output_sam, temp_dir,
                                 graphmap_path, threads, partial_ref, scoring_scheme):
    '''
    This function does the primary work of this module: aligning long reads to references in an
    end-gap-free, semi-global manner. It returns a list of LongRead objects which contain their
    alignments.
    If seqan_all is True, then every Alignment object will be refined by using Seqan.
    If seqan_all is False, then only the overlap alignments and a small set of long contained
    alignments will be run through Seqan.
    '''
    graphmap_sam = os.path.join(temp_dir, 'graphmap_alignments.sam')
    # run_graphmap(ref_fasta, reads_fastq, graphmap_sam, graphmap_path, threads, scoring_scheme)
    graphmap_alignments = load_sam_alignments(graphmap_sam, reads, references, scoring_scheme)

    if VERBOSITY > 2 and graphmap_alignments:
        print('All GraphMap alignments')
        print('-----------------------')
        for alignment in graphmap_alignments:
            print(alignment)
            if VERBOSITY > 3:
                print(alignment.cigar)
        print()

    if VERBOSITY > 0:
        print_alignment_summary_table(graphmap_alignments)

    semi_global_graphmap_alignments = extend_to_semi_global(graphmap_alignments,
                                                                       scoring_scheme)




    # Give the alignments to their corresponding reads.
    for alignment in semi_global_graphmap_alignments:
        reads[alignment.read.name].alignments.append(alignment)












    write_sam_file(semi_global_graphmap_alignments, graphmap_sam, output_sam) # TEMP














#     # We now choose some alignments to use for the identity mean and std dev calculation.
#     # We want long alignments that haven't been extended much (i.e. for which GraphMap already
#     # took the alignment close to the full length).
#     sorted_alignments = sorted(paf_alignments, reverse=True,
#                                key=lambda x: (x.ref_end_pos - x.ref_start_pos) /
#                                (x.extension_length + 1))
#     target_alignment_count = 20 # TO DO: MAKE THIS A PARAMETER?
#     alignments_for_identity = sorted_alignments[:target_alignment_count]

#     if VERBOSITY > 2 and alignments_for_identity:
#         print('Alignments to be used for identity mean and std dev')
#         print('---------------------------------------------------')
#         for alignment in alignments_for_identity:
#             print(alignment)
#         print()

#     # Determine which alignments to refine in Seqan.  If seqan_all is set, we refine them all. If
#     # not, then we only refine overlapping alignments and the alignments that we'll use to get an
#     # identity mean and std dev.
#     if seqan_all:
#         alignments_to_refine = paf_alignments
#         alignments_not_to_refine = []
#         if VERBOSITY > 0 and alignments_to_refine:
#             print('Refining all alignments with Seqan')
#             print('----------------------------------')
#     else:
#         alignments_to_refine = [x for x in paf_alignments
#                                 if x in alignments_for_identity or not x.is_whole_read()]
#         alignments_not_to_refine = [x for x in paf_alignments if x not in alignments_to_refine]
#         if VERBOSITY > 0 and alignments_to_refine:
#             count_str = str(len(alignments_to_refine))
#             print('Refining ' + count_str + ' alignments with Seqan')
#             print('-------------------------------' + '-' * len(count_str))

#     # Refine the appropriate alignments using Seqan.
#     alignments = []
#     refined_alignments_for_identity = []
#     completed_count = 0
#     if VERBOSITY == 1 and alignments_to_refine:
#         print_progress_line(completed_count, len(alignments_to_refine))

#     # If single-threaded, just do the work in a simple loop.
#     if threads == 1:
#         for alignment in alignments_to_refine:
#             refined_alignments, output = seqan_from_paf_alignment(alignment, reads, references,
#                                                                   median_ref_to_read_ratio)
#             alignments += refined_alignments
#             if alignment in alignments_for_identity:
#                 refined_alignments_for_identity += refined_alignments
#             completed_count += 1
#             if VERBOSITY == 1:
#                 print_progress_line(completed_count, len(alignments_to_refine))
#             if VERBOSITY > 1:
#                 print(output, end='')

#     # If multi-threaded, use a thread pool.
#     else:
#         pool = ThreadPool(threads)
#         arg_list = []
#         for alignment in alignments_to_refine:
#             arg_list.append((alignment, reads, references, median_ref_to_read_ratio))
#         for refined_alignments, output in pool.imap_unordered(seqan_from_paf_alignment_one_arg,
#                                                               arg_list, 1):
#             alignments += refined_alignments
#             if alignment in alignments_for_identity:
#                 refined_alignments_for_identity += refined_alignments
#             completed_count += 1
#             if VERBOSITY == 1:
#                 print_progress_line(completed_count, len(alignments_to_refine))
#             if VERBOSITY > 1:
#                 print(output, end='')

#     alignments += alignments_not_to_refine
#     if VERBOSITY == 1:
#         print('\n')

#     # Give the alignments to their corresponding reads.
#     for alignment in alignments:
#         reads[alignment.read.name].alignments.append(alignment)

#     # Filter the alignments based on conflicting read position.
#     if VERBOSITY > 0:
#         print('Filtering alignments')
#         print('--------------------')
#         print('Alignments before filtering:        ', int_to_str(len(alignments), max_v))
#     filtered_alignments = []
#     for read in reads.itervalues():
#         read.remove_conflicting_alignments()
#         filtered_alignments += read.alignments
#     if VERBOSITY > 0:
#         print('Alignments after conflict filtering:', int_to_str(len(filtered_alignments), max_v))

#     # If any of our alignments for identity mean/stddev were removed in conflict filtering then we
#     # don't want to use them.
#     refined_alignments_for_identity = [x for x in refined_alignments_for_identity \
#                                        if x in filtered_alignments]
#     mean_id, std_dev_id = get_mean_and_st_dev_identity(refined_alignments_for_identity)

#     # Filter the alignments based on identity.
#     if mean_id == 0.0 or std_dev_id == -1:
#         low_id_cutoff = 75.0
#         if VERBOSITY > 0:
#             print('Not enough alignments to automatically set a low identity cutoff. Using 75%.')
#     else:
#         std_devs_below_mean = 3.0 # TO DO: make this a parameter?
#         low_id_cutoff = mean_id - (std_devs_below_mean * std_dev_id)
#         if VERBOSITY > 0:
#             print('Complete alignment identity mean:   ', float_to_str(mean_id, max_v) + '%')
#             print('              standard deviation:   ', float_to_str(std_dev_id, max_v) + '%')
#             print('Low identity cutoff:                ', float_to_str(low_id_cutoff, max_v) + '%')
#     filtered_alignments = []
#     for read in reads.itervalues():
#         read.remove_low_id_alignments(low_id_cutoff)
#         filtered_alignments += read.alignments
#     if VERBOSITY > 0:
#         print('Alignments after identity filtering:', int_to_str(len(filtered_alignments), max_v))
#         print()
   
#     fully_aligned, partially_aligned, unaligned = group_reads_by_fraction_aligned(reads)

#     if VERBOSITY > 0:
#         print('Read summary')
#         print('------------')
#         max_v = len(reads)
#         print('Total read count:       ', int_to_str(len(reads), max_v))
#         print('Fully aligned reads:    ', int_to_str(len(fully_aligned), max_v))
#         print('Partially aligned reads:', int_to_str(len(partially_aligned), max_v))
#         print('Unaligned reads:        ', int_to_str(len(unaligned), max_v))
#         print()

#     # # Try to realign any reads which are not completely aligned.
#     # incomplete_reads = partially_aligned + unaligned
#     # if incomplete_reads:
#     #     completed_count = 0
#     #     if VERBOSITY > 0:
#     #         print('Attempting realignment of incomplete reads')
#     #         print('------------------------------------------')
#     #     if VERBOSITY == 1:
#     #         print_progress_line(completed_count, len(incomplete_reads))
#     #     new_alignments = []
#     #     for read in incomplete_reads:
#     #         new_alignments += seqan_alignment_one_read_all_refs(reads, references, read,
#     #                                                             median_ref_to_read_ratio,
#     #                                                             exhaustive)
#     #         completed_count += 1
#     #         if VERBOSITY == 1:
#     #             print_progress_line(completed_count, len(incomplete_reads))
#     #     if VERBOSITY == 1:
#     #         print('\n')

#     #     # Give the alignments to their corresponding reads.
#     #     for alignment in new_alignments:
#     #         reads[alignment.read.name].alignments.append(alignment)
#     #     all_alignments_count = sum([len(x.alignments) for x in reads.itervalues()])

#     #     # Filter the alignments based on conflicting read position.
#     #     if VERBOSITY > 0:
#     #         print('Filtering alignments')
#     #         print('--------------------')
#     #         print('Alignments before filtering:        ', int_to_str(all_alignments_count, max_v))
#     #     filtered_alignments = []
#     #     for read in reads.itervalues():
#     #         read.remove_conflicting_alignments()
#     #         filtered_alignments += read.alignments
#     #     if VERBOSITY > 0:
#     #         print('Alignments after conflict filtering:', int_to_str(len(filtered_alignments),
#     #                                                                  max_v))

#     #     fully_aligned, partially_aligned, unaligned = group_reads_by_fraction_aligned(reads)
#     #     if VERBOSITY == 1:
#     #         print()
#     #     if VERBOSITY > 0:
#     #         print('Updated read summary')
#     #         print('--------------------')
#     #         max_v = len(reads)
#     #         print('Total read count:       ', int_to_str(len(reads), max_v))
#     #         print('Fully aligned reads:    ', int_to_str(len(fully_aligned), max_v))
#     #         print('Partially aligned reads:', int_to_str(len(partially_aligned), max_v))
#     #         print('Unaligned reads:        ', int_to_str(len(unaligned), max_v))
#     #         print()

    # write_sam_file(filtered_alignments, output_sam)
    
    # os.remove(graphmap_sam)

    return reads

def summarise_errors(references, long_reads, table_prefix, percent_alignments_to_keep):
    '''
    Writes a table file summarising the alignment errors in terms of reference sequence position.
    Works in a brute force manner - could be made more efficient later if necessary.
    '''
    # If we are not making table files or printing a summary, then quit now because there's nothing
    # else to do.
    if not table_prefix and VERBOSITY == 0:
        return

    # We are be willing to throw out the worst alignments for any particular location. This value
    # specifies how many of our alignments we'll keep.
    # E.g. if 0.75, we'll throw out up to 25% of the alignments for each reference location.
    frac_alignments_to_keep = percent_alignments_to_keep / 100.0

    # Ensure the table prefix is nicely formatted and any necessary directories are made.
    if table_prefix:
        if os.path.isdir(table_prefix) and not table_prefix.endswith('/'):
            table_prefix += '/'
        if table_prefix.endswith('/') and not os.path.isdir(table_prefix):
            os.makedirs(table_prefix)
        if not table_prefix.endswith('/'):
            directory = os.path.dirname(table_prefix)
            if directory and not os.path.isdir(directory):
                os.makedirs(directory)

    if VERBOSITY > 0:
        # So the numbers align nicely, we look for and remember the largest number to be displayed.
        max_v = max(100,
                    sum([len(x.alignments) for x in long_reads.itervalues()]),
                    max([len(seq) for seq in references.itervalues()]))
        print('Alignment summaries per reference')
        print('---------------------------------')

    for header, seq in references.iteritems():
        # Create a table file for the reference.
        if table_prefix:
            header_for_filename = clean_str_for_filename(header)
            if table_prefix.endswith('/'):
                table_filename = table_prefix + header_for_filename + '.txt'
            else:
                table_filename = table_prefix + '_' + header_for_filename + '.txt'
            table = open(table_filename, 'w')
            table.write('\t'.join(['base',
                                   'total read depth',
                                   'filtered depth',
                                   'mismatches',
                                   'deletions',
                                   'insertions',
                                   'mismatch rate',
                                   'deletion rate',
                                   'insertion rate',
                                   'insertion sizes']) + '\n')

        # Gather up the alignments for this reference and count up the depth for each reference
        # position.
        seq_len = len(seq)
        total_depths = [0] * seq_len
        alignments = []
        for read in long_reads.itervalues():
            for alignment in read.alignments:
                if alignment.ref_name == header:
                    alignments.append(alignment)
                    for pos in xrange(alignment.ref_start_pos, alignment.ref_end_pos):
                        total_depths[pos] += 1

        # Discard as many of the worst alignments as we can while keeping a sufficient read depth
        # at each position.
        sufficient_depths = [math.ceil(frac_alignments_to_keep * x) for x in total_depths]
        filtered_depths = total_depths[:]
        filtered_alignments = []
        alignments = sorted(alignments, key=lambda x: x.scaled_score)
        for alignment in alignments:
            removal_okay = True
            for pos in xrange(alignment.ref_start_pos, alignment.ref_end_pos):
                if filtered_depths[pos] - 1 < sufficient_depths[pos]:
                    removal_okay = False
                    break
            if removal_okay:
                for pos in xrange(alignment.ref_start_pos, alignment.ref_end_pos):
                    filtered_depths[pos] -= 1
            else:
                filtered_alignments.append(alignment)

        # Determine how much of the reference sequence has at least one aligned read.
        aligned_len = sum([1 for x in filtered_depths if x > 0])
        aligned_percent = 100.0 * aligned_len / seq_len

        # Using the filtered alignments, add up mismatches, deletions and insertions and also get
        # their rates (divided by the depth).
        # Insertions are counted in two ways:
        #   1) With multi-base insertions totalled up. E.g. 3I results in 3 counts all in the same
        #      reference location.
        #   2) With multi-base insertions counted only once. E.g. 3I results in 1 count at the
        #      reference location.
        mismatches = [0] * seq_len
        deletions = [0] * seq_len
        insertions = [0] * seq_len
        insertion_sizes = {}

        for alignment in filtered_alignments:
            for pos in alignment.ref_mismatch_positions:
                mismatches[pos] += 1
            for pos in alignment.ref_deletion_positions:
                deletions[pos] += 1
            for pos, size in alignment.ref_insertion_positions_and_sizes:
                insertions[pos] += 1
                if pos not in insertion_sizes:
                    insertion_sizes[pos] = []
                insertion_sizes[pos].append(size)
        mismatch_rates = []
        deletion_rates = []
        insertion_rates = []
        insertion_rates_multi_base = []
        for i in xrange(seq_len):
            mismatch_rate = 0.0
            deletion_rate = 0.0
            insertion_rate = 0.0
            insertion_rate_multi_base = 0.0
            if filtered_depths[i] > 0:
                mismatch_rate = mismatches[i] / filtered_depths[i]
                deletion_rate = deletions[i] / filtered_depths[i]
                insertion_rate = insertions[i] / filtered_depths[i]
                if i in insertion_sizes:
                    insertion_rate_multi_base = sum(insertion_sizes[i]) / filtered_depths[i]
            mismatch_rates.append(mismatch_rate)
            deletion_rates.append(deletion_rate)
            insertion_rates.append(insertion_rate)
            insertion_rates_multi_base.append(insertion_rate_multi_base)

        if table_prefix:
            for i in xrange(seq_len):
                insertion_sizes_str = ''
                if i in insertion_sizes:
                    insertion_sizes_str = ', '.join([str(x) for x in insertion_sizes[i]])
                table.write('\t'.join([str(i+1),
                                       str(total_depths[i]),
                                       str(filtered_depths[i]),
                                       str(mismatches[i]),
                                       str(deletions[i]),
                                       str(insertions[i]),
                                       str(mismatch_rates[i]),
                                       str(deletion_rates[i]),
                                       str(insertion_rates[i]),
                                       insertion_sizes_str]))
                table.write('\n')
            table.close()

        if VERBOSITY > 0:
            contained_alignment_count = 0
            overlapping_alignment_count = 0
            for alignment in filtered_alignments:
                if alignment.is_whole_read():
                    contained_alignment_count += 1
                else:
                    overlapping_alignment_count += 1
            print(header)
            print('  Total length:       ', int_to_str(seq_len, max_v) + ' bp')
            if alignments:
                mean_depth = sum(filtered_depths) / seq_len
                mean_mismatch_rate = 100.0 * sum(mismatch_rates) / aligned_len
                mean_insertion_rate = 100.0 * sum(insertion_rates_multi_base) / aligned_len
                mean_deletion_rate = 100.0 * sum(deletion_rates) / aligned_len
                print('  Total alignments:   ', int_to_str(len(alignments), max_v))
                print('  Filtered alignments:', int_to_str(len(filtered_alignments), max_v))
                print('    Contained:        ', int_to_str(contained_alignment_count, max_v))
                print('    Overlapping:      ', int_to_str(overlapping_alignment_count, max_v))
                print('  Covered length:     ', int_to_str(aligned_len, max_v) + ' bp')
                print('  Covered fraction:   ', float_to_str(aligned_percent, 2, max_v) + '%')
                print('  Mean read depth:    ', float_to_str(mean_depth, 2, max_v))
                print('  Mismatch rate:      ', float_to_str(mean_mismatch_rate, 2, max_v) + '%')
                print('  Insertion rate:     ', float_to_str(mean_insertion_rate, 2, max_v) + '%')
                print('  Deletion rate:      ', float_to_str(mean_deletion_rate, 2, max_v) + '%')
            else:
                print('  Total alignments:   ', int_to_str(0, max_v))
            print()


def print_alignment_summary_table(graphmap_alignments):
    '''
    Prints a small table showing some details about the GraphMap alignments.
    '''
    read_to_refs = [x.get_read_to_ref_ratio() for x in graphmap_alignments]
    percent_ids = [x.percent_identity for x in graphmap_alignments]
    scores = [x.scaled_score for x in graphmap_alignments]
    read_to_ref_median, read_to_ref_mad = get_median_and_mad(read_to_refs)
    percent_id_median, percent_id_mad = get_median_and_mad(percent_ids)
    score_median, score_mad = get_median_and_mad(scores)

    print('Alignment summary')
    print('-----------------')
    print('Total GraphMap alignments:', int_to_str(len(graphmap_alignments)))
    print()

    table_lines = ['',
                   'Read / reference length:',
                   'Percent identity:',
                   'Score:']

    pad_length = max([len(x) for x in table_lines]) + 2
    table_lines = [x.ljust(pad_length) for x in table_lines]

    table_lines[0] += 'Median'
    table_lines[1] += float_to_str(read_to_ref_median, 3)
    table_lines[2] += float_to_str(percent_id_median, 2) + '%'
    table_lines[3] += float_to_str(score_median, 2)

    pad_length = max([len(x) for x in table_lines]) + 2
    table_lines = [x.ljust(pad_length) for x in table_lines]

    table_lines[0] += 'MAD'
    table_lines[1] += float_to_str(read_to_ref_mad, 3)
    table_lines[2] += float_to_str(percent_id_mad, 2) + '%'
    table_lines[3] += float_to_str(score_mad, 2)

    for line in table_lines:
        print(line)
    print()


def extend_to_semi_global(alignments, scoring_scheme):
    '''
    This function returns truly semi-global alignments made from the input alignments.
    '''
    allowed_missing_bases = 100 # TO DO: MAKE THIS A PARAMETER
    semi_global_alignments = []
    for alignment in alignments:
        total_missing_bases = alignment.get_total_missing_bases()

        # If an input alignment is already semi-global, then it's included in the output.
        if total_missing_bases == 0:
            semi_global_alignments.append(alignment)

        # If an input alignment is almost semi-global (below a threshold), and not too close to the
        # end of the reference, then it is extended to make it semi-global.
        elif total_missing_bases <= allowed_missing_bases:
            missing_start = alignment.get_missing_bases_at_start()
            missing_end = alignment.get_missing_bases_at_end()
            if missing_start and alignment.ref_start_pos > 2 * missing_start:
                alignment.extend_start(scoring_scheme)
            if missing_end and alignment.ref_end_gap > 2 * missing_end:
                alignment.extend_end(scoring_scheme)
            semi_global_alignments.append(alignment)

        # If an input alignment is above the threshold (not close to being semi-global), it is
        # discarded.
        else:
            pass

    return semi_global_alignments

def load_references(fasta_filename):
    '''
    This function loads in reference sequences from a FASTA file and returns a dictionary where
    key = nice reference name and value = reference sequence.
    '''
    references = {}
    total_bases = 0
    if VERBOSITY > 0:
        print('Loading references')
        print('------------------')
        num_refs = sum(1 for line in open(fasta_filename) if line.startswith('>'))
        if not num_refs:
            quit_with_error('There are no references sequences in ' + fasta_filename)
        print_progress_line(0, num_refs)
    fasta_file = open(fasta_filename, 'r')
    name = ''
    sequence = ''
    for line in fasta_file:
        line = line.strip()
        if not line:
            continue
        if line.startswith('>'): # Header line = start of new contig
            if name:
                references[get_nice_header(name)] = sequence
                total_bases += len(sequence)
                if VERBOSITY > 0:
                    print_progress_line(len(references), num_refs, total_bases)
                name = ''
                sequence = ''
            name = line[1:]
        else:
            sequence += line
    fasta_file.close()
    if name:
        references[get_nice_header(name)] = sequence
        total_bases += len(sequence)
        if VERBOSITY > 0:
            print_progress_line(len(references), num_refs, total_bases)
            print('\n')
    return references

def load_long_reads(fastq_filename):
    '''
    This function loads in long reads from a FASTQ file and returns a dictionary where key = read
    name and value = LongRead object.
    '''
    reads = {}
    total_bases = 0
    if VERBOSITY > 0:
        print('Loading reads')
        print('-------------')
        num_reads = sum(1 for line in open(fastq_filename)) // 4
        print_progress_line(0, num_reads)
    fastq = open(fastq_filename, 'r')
    for line in fastq:
        name = line.strip()[1:]
        sequence = next(fastq).strip()
        _ = next(fastq)
        qualities = next(fastq).strip()
        reads[name] = LongRead(name, sequence, qualities)
        total_bases += len(sequence)
        if VERBOSITY > 0:
            print_progress_line(len(reads), num_reads, total_bases)
    fastq.close()
    if VERBOSITY > 0:
        print('\n')
    return reads

# def graphmap_alignment(ref_fasta, long_reads_fastq, graphmap_sam, graphmap_path, working_dir,
#                        threads):
#     '''
#     This function runs GraphMap separately for each individual sequence in the reference.
#     Resulting alignments are collected in a single PAF file.
#     '''
#     references = load_fasta(ref_fasta)
#     if VERBOSITY > 0:
#         print('Raw GraphMap alignments per reference')
#         print('-------------------------------------')

#     final_paf = open(paf_file, 'w')
#     if VERBOSITY > 0 and references:
#         longest_header = max([len(get_nice_header_and_len(header, seq)) for \
#                                    header, seq in references])
#         read_count = line_count(long_reads_fastq) / 4
#     for header, seq in references:
#         nice_header = get_nice_header(header)
#         one_seq_fasta = os.path.join(working_dir, nice_header + '.fasta')
#         save_to_fasta(nice_header, seq, one_seq_fasta)
#         paf_filename = os.path.join(working_dir, nice_header + '.paf')
#         run_graphmap(one_seq_fasta, long_reads_fastq, paf_filename, graphmap_path, threads)

#         # Copy the segment's SAM alignments to the final SAM file.
#         alignment_count = 0
#         one_seq_paf = open(paf_filename, 'r')
#         for line in one_seq_paf:
#             final_paf.write(line)
#             alignment_count += 1
#         one_seq_paf.close()

#         if VERBOSITY > 0:
#             nice_header = get_nice_header_and_len(header, seq, pad_length=longest_header) + ':'
#             print(nice_header, int_to_str(alignment_count, read_count))
#             sys.stdout.flush()

#         # Clean up
#         os.remove(one_seq_fasta)
#         os.remove(paf_filename)

#     if VERBOSITY > 0:
#         print()

#     final_paf.close()

def run_graphmap(fasta, long_reads_fastq, sam_file, graphmap_path, threads, scoring_scheme):
    '''
    This function runs GraphMap for the given inputs and produces a SAM file at the given location.
    '''
    command = [graphmap_path,
               '-r', fasta,
               '-d', long_reads_fastq,
               '-o', sam_file,
               '-t', str(threads),
               '-a', 'anchorgotoh']
    command += scoring_scheme.get_graphmap_parameters()

    if VERBOSITY > 0:
        print('Running GraphMap')
        print('----------------')
        print(' '.join(command))
        print()

    # Print the GraphMap output as it comes. I gather up and display lines in a manual manner so
    # I can replace carriage returns with newlines, which makes the progress a bit cleaner.
    process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    line = ''
    while process.poll() is None:
        graphmap_output = process.stderr.read(1)
        if VERBOSITY > 0:
            line += graphmap_output
            if line.endswith('\r'):
                line = line[:-1] + '\n'
            if line.endswith('\n'):
                line = line.strip()
                if line:
                    print(line)
                line = ''
    if VERBOSITY > 0:
        print()

    # Clean up.
    if os.path.isfile(fasta + '.gmidx'):
        os.remove(fasta + '.gmidx')
    if os.path.isfile(fasta + '.gmidxsec'):
        os.remove(fasta + '.gmidxsec')

    if not os.path.isfile(sam_file):
        quit_with_error('GraphMap failure')

def load_sam_alignments(sam_filename, reads, references, scoring_scheme):
    '''
    This function returns a list of Alignment objects from the given SAM file.
    '''
    sam_alignments = []
    sam_file = open(sam_filename, 'r')
    for line in sam_file:
        line = line.strip()
        if line and not line.startswith('@') and line.split('\t', 3)[2] != '*':
            sam_alignments.append(Alignment(reads, references, scoring_scheme, sam_line=line))
    return sam_alignments

# def seqan_alignment_one_read_all_refs(reads, references, read, expected_ref_to_read_ratio):
#     '''
#     Aligns a single read against all reference sequences using Seqan. Both forward and reverse
#     complement alignments are tried.
#     Returns a list of all alignments found.
#     '''
#     output = ''
#     global VERBOSITY
#     if VERBOSITY > 1:
#         output += str(read) + '\n'
#     alignments = []

#     # For lower verbosity levels we will suppress output for the Seqan alignent.
#     verbosity_at_start = VERBOSITY
#     if VERBOSITY <= 2:
#         VERBOSITY = 0

#     for ref_name, ref_seq in references.iteritems():
#         if verbosity_at_start > 2:
#             output += 'Reference:', ref_name + '+\n'
#         forward_alignment, forward_alignment_output = \
#                         make_seqan_alignment_all_lines(reads, references, ref_name, ref_seq, read,
#                                                        False, expected_ref_to_read_ratio)
#         alignments += forward_alignment
#         output += forward_alignment_output

#         if verbosity_at_start > 2:
#             output += 'Reference:', ref_name + '-\n'
#         reverse_alignment, reverse_alignment_output = \
#                         make_seqan_alignment_all_lines(reads, references, ref_name, ref_seq, read,
#                                                        True, expected_ref_to_read_ratio)
#         alignments += reverse_alignment
#         output += reverse_alignment_output

#     if verbosity_at_start > 1:
#         if not alignments:
#             output += 'No alignments found for read ' + read.name + '\n'
#         else:
#             output += 'Alignments found for read ' + read.name + ':\n'
#             for alignment in alignments:
#                 output += '  ' + str(alignment) + '\n'
#         output += '\n'

#     # Restore the verbosity if it was suppressed.
#     if verbosity_at_start <= 2:
#         VERBOSITY = verbosity_at_start

#     return alignments

# def make_seqan_alignment_all_lines(reads, references, ref_name, ref_seq, read, rev_comp,
#                                    expected_ref_to_read_ratio):
#     '''
#     Runs an alignment using Seqan between one read and one reference.
#     Returns a list of Alignment objects: empty list means it did not succeed, a list of one means
#     it got one alignment (common) and a list of more than one means it got multiple alignments (not
#     common but possible).
#     '''
#     output = ''
#     if rev_comp:
#         read_seq = reverse_complement(read.sequence)
#     else:
#         read_seq = read.sequence

#     # Get the alignment line(s).
#     ptr = C_LIB.findAlignmentLines(read_seq, ref_seq, len(read_seq), len(ref_seq),
#                                    expected_ref_to_read_ratio, VERBOSITY)
#     line_result = cast(ptr, c_char_p).value
#     C_LIB.free_c_string(ptr)

#     line_finding_output, line_result = line_result.split(';', 1)
#     output += line_finding_output

#     if line_result.startswith('Fail'):
#         if VERBOSITY > 1:
#             output += '  No alignment lines found\n'
#         return [], output

#     # If the code got here, then we have at least one line to use in a banded alignment. Conduct
#     # an alignment for each line.
#     alignments = []
#     lines_info = line_result.split(';')
#     for line_info in lines_info:
#         line_info_parts = line_info.split(',', 3)
#         slope = float(line_info_parts[0])
#         intercept = float(line_info_parts[1])
#         k_size = int(line_info_parts[2])
#         kmer_locations = line_info_parts[3]
#         alignment, alignment_output = make_seqan_alignment_one_line(reads, references, ref_name,
#                                                                     ref_seq, read, rev_comp, slope,
#                                                                     intercept, k_size,
#                                                                     kmer_locations)
#         if alignment:
#             alignments.append(alignment)
#             output += alignment_output
#     return alignments, output

# def make_seqan_alignment_one_line(reads, references, ref_name, ref_seq, read, rev_comp,
#                                   slope, intercept, k_size, kmer_locations):
#     '''
#     Runs an alignment using Seqan between one read and one reference along one line.
#     It starts with a smallish band size (fast) and works up to larger ones to see if they improve
#     the alignment.
#     It returns either an Alignment object or None, depending on whether or not it was successful.
#     '''
#     output = ''
#     band_size = 20 # TO DO: make this a parameter?
#     max_band_size = 160 # TO DO: make this a parameter?

#     alignment, alignment_output = run_one_banded_seqan_alignment(reads, references, ref_name,
#                                                                  ref_seq, read, rev_comp,
#                                                                  band_size, slope, intercept,
#                                                                  k_size, kmer_locations)
#     output += alignment_output
#     if not alignment:
#         return None, output

#     # If our alignment succeeded, we try larger bands to see if that helps.
#     while True:
#         band_size *= 2

#         # If we've reached the max band size, then we return what we've got.
#         if band_size > max_band_size:
#             return alignment, output

#         new_alignment, alignment_output = run_one_banded_seqan_alignment(reads, references,
#                                                                          ref_name, ref_seq, read,
#                                                                          rev_comp, band_size,
#                                                                          slope, intercept, k_size,
#                                                                          kmer_locations)
#         output += alignment_output

#         # If our new alignment with a larger band size failed to improve upon our previous
#         # alignment with a smaller band size, then we don't bother trying for larger bands and just
#         # return our best alignment so far.
#         if new_alignment.scaled_score <= alignment.scaled_score:
#             return alignment, output
#         else:
#             alignment = new_alignment

# def run_one_banded_seqan_alignment(reads, references, ref_name, ref_seq, read, rev_comp,
#                                    band_size, slope, intercept, k_size, kmer_locations):
#     '''
#     Runs a single alignment using Seqan.
#     Since this does a banded alignment, it is efficient but may or may not be successful. Returns
#     either an Alignment object (if successful) or None (if not).
#     '''
#     output = ''
#     if rev_comp:
#         read_seq = reverse_complement(read.sequence)
#     else:
#         read_seq = read.sequence
#     ptr = C_LIB.bandedSemiGlobalAlignment(read_seq, ref_seq, len(read_seq), len(ref_seq), slope,
#                                           intercept, k_size, band_size, VERBOSITY,
#                                           3, -6, -5, -2, kmer_locations)
#     alignment_result = cast(ptr, c_char_p).value    
#     C_LIB.free_c_string(ptr)

#     alignment_output, alignment_result = alignment_result.split(';', 1)
#     output += alignment_output

#     if alignment_result.startswith('Failed'):
#         output += alignment_result + '\n'
#         return None, output

#     alignment = Alignment(reads, references, seqan_output=alignment_result, read_name=read.name,
#                           ref_name=ref_name, rev_comp=rev_comp)
#     if VERBOSITY > 1:
#         output += 'Seqan alignment, bandwidth = ' + str(band_size) + ': ' + str(alignment) + '\n'
#     if VERBOSITY > 2:
#         output += alignment.cigar + '\n'
#     return alignment, output

def get_ref_shift_from_cigar_part(cigar_part):
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

# def simplify_ranges(ranges):
#     '''
#     Collapses overlapping ranges together. Input ranges are tuples of (start, end) in the normal
#     Python manner where the end isn't included.
#     '''
#     fixed_ranges = []
#     for int_range in ranges:
#         if int_range[0] > int_range[1]:
#             fixed_ranges.append((int_range[1], int_range[0]))
#         elif int_range[0] < int_range[1]:
#             fixed_ranges.append(int_range)
#     starts_ends = [(x[0], 1) for x in fixed_ranges]
#     starts_ends += [(x[1], -1) for x in fixed_ranges]
#     starts_ends.sort(key=lambda x: x[0])
#     current_sum = 0
#     cumulative_sum = []
#     for start_end in starts_ends:
#         current_sum += start_end[1]
#         cumulative_sum.append((start_end[0], current_sum))
#     prev_depth = 0
#     start = 0
#     combined = []
#     for pos, depth in cumulative_sum:
#         if prev_depth == 0:
#             start = pos
#         elif depth == 0:
#             combined.append((start, pos))
#         prev_depth = depth
#     return combined

# def range_is_contained(test_range, other_ranges):
#     '''
#     Returns True if test_range is entirely contained within any range in other_ranges.
#     '''
#     start, end = test_range
#     for other_range in other_ranges:
#         if other_range[0] <= start and other_range[1] >= end:
#             return True
#     return False

def write_sam_file(alignments, graphmap_sam, sam_filename):
    '''
    Writes the given alignments to a SAM file.
    '''
    sam_file = open(sam_filename, 'w')

    # Copy the header from the GraphMap SAM file.
    graphmap_sam_file = open(graphmap_sam, 'r')
    for line in graphmap_sam_file:
        if line.startswith('@'):
            sam_file.write(line)
        else:
            break

    for alignment in alignments:
        sam_file.write(alignment.get_sam_line())
        sam_file.write('\n')
    sam_file.close()

# def load_fasta(filename): # type: (str) -> list[tuple[str, str]]
#     '''
#     Returns the names and sequences for the given fasta file.
#     '''
#     fasta_seqs = []
#     fasta_file = open(filename, 'r')
#     name = ''
#     sequence = ''
#     for line in fasta_file:
#         line = line.strip()
#         if not line:
#             continue
#         if line[0] == '>': # Header line = start of new contig
#             if name:
#                 fasta_seqs.append((name.split()[0], sequence))
#                 name = ''
#                 sequence = ''
#             name = line[1:]
#         else:
#             sequence += line
#     if name:
#         fasta_seqs.append((name.split()[0], sequence))
#     return fasta_seqs

def is_header_spades_format(contig_name):
    '''
    Returns whether or not the header appears to be in the SPAdes/Velvet format.
    Example: NODE_5_length_150905_cov_4.42519
    '''
    contig_name_parts = contig_name.split('_')
    return len(contig_name_parts) > 5 and \
           (contig_name_parts[0] == 'NODE' or contig_name_parts[0] == 'EDGE') and \
           contig_name_parts[2] == 'length' and contig_name_parts[4] == 'cov'

def get_nice_header(header):
    '''
    For a header with a SPAdes/Velvet format, this function returns a simplified string that is
    just NODE_XX where XX is the contig number.
    For any other format, this function trims off everything following the first whitespace.
    '''
    if is_header_spades_format(header):
        return 'NODE_' + header.split('_')[1]
    else:
        return header.split()[0]

# def get_nice_header_and_len(header, seq, pad_length=0):
#     '''
#     Add the length in base pairs to the nice header. If there is a pad length, it will add spaces
#     in between the header and length so things can line up nicely.
#     '''
#     part_1 = get_nice_header(header) + ' '
#     part_2 = '(' + '{:,}'.format(len(seq)) + ' bp)'
#     if len(part_1) + len(part_2) < pad_length:
#         spaces = ' ' * (pad_length - len(part_1) - len(part_2))
#     else:
#         spaces = ''
#     return part_1 + spaces + part_2

# def save_to_fasta(header, sequence, filename):
#     '''
#     Saves the header/sequence to FASTA file.
#     '''
#     fasta = open(filename, 'w')
#     fasta.write('>' + header + '\n')
#     fasta.write(add_line_breaks_to_sequence(sequence, 60))
#     fasta.close()

# def save_reads_to_fastq(reads, fastq_filename):
#     '''
#     Writes the given reads to a FASTQ file.
#     '''
#     fastq = open(fastq_filename, 'w')
#     for read in reads:
#         fastq.write(read.get_fastq())

# def add_line_breaks_to_sequence(sequence, length):
#     '''
#     Wraps sequences to the defined length.  All resulting sequences end in a line break.
#     '''
#     seq_with_breaks = ''
#     while len(sequence) > length:
#         seq_with_breaks += sequence[:length] + '\n'
#         sequence = sequence[length:]
#     if len(sequence) > 0:
#         seq_with_breaks += sequence
#         seq_with_breaks += '\n'
#     return seq_with_breaks

def check_file_exists(filename): # type: (str) -> bool
    '''
    Checks to make sure the single given file exists.
    '''
    if not os.path.isfile(filename):
        quit_with_error('could not find ' + filename)

def quit_with_error(message): # type: (str) -> None
    '''
    Displays the given message and ends the program's execution.
    '''
    print('Error:', message, file=sys.stderr)
    sys.exit(1)

def float_to_str(num, decimals, max_num=0):
    '''
    Converts a number to a string. Will add left padding based on the max value to ensure numbers
    align well.
    '''
    num_str = '%.' + str(decimals) + 'f'
    num_str = num_str % num
    after_decimal = num_str.split('.')[1]
    num_str = int_to_str(int(num)) + '.' + after_decimal
    if max_num > 0:
        max_str = float_to_str(max_num, decimals)
        num_str = num_str.rjust(len(max_str))
    return num_str

def int_to_str(num, max_num=0):
    '''
    Converts a number to a string. Will add left padding based on the max value to ensure numbers
    align well.
    '''
    num_str = '{:,}'.format(num)
    max_str = '{:,}'.format(int(max_num))
    return num_str.rjust(len(max_str))

# def line_count(filename):
#     '''
#     Counts the lines in the given file.
#     '''
#     i = 0
#     with open(filename) as file_to_count:
#         for i, _ in enumerate(file_to_count):
#             pass
#     return i + 1

def clean_str_for_filename(filename):
    '''
    This function removes characters from a string which would not be suitable in a filename.
    It also turns spaces into underscores, because filenames with spaces can occasionally cause
    issues.
    http://stackoverflow.com/questions/295135/turn-a-string-into-a-valid-filename-in-python
    '''
    valid_chars = "-_.() %s%s" % (string.ascii_letters, string.digits)
    filename_valid_chars = ''.join(c for c in filename if c in valid_chars)
    return filename_valid_chars.replace(' ', '_')

def reverse_complement(seq):
    '''Given a DNA sequences, this function returns the reverse complement sequence.'''
    rev_comp = ''
    for i in reversed(range(len(seq))):
        rev_comp += complement_base(seq[i])
    return rev_comp

def complement_base(base):
    '''Given a DNA base, this returns the complement.'''
    forward = 'ATGCatgcRYSWKMryswkmBDHVbdhvNn.-?'
    reverse = 'TACGtacgYRSWMKyrswmkVHDBvhdbNn.-?N'
    return reverse[forward.find(base)]


# def seqan_from_paf_alignment_one_arg(all_args):
#     '''
#     This function is just a convenience for calling seqan_from_paf_alignment with a single argument
#     (makes multithreading easier).
#     '''
#     paf_alignment, reads, references, expected_ref_to_read_ratio = all_args
#     return seqan_from_paf_alignment(paf_alignment, reads, references, expected_ref_to_read_ratio)

# def seqan_from_paf_alignment(paf_alignment, reads, references, expected_ref_to_read_ratio):
#     '''
#     This function is used on PAF alignments to turn them into Seqan alignments. This takes them
#     from basic alignments with only approximate start/end positions into full alignments.
#     '''
#     output = ''
#     assert paf_alignment.alignment_type == 'PAF'
#     if VERBOSITY > 1:
#         output += 'PAF alignment before Seqan: ' + str(paf_alignment) + '\n'

#     # We can use the PAF alignment's range to trim the sequences down to the alignment region.
#     # This is especially useful if the reference sequence is very large (like a whole bacterial
#     # genome).
#     # TO DO
#     # TO DO
#     # TO DO
#     # TO DO
#     # TO DO
#     # TO DO
#     # TO DO
#     # TO DO
#     # TO DO
#     # TO DO
#     # TO DO
#     # TO DO
#     # TO DO
#     # TO DO
#     # TO DO

#     alignments, alignment_output = make_seqan_alignment_all_lines(reads, references,
#                                                                   paf_alignment.ref_name,
#                                                                   paf_alignment.full_ref_sequence,
#                                                                   paf_alignment.read,
#                                                                   paf_alignment.rev_comp,
#                                                                   expected_ref_to_read_ratio)
#     output += alignment_output
#     if VERBOSITY > 1:
#         output += '\n'
#     return alignments, output

def get_median(num_list):
    '''
    Returns the median of the given list of numbers.
    '''
    count = len(num_list)
    if count == 0:
        return 0.0
    sorted_list = sorted(num_list)
    if count % 2 == 0:
        return (sorted_list[count // 2 - 1] + sorted_list[count // 2]) / 2.0
    else:
        return sorted_list[count // 2]

def get_median_and_mad(num_list):
    '''
    Returns the median and MAD of the given list of numbers.
    '''
    median = get_median(num_list)
    absolute_deviations = [abs(x - median) for x in num_list]
    mad = 1.4826 * get_median(absolute_deviations)
    return median, mad


def get_mean_and_st_dev(num_list):
    '''
    This function returns the mean and standard deviation of the given list of numbers.
    '''
    num = len(num_list)
    if num == 0:
        return 0.0, 0.0
    mean = sum(num_list) / num
    if num == 1:
        return mean, -1
    sum_squares = sum((x - mean) ** 2 for x in num_list)
    st_dev = (sum_squares / (num - 1)) ** 0.5
    return mean, st_dev

def print_progress_line(completed, total, base_pairs=None):
    '''
    Prints a progress line to the screen using a carriage return to overwrite the previous progress
    line.
    '''
    progress_str = int_to_str(completed) + ' / ' + int_to_str(total)
    progress_str += ' (' + '%.1f' % (100.0 * completed / total) + '%)'
    if base_pairs is not None:
        progress_str += ' - ' + int_to_str(base_pairs) + ' bp'
    print('\r' + progress_str, end='')
    sys.stdout.flush()

# def group_reads_by_fraction_aligned(reads):
#     '''
#     Groups reads into three groups:
#       1) Fully aligned
#       2) Partially aligned
#       3) Unaligned
#     '''
#     fully_aligned_reads = []
#     partially_aligned_reads = []
#     unaligned_reads = []
#     for read in reads.itervalues():
#         fraction_aligned = read.get_fraction_aligned()
#         if fraction_aligned == 1.0:
#             fully_aligned_reads.append(read)
#         elif fraction_aligned == 0.0:
#             unaligned_reads.append(read)
#         else:
#             partially_aligned_reads.append(read)
#     return fully_aligned_reads, partially_aligned_reads, unaligned_reads



class AlignmentScoringScheme(object):
    '''
    This class holds an alignment scoring scheme.
    '''
    def __init__(self, scheme_string):
        scheme_parts = scheme_string.split(',')

        # Default scoring scheme
        self.match = 3
        self.mismatch = -6
        self.gap_open = -5
        self.gap_extend = -2

        if len(scheme_parts) == 4:
            self.match = int(scheme_parts[0])
            self.mismatch = int(scheme_parts[1])
            self.gap_open = int(scheme_parts[2])
            self.gap_extend = int(scheme_parts[3])

        assert self.match > 0
        assert self.mismatch < 0
        assert self.gap_open < 0
        assert self.gap_extend < 0

    def get_graphmap_parameters(self):
        '''
        Returns the scoring scheme in the form of GraphMap parameters for subprocess.
        '''
        return ['-M', str(self.match),
                '-X', str(-self.mismatch),
                '-G', str(-self.gap_open),
                '-E', str(-self.gap_extend)]



class LongRead(object):
    '''
    This class holds a long read, e.g. from PacBio or Oxford Nanopore.
    '''
    def __init__(self, name, sequence, qualities):
        self.name = name
        self.sequence = sequence
        self.rev_comp_seq = reverse_complement(sequence)
        self.qualities = qualities
        self.alignments = []

#     def __repr__(self):
#         return self.name + ' (' + str(len(self.sequence)) + ' bp)'

    def get_length(self):
        '''
        Returns the sequence length.
        '''
        return len(self.sequence)

#     def remove_conflicting_alignments(self):
#         '''
#         This function removes alignments from the read which are likely to be spurious. It sorts
#         alignments by identity and works through them from highest identity to lowest identity,
#         only keeping alignments that cover new parts of the read.
#         It also uses an identity threshold to remove alignments with very poor identity.
#         '''
#         self.alignments = sorted(self.alignments, reverse=True,
#                                  key=lambda x: (x.percent_identity, random.random()))
#         kept_alignments = []
#         read_ranges = []
#         for alignment in self.alignments:
#             read_range = alignment.read_start_end_positive_strand()
#             if not range_is_contained(read_range, read_ranges):
#                 read_ranges.append(read_range)
#                 read_ranges = simplify_ranges(read_ranges)
#                 kept_alignments.append(alignment)
#         self.alignments = kept_alignments

#     def remove_low_id_alignments(self, id_threshold):
#         '''
#         This function removes alignments with identity below the cutoff.
#         '''
#         self.alignments = [x for x in self.alignments if x.percent_identity >= id_threshold]

#     def get_fastq(self):
#         '''
#         Returns a string for the read in FASTQ format. It contains four lines and ends in a line
#         break.
#         '''
#         return '@' + self.name + '\n' + \
#                self.sequence + '\n' + \
#                '+' + self.name + '\n' + \
#                self.qualities + '\n'

#     def get_descriptive_string(self):
#         '''
#         Returns a multi-line string that describes the read and its alignments.
#         '''
#         header = self.name + ' (' + str(len(self.sequence)) + ' bp)'
#         line = '-' * len(header)
#         description = header + '\n' + line + '\n'
#         if not self.alignments:
#             description += 'no alignments'
#         else:
#             description += '%.2f' % (100.0 * self.get_fraction_aligned()) + '% aligned\n'
#             description += '\n'.join([str(x) for x in self.alignments])
#         return description + '\n\n'

#     def get_fraction_aligned(self):
#         '''
#         This function returns the fraction of the read which is covered by any of the read's
#         alignments.
#         '''
#         read_ranges = [x.read_start_end_positive_strand() for x in self.alignments]
#         read_ranges = simplify_ranges(read_ranges)
#         aligned_length = sum([x[1] - x[0] for x in read_ranges])
#         return aligned_length / len(self.sequence)


class Alignment(object):
    '''
    This class describes an alignment between a long read and a contig.
    It can either describe a full Alignment generated by Seqan or a basic alignment from a PAF
    file made by GraphMap in owler mode.
    Both types require a dictionary of reads and a dictionary of references.
    To make a PAF alignment, the only additional requirement is the PAF line.
    To make a Seqan alignment, the additional requirements are the seqan output, read name,
    reference name and whether the alignment was on the reverse complement strand.
    '''
    def __init__(self, reads, references, scoring_scheme, sam_line=None, seqan_output=None,
                 read_name=None, ref_name=None, rev_comp=None):

        assert sam_line or (seqan_output and read_name and ref_name)

        # Read details
        self.read = None
        self.read_start_pos = None
        self.read_end_pos = None
        self.read_end_gap = None
        self.aligned_read_seq = None

        # Reference details
        self.ref_name = None
        self.full_ref_sequence = None
        self.ref_start_pos = None
        self.ref_end_pos = None
        self.ref_end_gap = None
        self.aligned_ref_seq = None

        # Alignment details
        self.alignment_type = None
        self.rev_comp = None
        self.cigar = None
        self.cigar_parts = None
        self.match_count = None
        self.mismatch_count = None
        self.insertion_count = None
        self.deletion_count = None
        self.alignment_length = None
        self.edit_distance = None
        self.percent_identity = None
        self.ref_mismatch_positions = None
        self.ref_insertion_positions_and_sizes = None
        self.ref_deletion_positions = None
        self.raw_score = None
        self.scaled_score = None
        self.milliseconds = None

        if seqan_output and read_name and ref_name:
            self.setup_using_seqan_output(seqan_output, reads, references, read_name, ref_name,
                                          rev_comp)
        elif sam_line:
            self.setup_using_graphmap_sam(sam_line, reads, references, scoring_scheme)


    # def setup_using_seqan_output(self, seqan_output, reads, references, read_name, ref_name,
    #                              rev_comp):
    #     '''
    #     This function sets up the Alignment using the Seqan results. This kind of alignment has
    #     complete details about the alignment.
    #     '''
    #     self.alignment_type = 'Seqan'
    #     seqan_parts = seqan_output.split(';')
    #     assert len(seqan_parts) >= 17

    #     # Read details
    #     self.read = reads[read_name]
    #     self.read_start_pos = int(seqan_parts[1])
    #     self.read_end_pos = int(seqan_parts[2])
    #     self.read_end_gap = self.read.get_length() - self.read_end_pos
    #     self.aligned_read_seq = self.read.sequence[self.read_start_pos:self.read_end_pos]

    #     # Reference details
    #     self.ref_name = ref_name
    #     self.full_ref_sequence = references[self.ref_name]
    #     self.ref_start_pos = int(seqan_parts[3])
    #     self.ref_end_pos = int(seqan_parts[4])
    #     self.ref_end_gap = len(self.full_ref_sequence) - self.ref_end_pos
    #     self.aligned_ref_seq = self.full_ref_sequence[self.ref_start_pos:self.ref_end_pos]

    #     # Alignment details
    #     self.rev_comp = rev_comp
    #     self.cigar = seqan_parts[0]
    #     self.cigar_parts = re.findall(r'\d+\w', self.cigar)
    #     self.match_count = int(seqan_parts[6])
    #     self.mismatch_count = int(seqan_parts[7])
    #     self.insertion_count = int(seqan_parts[9])
    #     self.deletion_count = int(seqan_parts[11])
    #     self.alignment_length = int(seqan_parts[5])
    #     self.edit_distance = int(seqan_parts[13])
    #     self.percent_identity = float(seqan_parts[14])
    #     self.ref_mismatch_positions = []
    #     self.ref_insertion_positions = []
    #     self.ref_deletion_positions = []
    #     if seqan_parts[8]:
    #         self.ref_mismatch_positions = [int(x) for x in seqan_parts[8].split(',')]
    #     if seqan_parts[10]:
    #         self.ref_insertion_positions = [int(x) for x in seqan_parts[10].split(',')]
    #     if seqan_parts[12]:
    #         self.ref_deletion_positions = [int(x) for x in seqan_parts[12].split(',')]
    #     self.raw_score = int(seqan_parts[15])
    #     self.scaled_score = float(seqan_parts[16])
    #     self.milliseconds = int(seqan_parts[17])

    def setup_using_graphmap_sam(self, sam_line, reads, references, scoring_scheme):
        '''
        This function sets up the Alignment using a SAM line.
        '''
        self.alignment_type = 'GraphMap'
        sam_parts = sam_line.split('\t')
        self.rev_comp = bool(int(sam_parts[1]) & 0x10)

        self.cigar = sam_parts[5]
        self.cigar_parts = re.findall(r'\d+\w', self.cigar)

        self.read = reads[sam_parts[0]]
        self.read_start_pos = self.get_start_soft_clips()
        self.read_end_pos = self.read.get_length() - self.get_end_soft_clips()
        self.read_end_gap = self.get_end_soft_clips()
        if self.rev_comp:
            self.aligned_read_seq = self.read.rev_comp_seq[self.read_start_pos:self.read_end_pos]
        else:
            self.aligned_read_seq = self.read.sequence[self.read_start_pos:self.read_end_pos]

        self.ref_name = get_nice_header(sam_parts[2])
        self.full_ref_sequence = references[self.ref_name]
        self.ref_start_pos = int(sam_parts[3]) - 1
        self.ref_end_pos = self.ref_start_pos
        for cigar_part in self.cigar_parts:
            self.ref_end_pos += get_ref_shift_from_cigar_part(cigar_part)
        self.ref_end_gap = len(self.full_ref_sequence) - self.ref_end_pos
        self.aligned_ref_seq = self.full_ref_sequence[self.ref_start_pos:self.ref_end_pos]

        self.tally_up_score_and_errors(scoring_scheme)


    def tally_up_score_and_errors(self, scoring_scheme):
        '''
        This function steps through the CIGAR string for the alignment to get the score, identity
        and count/locations of errors.
        '''
        # Clear any existing tallies.
        self.match_count = 0
        self.mismatch_count = 0
        self.insertion_count = 0
        self.deletion_count = 0
        self.percent_identity = 0.0
        self.ref_mismatch_positions = []
        self.ref_insertion_positions_and_sizes = []
        self.ref_deletion_positions = []
        
        # Remove the soft clipping parts of the CIGAR string for tallying.
        cigar_parts = self.cigar_parts[:]
        if cigar_parts[0][-1] == 'S':
            cigar_parts.pop(0)
        if cigar_parts[-1][-1] == 'S':
            cigar_parts.pop()

        read_i = 0
        ref_i = 0
        align_i = 0
        self.raw_score = 0
        for cigar_part in cigar_parts:
            cigar_count = int(cigar_part[:-1])
            cigar_type = cigar_part[-1]
            if cigar_type == 'I':
                self.insertion_count += cigar_count
                self.ref_insertion_positions_and_sizes.append((ref_i + self.ref_start_pos,
                                                               cigar_count))
                read_i += cigar_count
                self.raw_score += scoring_scheme.gap_open + \
                                  ((cigar_count - 1) * scoring_scheme.gap_extend)
            elif cigar_type == 'D':
                self.deletion_count += cigar_count
                for i in xrange(cigar_count):
                    self.ref_deletion_positions.append(ref_i + self.ref_start_pos + i)
                ref_i += cigar_count
                self.raw_score += scoring_scheme.gap_open + \
                                  ((cigar_count - 1) * scoring_scheme.gap_extend)
            else: # match/mismatch
                for _ in xrange(cigar_count):
                    read_base = self.aligned_read_seq[read_i]
                    ref_base = self.aligned_ref_seq[ref_i]
                    if read_base == ref_base:
                        self.match_count += 1
                        self.raw_score += scoring_scheme.match
                    else:
                        self.mismatch_count += 1
                        self.ref_mismatch_positions.append(ref_i + self.ref_start_pos)
                        self.raw_score += scoring_scheme.mismatch
                    read_i += 1
                    ref_i += 1
            align_i += cigar_count
        self.percent_identity = 100.0 * self.match_count / align_i
        self.edit_distance = self.mismatch_count + self.insertion_count + self.deletion_count
        self.alignment_length = align_i
        perfect_score = scoring_scheme.match * self.alignment_length
        self.scaled_score = 100.0 * self.raw_score / perfect_score

    # def setup_using_paf_line(self, paf_line, reads, references):
    #     self.alignment_type = 'PAF'
    #     paf_parts = paf_line.split()

    #     # Read details
    #     self.read = reads[paf_parts[0]]
    #     self.read_start_pos = int(paf_parts[2])
    #     self.read_end_pos = int(paf_parts[3])
    #     self.read_end_gap = self.read.get_length() - self.read_end_pos

    #     # Reference details
    #     self.ref_name = paf_parts[5]
    #     self.full_ref_sequence = references[self.ref_name]
    #     self.ref_start_pos = int(paf_parts[7])
    #     self.ref_end_pos = int(paf_parts[8])
    #     self.ref_end_gap = len(self.full_ref_sequence) - self.ref_end_pos

    #     # Alignment details
    #     self.rev_comp = (paf_parts[4] == '-')
    #     self.alignment_length = int(paf_parts[10])

    #     # Extend the alignment to reach either the end of the read or reference, whichever comes
    #     # first.
    #     slope = (self.ref_end_pos - self.ref_start_pos) / (self.read_end_pos - self.read_start_pos)
    #     intercept = self.ref_start_pos - (slope * self.read_start_pos)
    #     missing_bases_at_start = min(self.read_start_pos, self.ref_start_pos)
    #     missing_bases_at_end = min(self.read_end_gap, self.ref_end_gap)
    #     old_read_start = self.read_start_pos
    #     old_ref_start = self.ref_start_pos
    #     old_read_end = self.read_end_pos
    #     old_ref_end = self.ref_end_pos
    #     if missing_bases_at_start:
    #         if intercept >= 0:
    #             self.read_start_pos = 0
    #             self.ref_start_pos = int(round(intercept))
    #         else:
    #             self.read_start_pos = int(round(-intercept / slope))
    #             self.ref_start_pos = 0
    #     if missing_bases_at_end:
    #         read_end_intercept = (slope * len(self.read.sequence)) + intercept
    #         if read_end_intercept <= len(self.full_ref_sequence):
    #             self.read_end_pos = len(self.read.sequence)
    #             self.ref_end_pos = int(round(read_end_intercept))
    #         else:
    #             self.read_end_pos = int(round((len(self.full_ref_sequence) - intercept) / slope))
    #             self.ref_end_pos = len(self.full_ref_sequence)
    #     self.read_end_gap = self.read.get_length() - self.read_end_pos
    #     self.ref_end_gap = len(self.full_ref_sequence) - self.ref_end_pos
    #     self.extension_length = (old_read_start - self.read_start_pos) + \
    #                             (old_ref_start - self.ref_start_pos) + \
    #                             (self.read_end_pos - old_read_end) + \
    #                             (self.ref_end_pos - old_ref_end)

    #     # The alignment extension should not have made much of a change to the slope.
    #     # TEMP CHECKING CODE - CAN BE REMOVED LATER
    #     new_slope = (self.ref_end_pos - self.ref_start_pos) / \
    #                 (self.read_end_pos - self.read_start_pos)
    #     assert new_slope / slope > 0.995
    #     assert new_slope / slope < 1.005

    def extend_start(self, scoring_scheme):
        '''
        This function extends the start of the alignment to remove any missing start bases.
        '''
        missing_start_bases = self.get_missing_bases_at_start()
        realigned_bases = 2 * missing_start_bases

        realigned_read_end = self.read_start_pos
        realigned_read_start = max(0, realigned_read_end - realigned_bases)

        realigned_ref_end = self.ref_start_pos
        realigned_ref_start = max(0, realigned_ref_end - realigned_bases)

        if self.rev_comp:
            realigned_read_seq = \
                reverse_complement(self.read.sequence)[realigned_read_start:realigned_read_end]
        else:
            realigned_read_seq = self.read.sequence[realigned_read_start:realigned_read_end]
        realigned_ref_seq = self.full_ref_sequence[realigned_ref_start:realigned_ref_end]

        assert len(realigned_ref_seq) >= len(realigned_read_seq)




        print('BEFORE')
        print('------')
        print(self)
        print('CIGAR:', self.cigar[:50] + '...')
        print('missing_start_bases:', missing_start_bases)
        print('realigned_bases:', realigned_bases)
        print('realigned_read_start:', realigned_read_start)
        print('realigned_read_end:', realigned_read_end)
        print('realigned_ref_start:', realigned_ref_start)
        print('realigned_ref_end:', realigned_ref_end)
        print('realigned_read_seq:', realigned_read_seq)
        print('realigned_ref_seq:', realigned_ref_seq)
        print('')






        # Call the C++ function to do the actual alignment.
        ptr = C_LIB.startExtensionAlignment(realigned_read_seq, realigned_ref_seq,
                                            len(realigned_read_seq), len(realigned_ref_seq),
                                            VERBOSITY, scoring_scheme.match,
                                            scoring_scheme.mismatch, scoring_scheme.gap_open,
                                            scoring_scheme.gap_extend)
        alignment_result = cast(ptr, c_char_p).value    
        C_LIB.free_c_string(ptr)

        seqan_parts = alignment_result.split(';')
        assert len(seqan_parts) >= 18


        print('alignment_result:', alignment_result)
        print('')

        # Set the new read start.
        self.read_start_pos = int(seqan_parts[2])
        assert self.read_start_pos == 0
        if self.rev_comp:
            self.aligned_read_seq = self.read.rev_comp_seq[self.read_start_pos:self.read_end_pos]
        else:
            self.aligned_read_seq = self.read.sequence[self.read_start_pos:self.read_end_pos]

        # Set the new reference start
        self.ref_start_pos = realigned_ref_start + int(seqan_parts[4])
        self.aligned_ref_seq = self.full_ref_sequence[self.ref_start_pos:self.ref_end_pos]

        # Replace the S part at the beginning the alignment's CIGAR with the CIGAR just made. If
        # the last part of the new CIGAR is of the same type as the first part of the existing
        # CIGAR, they will need to be merged.
        new_cigar_parts = re.findall(r'\d+\w', seqan_parts[1])
        old_cigar_parts = self.cigar_parts[1:]
        if new_cigar_parts[-1][-1] == old_cigar_parts[0][-1]:
            part_sum = int(new_cigar_parts[-1][:-1]) + int(old_cigar_parts[0][:-1])
            merged_part = str(part_sum) + new_cigar_parts[-1][-1]
            new_cigar_parts = new_cigar_parts[:-1] + [merged_part]
            old_cigar_parts = old_cigar_parts[1:]
        self.cigar_parts = new_cigar_parts + old_cigar_parts
        self.cigar = ''.join(self.cigar_parts)

        self.tally_up_score_and_errors(scoring_scheme)


        print('AFTER')
        print('-----')
        print(self)
        print('CIGAR:', self.cigar[:50] + '...')
        print('missing_start_bases:', self.get_missing_bases_at_start())
        print('\n\n\n\n')

    def extend_end(self, scoring_scheme):
        '''
        This function extends the end of the alignment to remove any missing end bases.
        '''
        # TO DO
        # TO DO
        # TO DO
        # TO DO
        # TO DO
        # TO DO
        # TO DO
        # TO DO
        # TO DO
        # TO DO
        # TO DO
        # TO DO

    def __repr__(self):
        read_start, read_end = self.read_start_end_positive_strand()
        return_str = self.read.name + ' (' + str(read_start) + '-' + str(read_end) + ', '
        if self.rev_comp:
            return_str += 'strand: -), '
        else:
            return_str += 'strand: +), '
        return_str += self.ref_name + ' (' + str(self.ref_start_pos) + '-' + \
                      str(self.ref_end_pos) + ')'
        return_str += ', ' + '%.2f' % self.percent_identity + '% ID'
        return_str += ', score = ' + '%.2f' % self.scaled_score
        return_str += ', longest indel: ' + str(self.get_longest_indel_run())
        if self.alignment_type == 'Seqan':
            return_str += ', ' + str(self.milliseconds) + ' ms'
        return return_str

    def get_ref_to_read_ratio(self):
        '''
        Returns the length ratio between the aligned parts of the reference and read.
        '''
        return (self.ref_end_pos - self.ref_start_pos) / (self.read_end_pos - self.read_start_pos)

    def get_read_to_ref_ratio(self):
        '''
        Returns the length ratio between the aligned parts of the read and reference.
        '''
        return 1.0 / self.get_ref_to_read_ratio()

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

    def get_sam_line(self):
        '''
        Returns a SAM alignment line.
        '''
        edit_distance = self.mismatch_count + self.insertion_count + self.deletion_count
        sam_flag = 0 #TEMP
        return '\t'.join([self.read.name, str(sam_flag), self.ref_name,
                          str(self.ref_start_pos + 1), '255', self.cigar,
                          '*', '0', '0', self.read.sequence, self.read.qualities,
                          'NM:i:' + str(edit_distance), 'AS:i:' + str(self.raw_score)])

    def is_whole_read(self):
        '''
        Returns True if the alignment covers the entirety of the read.
        '''
        return self.read_start_pos == 0 and self.read_end_gap == 0

    def get_longest_indel_run(self):
        '''
        Returns the longest indel in the alignment.
        '''
        longest_indel_run = 0
        for cigar_part in self.cigar_parts:
            cigar_type = cigar_part[-1]
            if cigar_type == 'I' or cigar_type == 'D':
                longest_indel_run = max(longest_indel_run, int(cigar_part[:-1]))
        return longest_indel_run

    def get_missing_bases_at_start(self):
        '''
        Returns the number of bases at the start of the alignment which are missing in both the
        read and the reference (preventing the alignment from being semi-global).
        '''
        return min(self.read_start_pos, self.ref_start_pos)

    def get_missing_bases_at_end(self):
        '''
        Returns the number of bases at the end of the alignment which are missing in both the read
        and the reference (preventing the alignment from being semi-global).
        '''
        return min(self.read_end_gap, self.ref_end_gap)

    def get_total_missing_bases(self):
        '''
        Returns the number of bases at the start and end of the alignment which are missing in both
        the read and the reference (preventing the alignment from being semi-global).
        '''
        return self.get_missing_bases_at_start() + self.get_missing_bases_at_end()


if __name__ == '__main__':
    main()
