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
import random
import argparse
from multiprocessing.dummy import Pool as ThreadPool
from multiprocessing import cpu_count

'''
VERBOSITY controls how much the script prints to the screen.
'''
VERBOSITY = 0


def main():
    '''
    If this script is run on its own, execution starts here.
    '''
    full_command = ' '.join(sys.argv)
    args = get_arguments()

    # LOAD IN THE REFERENCES FROM THE SAM HEADER. JUST NEED NAME AND LENGTH.

    # FOR EACH REFERENCE, CREATE TWO LISTS: DEPTH AND ERROR COUNT

    # LOAD IN THE ALIGNMENTS FROM THE SAM, ADDING ERRORS

    summarise_errors(references, reads, args.table)

    sys.exit(0)

def get_arguments():
    '''
    Specifies the command line arguments required by the script.
    '''
    parser = argparse.ArgumentParser(description='Long read assembly checker',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('--sam', type=str, required=True, default=argparse.SUPPRESS,
                        help='Input SAM file of alignments')
    parser.add_argument('--window', type=int, required=False, default=1000,
                        help='Window size for error summaries')
    parser.add_argument('--tables', type=str, required=False, default=argparse.SUPPRESS,
                        help='Path and/or prefix for table files summarising reference errors at '
                             'each base (default: do not save base tables)')
    parser.add_argument('--html', type=str, required=False, default=argparse.SUPPRESS,
                        help='Path and/or prefix for html files with plots (default: do not save '
                             'html files)')
    parser.add_argument('--base_tables', type=str, required=False, default=argparse.SUPPRESS,
                        help='Path and/or prefix for table files summarising reference errors at '
                             'each base (default: do not save base tables)')
    parser.add_argument('--threads', type=int, required=False, default=argparse.SUPPRESS,
                        help='Number of CPU threads used to align (default: the number of '
                             'available CPUs)')
    parser.add_argument('--verbosity', type=int, required=False, default=1,
                        help='Level of stdout information (0 to 4)')

    args = parser.parse_args()

    global VERBOSITY
    VERBOSITY = args.verbosity

    # If some arguments weren't set, set them to None/False. We don't use None/False as a default
    # in add_argument because it makes the help text look weird.
    try:
        args.table
    except AttributeError:
        args.table = None
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

def summarise_errors(references, long_reads, table_prefix):
    '''
    Writes a table file summarising the alignment errors in terms of reference sequence position.
    Works in a brute force manner - could be made more efficient later if necessary.
    '''
    # If we are not making table files or printing a summary, then quit now because there's nothing
    # else to do.
    if not table_prefix and VERBOSITY == 0:
        return

    # # We are be willing to throw out the worst alignments for any particular location. This value
    # # specifies how many of our alignments we'll keep.
    # # E.g. if 0.75, we'll throw out up to 25% of the alignments for each reference location.
    # frac_alignments_to_keep = percent_alignments_to_keep / 100.0

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
                    max([len(x.sequence) for x in references]))
        print('Alignment summaries per reference')
        print('---------------------------------')

    for reference in references:
        # Create a table file for the reference.
        if table_prefix:
            header_for_filename = clean_str_for_filename(reference.name)
            if table_prefix.endswith('/'):
                table_filename = table_prefix + header_for_filename + '.txt'
            else:
                table_filename = table_prefix + '_' + header_for_filename + '.txt'
            table = open(table_filename, 'w')
            table.write('\t'.join(['base',
                                   'read depth',
                                   'mismatches',
                                   'deletions',
                                   'insertions',
                                   'mismatch rate',
                                   'deletion rate',
                                   'insertion rate',
                                   'insertion sizes']) + '\n')

        # Gather up the alignments for this reference and count up the depth for each reference
        # position.
        seq_len = len(reference.sequence)
        depths = [0] * seq_len
        alignments = []
        for read in long_reads.itervalues():
            for alignment in read.alignments:
                if alignment.ref.name == reference.name:
                    alignments.append(alignment)
                    for pos in xrange(alignment.ref_start_pos, alignment.ref_end_pos):
                        depths[pos] += 1

        # # Discard as many of the worst alignments as we can while keeping a sufficient read depth
        # # at each position.
        # sufficient_depths = [math.ceil(frac_alignments_to_keep * x) for x in total_depths]
        # filtered_depths = total_depths[:]
        # filtered_alignments = []
        # alignments = sorted(alignments, key=lambda x: x.scaled_score)
        # for alignment in alignments:
        #     removal_okay = True
        #     for pos in xrange(alignment.ref_start_pos, alignment.ref_end_pos):
        #         if filtered_depths[pos] - 1 < sufficient_depths[pos]:
        #             removal_okay = False
        #             break
        #     if removal_okay:
        #         for pos in xrange(alignment.ref_start_pos, alignment.ref_end_pos):
        #             filtered_depths[pos] -= 1
        #     else:
        #         filtered_alignments.append(alignment)

        # Determine how much of the reference sequence has at least one aligned read.
        aligned_len = sum([1 for x in depths if x > 0])
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

        for alignment in alignments:
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
            if depths[i] > 0:
                mismatch_rate = mismatches[i] / depths[i]
                deletion_rate = deletions[i] / depths[i]
                insertion_rate = insertions[i] / depths[i]
                if i in insertion_sizes:
                    insertion_rate_multi_base = sum(insertion_sizes[i]) / depths[i]
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
                                       str(depths[i]),
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
            for alignment in alignments:
                if alignment.is_whole_read():
                    contained_alignment_count += 1
                else:
                    overlapping_alignment_count += 1
            print(reference.name)
            print('  Length:       ', int_to_str(seq_len, max_v) + ' bp')
            if alignments:
                mean_depth = sum(depths) / seq_len
                mean_mismatch_rate = 100.0 * sum(mismatch_rates) / aligned_len
                mean_insertion_rate = 100.0 * sum(insertion_rates_multi_base) / aligned_len
                mean_deletion_rate = 100.0 * sum(deletion_rates) / aligned_len
                if VERBOSITY > 0:
                    print('  Alignments:         ', int_to_str(len(alignments), max_v))
                if VERBOSITY > 1:
                    print('    Contained:        ', int_to_str(contained_alignment_count, max_v))
                    print('    Overlapping:      ', int_to_str(overlapping_alignment_count, max_v))
                print('  Covered length:     ', int_to_str(aligned_len, max_v) + ' bp')
                print('  Covered fraction:   ', float_to_str(aligned_percent, 2, max_v) + '%')
                if VERBOSITY > 1:
                    print('  Mean read depth:    ', float_to_str(mean_depth, 2, max_v))
                    print('  Mismatch rate:      ',
                          float_to_str(mean_mismatch_rate, 2, max_v) +'%')
                    print('  Insertion rate:     ',
                          float_to_str(mean_insertion_rate, 2, max_v) + '%')
                    print('  Deletion rate:      ',
                          float_to_str(mean_deletion_rate, 2, max_v) + '%')
            else:
                print('  Alignments:         ', int_to_str(0, max_v))
            print()


if __name__ == '__main__':
    main()
