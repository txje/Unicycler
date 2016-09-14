#!/usr/bin/env python3
"""
Unicycler align - a sensitive semi-global long read aligner

This is a script to align error-prone long reads (e.g. PacBio or Nanopore) to one or more
references in a semi-global manner. Semi-global alignment does not penalise end gaps, but the
alignment will continue until one of the two sequences ends. This includes cases where the two
sequences overlap and cases where one sequence is contained within the other:

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

Output: SAM file of alignments

Author: Ryan Wick
email: rrwick@gmail.com
"""

import subprocess
import sys
import os
import argparse
import time
import random
import shutil
import math
from multiprocessing.dummy import Pool as ThreadPool
from multiprocessing import cpu_count
import threading
from .misc import int_to_str, float_to_str, check_file_exists, quit_with_error, check_graphmap, \
    get_mean_and_st_dev, print_progress_line, print_section_header, \
    weighted_average_list
from .cpp_function_wrappers import semi_global_alignment, new_kmer_positions, add_kmer_positions, \
    delete_all_kmer_positions, \
    get_random_sequence_alignment_mean_and_std_dev
from .read_ref import load_references, load_long_reads
from .alignment import Alignment, AlignmentScoringScheme
from . import settings

# Used to ensure that multiple threads writing to the same SAM file don't write at the same time.
SAM_WRITE_LOCK = threading.Lock()

# VERBOSITY controls how much the script prints to the screen.
# 0 = nothing is printed
# 1 = a relatively simple output is printed
# 2 = a more thorough output is printed, including details on each Seqan alignment
# 3 = even more output is printed, including stuff from the C++ code
# 4 = tons of stuff is printed, including all k-mer positions in each Seqan alignment
VERBOSITY = 0

# EXPECTED_SLOPE is the anticipated reference to read ratio. It is used by the C++ Seqan code to
# rotate the common k-mer rectangles when looking for alignment lines. It is a global because it
# will be constantly updated as reads are aligned.
# TOTAL_REF_LENGTH and TOTAL_READ_LENGTH are the totals used to calculate EXPECTED_SLOPE. They
# start at 10000 (not the more literal value of 0) so EXPECTED_SLOPE isn't too prone to fluctuation
# at the start.
EXPECTED_SLOPE = 1.0
TOTAL_REF_LENGTH = 10000
TOTAL_READ_LENGTH = 10000


def main():
    """
    If this script is run on its own, execution starts here.
    """
    # Fix the random seed so the program produces the same output every time it's run.
    random.seed(0)

    full_command = ' '.join(sys.argv)
    args = get_arguments()
    check_file_exists(args.ref)
    check_file_exists(args.reads)
    if not args.no_graphmap:
        check_graphmap(args.graphmap_path)

    references = load_references(args.ref, VERBOSITY)
    read_dict, read_names, read_filename = load_long_reads(args.reads, VERBOSITY)
    scoring_scheme = AlignmentScoringScheme(args.scores)

    semi_global_align_long_reads(references, args.ref, read_dict, read_names, read_filename,
                                 args.temp_dir, args.graphmap_path, args.threads, scoring_scheme,
                                 [args.low_score], not args.no_graphmap, args.keep_bad, args.kmer,
                                 args.min_len, args.sam, full_command, args.allowed_overlap,
                                 args.extra_sensitive, VERBOSITY)
    sys.exit(0)


def get_arguments():
    """
    Specifies the command line arguments required by the script.
    """
    terminal_width = shutil.get_terminal_size().columns
    os.environ['COLUMNS'] = str(terminal_width)

    parser = argparse.ArgumentParser(description='Unicycler align - a sensitive semi-global long '
                                                 'read aligner',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('--ref', type=str, required=True, default=argparse.SUPPRESS,
                        help='FASTA file containing one or more reference sequences')
    parser.add_argument('--reads', type=str, required=True, default=argparse.SUPPRESS,
                        help='FASTQ or FASTA file of long reads')
    parser.add_argument('--sam', type=str, required=True, default=argparse.SUPPRESS,
                        help='SAM file of resulting alignments')

    add_aligning_arguments(parser, True)

    parser.add_argument('--extra_sensitive', action='store_true',
                        help='Perform slow but very sensitive alignment')
    parser.add_argument('--threads', type=int, required=False, default=argparse.SUPPRESS,
                        help='Number of CPU threads used to align (default: the number of '
                             'available CPUs)')
    parser.add_argument('--verbosity', type=int, required=False, default=1,
                        help='Level of stdout information (0 to 4)')

    args = parser.parse_args()

    global VERBOSITY
    VERBOSITY = args.verbosity

    fix_up_arguments(args)

    return args


def add_aligning_arguments(parser, show_help):
    """
    Adds the aligning-specific arguments to the parser.
    """
    parser.add_argument('--temp_dir', type=str, required=False, default='align_temp_PID',
                        help='Temp directory for working files ("PID" will be replaced with the '
                             'process ID)'
                             if show_help else argparse.SUPPRESS)
    parser.add_argument('--no_graphmap', action='store_true', default=argparse.SUPPRESS,
                        help='Do not use GraphMap as a first-pass aligner (default: GraphMap is '
                             'used)'
                             if show_help else argparse.SUPPRESS)
    parser.add_argument('--graphmap_path', type=str, required=False, default='graphmap',
                        help='Path to the GraphMap executable'
                             if show_help else argparse.SUPPRESS)
    parser.add_argument('--scores', type=str, required=False, default='3,-6,-5,-2',
                        help='Comma-delimited string of alignment scores: match, mismatch, '
                             'gap open, gap extend'
                             if show_help else argparse.SUPPRESS)
    parser.add_argument('--low_score', type=float, required=False, default=argparse.SUPPRESS,
                        help='Score threshold - alignments below this are considered poor '
                             '(default: set threshold automatically)'
                             if show_help else argparse.SUPPRESS)
    parser.add_argument('--min_len', type=float, required=False, default=100,
                        help='Minimum alignment length (bp) - exclude alignments shorter than this '
                             'length'
                             if show_help else argparse.SUPPRESS)
    parser.add_argument('--keep_bad', action='store_true', default=argparse.SUPPRESS,
                        help='Include alignments in the results even if they are below the low '
                             'score threshold (default: low-scoring alignments are discarded)'
                             if show_help else argparse.SUPPRESS)
    parser.add_argument('--allowed_overlap', type=int, required=False, default=100,
                        help='Allow this much overlap between alignments in a single read'
                             if show_help else argparse.SUPPRESS)
    parser.add_argument('--kmer', type=int, required=False, default=7,
                        help='K-mer size used for seeding alignments'
                             if show_help else argparse.SUPPRESS)


def fix_up_arguments(args):
    """
    Repairs issues with the arguments, like not existing. We don't use None/False as a default
    in add_argument because it makes the help text look weird.
    """
    try:
        args.low_score
    except AttributeError:
        args.low_score = None
    try:
        args.no_graphmap
    except AttributeError:
        args.no_graphmap = False
    try:
        args.keep_bad
    except AttributeError:
        args.keep_bad = False
    try:
        args.threads
    except AttributeError:
        args.threads = cpu_count()
        if VERBOSITY > 2:
            print('\nThread count set to', args.threads)

    # Add the process ID to the default temp directory so multiple instances can run at once in the
    # same directory.
    args.temp_dir = args.temp_dir.replace('PID', str(os.getpid()))


def semi_global_align_long_reads(references, ref_fasta, read_dict, read_names, reads_fastq,
                                 temp_dir, graphmap_path, threads, scoring_scheme,
                                 low_score_threshold_list, use_graphmap, keep_bad, kmer_size,
                                 min_align_length, sam_filename, full_command, allowed_overlap,
                                 extra_sensitive, verbosity=None, stdout_header='Aligning reads',
                                 display_low_score=True):
    """
    This function does the primary work of this module: aligning long reads to references in an
    end-gap-free, semi-global manner. It returns a list of Read objects which contain their
    alignments.
    If seqan_all is True, then every Alignment object will be refined by using Seqan.
    If seqan_all is False, then only the overlap alignments and a small set of long contained
    alignments will be run through Seqan.
    The low score threshold is taken as a list so the function can alter it and the caller can
    get the altered value.
    """
    if verbosity:
        global VERBOSITY
        VERBOSITY = verbosity

    # If the user supplied a low score threshold, we use that. Otherwise, we'll use the median
    # score minus three times the MAD.
    if display_low_score:
        print_section_header('Determining low-score threshold', VERBOSITY)
    low_score_threshold = low_score_threshold_list[0]
    if low_score_threshold is not None:
        if VERBOSITY > 0 and display_low_score:
            print('Using user-supplied threshold: ' + float_to_str(low_score_threshold, 2))
    else:
        if VERBOSITY > 0 and display_low_score:
            print('Automatically choosing a threshold using random alignments.\n')
        std_devs_over_mean = settings.AUTO_SCORE_STDEV_ABOVE_RANDOM_ALIGNMENT_MEAN
        low_score_threshold, rand_mean, rand_std_dev = get_auto_score_threshold(scoring_scheme,
                                                                                std_devs_over_mean)
        low_score_threshold_list[0] = low_score_threshold
        if VERBOSITY > 0 and display_low_score:
            print('Random alignment mean score: ' + float_to_str(rand_mean, 2))
            print('         standard deviation: ' + float_to_str(rand_std_dev, 2, rand_mean))
            print()
            print('Low score threshold = ' + float_to_str(rand_mean, 2) + ' + ' +
                  str(std_devs_over_mean) + ' x ' + float_to_str(rand_std_dev, 2) + ' = ' +
                  float_to_str(low_score_threshold, 2))

    reference_dict = {x.name: x for x in references}

    # Create the SAM file.
    if sam_filename:
        sam_file = open(sam_filename, 'w')

        # Header line.
        sam_file.write('@HD' + '\t')
        sam_file.write('VN:1.5' + '\t')
        sam_file.write('SO:unknown' + '\n')

        # Reference lines.
        for ref in references:
            sam_file.write('@SQ' + '\t')
            sam_file.write('SN:' + ref.name + '\t')
            sam_file.write('LN:' + str(ref.get_length()) + '\n')

        # Program line.
        sam_file.write('@PG' + '\t')
        sam_file.write('ID:' + 'unicycler_align')
        if full_command:
            sam_file.write('\tCL:' + full_command + '\t')
        sam_file.write('SC:' + str(scoring_scheme) + '\n')
        sam_file.close()

    # GraphMap can be used as a first-pass aligner. This has the advantage of saving time (GraphMap
    # is probably faster than the Seqan alignment) and it gives a good initial expected slope.
    if use_graphmap:

        # Make the temp directory, if necessary.
        temp_dir_exist_at_start = os.path.exists(temp_dir)
        if not temp_dir_exist_at_start:
            os.makedirs(temp_dir)

        # Run GraphMap and load in the resulting SAM.
        graphmap_sam = os.path.join(temp_dir, 'graphmap_alignments.sam')
        run_graphmap(ref_fasta, reads_fastq, graphmap_sam, graphmap_path, threads, scoring_scheme)
        graphmap_alignments = load_sam_alignments(graphmap_sam, read_dict, reference_dict,
                                                  scoring_scheme, threads, VERBOSITY)
        # Clean up files and directories.
        if use_graphmap:
            os.remove(graphmap_sam)
            if not temp_dir_exist_at_start:
                os.rmdir(temp_dir)

        if VERBOSITY > 3 and graphmap_alignments:
            print_section_header('GraphMap alignments before extension', VERBOSITY)
            for alignment in graphmap_alignments:
                print(alignment)
                if VERBOSITY > 4:
                    print(alignment.cigar)

        # Use Seqan to extend the GraphMap alignments so they are fully semi-global. In this
        # process, some alignments will be discarded (those too far from being semi-global).
        semi_global_graphmap_alignments = extend_to_semi_global(graphmap_alignments, scoring_scheme)
        if VERBOSITY > 3 and semi_global_graphmap_alignments:
            print_section_header('GraphMap alignments after extension', VERBOSITY)
            for alignment in semi_global_graphmap_alignments:
                print(alignment)
                if VERBOSITY > 4:
                    print(alignment.cigar)

        # Gather some statistics about the alignments.
        percent_ids = [x.percent_identity for x in semi_global_graphmap_alignments]
        scores = [x.scaled_score for x in semi_global_graphmap_alignments]
        percent_id_mean, percent_id_std_dev = get_mean_and_st_dev(percent_ids)
        score_mean, score_std_dev = get_mean_and_st_dev(scores)

        # Give the alignments to their corresponding reads.
        for alignment in semi_global_graphmap_alignments:
            read_dict[alignment.read.name].alignments.append(alignment)

        # We can now sort our reads into two different categories for further action:
        #   1) Reads with a single, high quality alignment in the middle of a reference. These reads
        #      are done!
        #   2) Reads that are incompletely (or not at all) aligned, have overlapping alignments or
        #      low quality alignments. These reads will be aligned again using Seqan.
        completed_reads = []
        reads_to_align = []
        for read_name in read_names:
            read = read_dict[read_name]
            update_expected_slope(read, low_score_threshold)
            if read.needs_seqan_realignment(low_score_threshold):
                reads_to_align.append(read)
            else:
                completed_reads.append(read)

        if VERBOSITY > 0:
            print_graphmap_summary_table(graphmap_alignments,
                                         percent_id_mean, percent_id_std_dev,
                                         score_mean, score_std_dev)
            max_v = len(read_dict)
            print()
            print('Reads completed by GraphMap:', int_to_str(len(completed_reads), max_v))
            print('Reads to be realigned:      ', int_to_str(len(reads_to_align), max_v))

        # Write GraphMap alignments to SAM.
        if sam_filename:
            sam_file = open(sam_filename, 'a')
            for read in completed_reads:
                for alignment in read.alignments:
                    sam_file.write(alignment.get_sam_line())
            sam_file.close()

            # OPTIONAL TO DO: for reads which are completed, I could still try to refine GraphMap
            #                 alignments using my Seqan code. I'm not sure if it's worth it, so I
            #                 should give it a try to see what kind of difference it makes.

    # If we aren't using GraphMap as a first pass, then we align every read using Seqan.
    else:
        reads_to_align = [read_dict[x] for x in read_names]

    num_alignments = len(reads_to_align)
    if VERBOSITY > 0:
        if VERBOSITY < 3:
            print_section_header(stdout_header, VERBOSITY, last_newline=(VERBOSITY == 1))
    if VERBOSITY == 1:
        print_progress_line(0, num_alignments, prefix='Read: ')
    completed_count = 0

    # Create a C++ KmerPositions object and add each reference sequence.
    kmer_positions_ptr = new_kmer_positions()
    for ref in references:
        add_kmer_positions(kmer_positions_ptr, ref.name, ref.sequence, kmer_size)

    # If single-threaded, just do the work in a simple loop.
    if threads == 1:
        for read in reads_to_align:
            output = seqan_alignment(read, reference_dict, scoring_scheme, kmer_positions_ptr,
                                     low_score_threshold, keep_bad, kmer_size, min_align_length,
                                     sam_filename, allowed_overlap, use_graphmap, extra_sensitive)
            completed_count += 1
            if VERBOSITY == 1:
                print_progress_line(completed_count, num_alignments, prefix='Read: ')
            if VERBOSITY > 1:
                fraction = str(completed_count) + '/' + str(num_alignments) + ': '
                print('\n' + fraction + output[1:], end='', flush=True)

    # If multi-threaded, use a thread pool.
    else:
        pool = ThreadPool(threads)
        arg_list = []
        for read in reads_to_align:
            arg_list.append((read, reference_dict, scoring_scheme, kmer_positions_ptr,
                             low_score_threshold, keep_bad, kmer_size, min_align_length,
                             sam_filename, allowed_overlap, use_graphmap, extra_sensitive))

        # If the verbosity is 1, then the order doesn't matter, so use imap_unordered to deliver
        # the results evenly. If the verbosity is higher, deliver the results in order with imap.
        if VERBOSITY > 1:
            imap_function = pool.imap
        else:
            imap_function = pool.imap_unordered

        for output in imap_function(seqan_alignment_one_arg, arg_list):
            completed_count += 1
            if VERBOSITY == 1:
                print_progress_line(completed_count, num_alignments, prefix='Read: ')
            if VERBOSITY > 1:
                fraction = str(completed_count) + '/' + str(num_alignments) + ': '
                print('\n' + fraction + output[1:], end='', flush=True)

    # We're done with the C++ KmerPositions object, so delete it now.
    delete_all_kmer_positions(kmer_positions_ptr)

    if VERBOSITY == 1:
        print()

    print_alignment_summary_table(read_dict, VERBOSITY)
    return read_dict


def print_alignment_summary_table(read_dict, verbosity):
    """
    Outputs a summary of the reads' alignments, grouping them by fully aligned, partially aligned
    and unaligned.
    """
    if verbosity == 0:
        return

    fully_aligned, partially_aligned, unaligned = group_reads_by_fraction_aligned(read_dict)
    ref_bases_aligned = 0
    for read in read_dict.values():
        ref_bases_aligned += read.get_reference_bases_aligned()
    print_section_header('Read alignment summary', verbosity)
    max_v = max(len(read_dict), ref_bases_aligned)
    print('Total read count:       ', int_to_str(len(read_dict), max_v))
    print('Fully aligned reads:    ', int_to_str(len(fully_aligned), max_v))
    print('Partially aligned reads:', int_to_str(len(partially_aligned), max_v))
    if VERBOSITY > 2 and partially_aligned:
        print('    ' + ', '.join([x.name for x in partially_aligned]))
    print('Unaligned reads:        ', int_to_str(len(unaligned), max_v))
    if VERBOSITY > 2 and unaligned:
        print('    ' + ', '.join([x.name for x in unaligned]))
    print('Total bases aligned:    ', int_to_str(ref_bases_aligned, max_v) + ' bp')

    identities = []
    lengths = []
    for read in fully_aligned + partially_aligned:
        identities += [x.percent_identity for x in read.alignments]
        lengths += [x.get_aligned_ref_length() for x in read.alignments]
    mean_identity = weighted_average_list(identities, lengths)
    print('Mean alignment identity:', float_to_str(mean_identity, 1, max_v) + '%')


def print_graphmap_summary_table(graphmap_alignments, percent_id_mean, percent_id_std_dev,
                                 score_mean, score_std_dev):
    """
    Prints a small table showing some details about the GraphMap alignments.
    """
    print_section_header('Graphmap alignment summary', VERBOSITY)
    print('Total alignments:', int_to_str(len(graphmap_alignments)))
    print()

    table_lines = ['',
                   'Identity:',
                   'Score:']

    pad_length = max([len(x) for x in table_lines]) + 2
    table_lines = [x.ljust(pad_length) for x in table_lines]

    table_lines[0] += 'Mean'
    table_lines[1] += float_to_str(percent_id_mean, 2)
    if percent_id_mean:
        table_lines[1] += '%'
    table_lines[2] += float_to_str(score_mean, 2)

    pad_length = max([len(x) for x in table_lines]) + 2
    table_lines = [x.ljust(pad_length) for x in table_lines]

    table_lines[0] += 'Std dev'
    table_lines[1] += float_to_str(percent_id_std_dev, 2)
    if percent_id_std_dev:
        table_lines[1] += '%'
    table_lines[2] += float_to_str(score_std_dev, 2)

    for line in table_lines:
        print(line)
    print()
    print('Current mean reference length / read length:', float_to_str(EXPECTED_SLOPE, 5))


def extend_to_semi_global(alignments, scoring_scheme):
    """
    This function returns truly semi-global alignments made from the input alignments.
    """
    if VERBOSITY > 3 and alignments:
        print_section_header('Extending alignments', VERBOSITY, last_newline=False)

    semi_global_alignments = []
    for alignment in alignments:
        total_missing_bases = alignment.get_total_missing_bases()

        # If an input alignment is already semi-global, then it's included in the output.
        if total_missing_bases == 0:
            semi_global_alignments.append(alignment)

        # If an input alignment is almost semi-global (below a threshold), and not too close to the
        # end of the reference, then it is extended to make it semi-global.
        elif total_missing_bases <= settings.ALLOWED_MISSING_GRAPHMAP_BASES:
            missing_start = alignment.get_missing_bases_at_start()
            missing_end = alignment.get_missing_bases_at_end()
            if missing_start and alignment.ref_start_pos >= 2 * missing_start:
                alignment.extend_start(scoring_scheme, VERBOSITY)
            if missing_end and alignment.ref_end_gap >= 2 * missing_end:
                alignment.extend_end(scoring_scheme, VERBOSITY)
            semi_global_alignments.append(alignment)

        # If an input alignment is above the threshold (not close to being semi-global), it is
        # discarded.
        else:
            pass

    return semi_global_alignments


def run_graphmap(fasta, long_reads_fastq, sam_file, graphmap_path, threads, scoring_scheme):
    """
    This function runs GraphMap for the given inputs and produces a SAM file at the given location.
    """
    graphmap_version = get_graphmap_version(graphmap_path)

    # Build the GraphMap command. There is a bit of difference if we're using a version before or
    # after v0.3.0.
    command = [graphmap_path]
    if graphmap_version >= 0.3:
        command.append('align')
    command += ['-r', fasta,
                '-d', long_reads_fastq,
                '-o', sam_file,
                '-t', str(threads),
                '-a', 'anchorgotoh']
    command += scoring_scheme.get_graphmap_parameters()

    print_section_header('Aligning with GraphMap', VERBOSITY)
    if VERBOSITY > 0:
        print(' '.join(command))

    # Print the GraphMap output as it comes. I gather up and display lines so I can display fewer
    # progress lines (only at every 0.1% of progress, instead of for every read). This helps when
    # piping the output to file (otherwise the output can be excessively large).
    process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    line = ''
    last_progress = -1.0
    step = settings.GRAPHMAP_PROGRESS_STEP
    read_progress_started = False
    read_progress_finished = False
    while process.poll() is None:
        graphmap_output = process.stderr.read(1).decode()
        if VERBOSITY > 0:
            line += graphmap_output
            if line.endswith('\n') or line.endswith('\r'):
                if line.strip():
                    if 'CPU time' in line:
                        read_progress_started = True
                        trimmed_line = line.strip().split('] ')[2].split(', l')[0]
                        progress = float(trimmed_line.split('(')[1].split(')')[0][:-1])
                        progress_rounded_down = math.floor(progress / step) * step
                        if progress == 100.0 or progress_rounded_down > last_progress:
                            print('\r' + trimmed_line, end='')
                            last_progress = progress_rounded_down
                    elif VERBOSITY > 1:
                        if read_progress_started and not read_progress_finished:
                            print()
                            read_progress_finished = True
                        print(line, end='')
                line = ''
    if VERBOSITY == 1:
        print('')

    # Clean up.
    if os.path.isfile(fasta + '.gmidx'):
        os.remove(fasta + '.gmidx')
    if os.path.isfile(fasta + '.gmidxsec'):
        os.remove(fasta + '.gmidxsec')

    if not os.path.isfile(sam_file):
        quit_with_error('GraphMap failure')


def get_graphmap_version(graphmap_path):
    """
    Returns the version of GraphMap.
    """
    command = [graphmap_path, '-h']
    process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out, err = process.communicate()
    all_out = (out + err).decode()
    if 'Version: v' not in all_out:
        command = [graphmap_path, 'align', '-h']
        process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        out, err = process.communicate()
        all_out = (out + err).decode()
    if 'Version: v' not in all_out:
        return 0.0
    version_i = all_out.find('Version: v')
    version = all_out[version_i + 10:]
    version = version.split()[0]
    version = '.'.join(version.split('.')[0:2])
    return float(version)


def load_sam_alignments(sam_filename, read_dict, reference_dict, scoring_scheme, threads,
                        verbosity):
    """
    This function returns a list of Alignment objects from the given SAM file.
    """
    print_section_header('Loading alignments', verbosity)

    # Load the SAM lines into a list.
    sam_lines = []
    sam_file = open(sam_filename, 'rt')
    for line in sam_file:
        line = line.strip()
        if line and not line.startswith('@') and line.split('\t', 3)[2] != '*':
            sam_lines.append(line)
    num_alignments = sum(1 for line in open(sam_filename) if not line.startswith('@'))
    if not num_alignments:
        return []
    if verbosity > 0:
        print_progress_line(0, num_alignments)

    if not sam_lines:
        if verbosity > 0:
            print('No alignments to load')
        return []

    # If single-threaded, just do the work in a simple loop.
    threads = 1  # TEMP
    sam_alignments = []
    last_progress = 0.0
    step = settings.LOADING_ALIGNMENTS_PROGRESS_STEP
    if threads == 1:
        for line in sam_lines:
            sam_alignments.append(Alignment(sam_line=line, read_dict=read_dict,
                                            reference_dict=reference_dict,
                                            scoring_scheme=scoring_scheme))
            if verbosity > 0:
                progress = 100.0 * len(sam_alignments) / num_alignments
                progress_rounded_down = math.floor(progress / step) * step
                if progress == 100.0 or progress_rounded_down > last_progress:
                    print_progress_line(len(sam_alignments), num_alignments)
                    last_progress = progress_rounded_down

    # # If multi-threaded, use processes.
    # else:
    #     sam_line_groups = chunkify(sam_lines, threads)
    #     manager = Manager()
    #     workers = []
    #     sam_alignments = manager.list([])
    #     for sam_line_group in sam_line_groups:
    #         child = Process(target=make_alignments, args=(sam_line_group, read_dict,
    #                                                       reference_dict, scoring_scheme,
    #                                                       sam_alignments))
    #         child.start()
    #         workers.append(child)
    #     while any(i.is_alive() for i in workers):
    #         time.sleep(0.1)
    #         if verbosity > 0:
    #             print_progress_line(len(sam_alignments), num_alignments)
    #     for worker in workers:
    #         worker.join()
    #     sam_alignments = sam_alignments._getvalue()

    # At this point, we should have loaded num_alignments alignments. But check to make sure and
    # fix up the progress line if any didn't load.
    if verbosity > 0:
        if len(sam_alignments) < num_alignments:
            print_progress_line(len(sam_alignments), len(sam_alignments))
        print()

    return sam_alignments


# def chunkify(full_list, pieces):
#     '''
#     http://stackoverflow.com/questions/2130016/
#     splitting-a-list-of-arbitrary-size-into-only-roughly-n-equal-parts
#     '''
#     return [full_list[i::pieces] for i in range(pieces)]

# def make_alignments(sam_lines, read_dict, reference_dict, scoring_scheme, alignments):
#     '''
#     Produces alignments from SAM lines and deposits them in a managed list.
#     '''
#     for line in sam_lines:
#         alignments.append(Alignment(line, read_dict, reference_dict, scoring_scheme))

def seqan_alignment_one_arg(all_args):
    """
    This is just a one-argument version of seqan_alignment to make it easier to use that function
    in a thread pool.
    """
    read, reference_dict, scoring_scheme, kmer_positions_ptr, low_score_threshold, keep_bad, \
        kmer_size, min_align_length, sam_filename, allowed_overlap, used_graphmap, \
        extra_sensitive = all_args
    return seqan_alignment(read, reference_dict, scoring_scheme, kmer_positions_ptr,
                           low_score_threshold, keep_bad, kmer_size, min_align_length,
                           sam_filename, allowed_overlap, used_graphmap, extra_sensitive)


def seqan_alignment(read, reference_dict, scoring_scheme, kmer_positions_ptr, low_score_threshold,
                    keep_bad, kmer_size, min_align_length, sam_filename, allowed_overlap,
                    used_graphmap, extra_sensitive):
    """
    Aligns a single read against all reference sequences using Seqan.
    """
    start_time = time.time()
    output = '\n'
    if VERBOSITY > 2:
        output += '\n'
    if VERBOSITY > 1:
        output += str(read) + '\n'
    if VERBOSITY > 2:
        output += '-' * len(str(read)) + '\n'

    if VERBOSITY > 2 and used_graphmap:
        starting_graphmap_alignments = len(read.alignments)
        output += 'Graphmap alignments:\n'
        if starting_graphmap_alignments:
            for alignment in read.alignments:
                output += '  ' + str(alignment) + '\n'
        else:
            output += '  None\n'

    # Don't bother trying to align reads too short to have a good alignment.
    if read.get_length() < min_align_length:
        if VERBOSITY > 1:
            output += '  too short to align\n'
        return output

    # print(read, EXPECTED_SLOPE, flush=True) # TEMP

    results = semi_global_alignment(read.name, read.sequence, VERBOSITY,
                                    EXPECTED_SLOPE, kmer_positions_ptr,
                                    scoring_scheme.match, scoring_scheme.mismatch,
                                    scoring_scheme.gap_open, scoring_scheme.gap_extend,
                                    low_score_threshold, keep_bad, kmer_size,
                                    extra_sensitive).split(';')
    alignment_strings = results[:-1]
    output += results[-1]

    for alignment_string in alignment_strings:
        alignment = Alignment(seqan_output=alignment_string, read=read,
                              reference_dict=reference_dict,
                              scoring_scheme=scoring_scheme)
        read.alignments.append(alignment)

    if VERBOSITY > 2:
        if not alignment_strings:
            output += '  None\n'
        else:
            output += 'All Seqan alignments (time to align = ' + \
                      float_to_str(time.time() - start_time, 3) + ' s):\n'
            for alignment in read.alignments:
                if alignment.alignment_type != 'SAM':
                    output += '  ' + alignment.get_str_no_read_name() + '\n'
                    if VERBOSITY > 3:
                        output += ''.join(alignment.cigar_parts) + '\n'

    read.remove_conflicting_alignments(allowed_overlap)
    if not keep_bad:
        read.remove_low_score_alignments(low_score_threshold)
    read.remove_short_alignments(min_align_length)

    if VERBOSITY > 2:
        output += 'Final alignments:\n'
    if VERBOSITY > 1:
        if read.alignments:
            for alignment in read.alignments:
                output += '  ' + alignment.get_str_no_read_name() + '\n'
        else:
            output += '  None\n'

    # Write alignments to SAM.
    if sam_filename and read.alignments:
        SAM_WRITE_LOCK.acquire()
        sam_file = open(sam_filename, 'a')
        for alignment in read.alignments:
            sam_file.write(alignment.get_sam_line())
        sam_file.close()
        SAM_WRITE_LOCK.release()

    update_expected_slope(read, low_score_threshold)
    return output


def group_reads_by_fraction_aligned(read_dict):
    """
    Groups reads into three lists:
      1) Fully aligned
      2) Partially aligned
      3) Unaligned
    """
    fully_aligned_reads = []
    partially_aligned_reads = []
    unaligned_reads = []
    for read in read_dict.values():
        fraction_aligned = read.get_fraction_aligned()
        if fraction_aligned == 1.0:
            fully_aligned_reads.append(read)
        elif fraction_aligned == 0.0:
            unaligned_reads.append(read)
        else:
            partially_aligned_reads.append(read)
    return fully_aligned_reads, partially_aligned_reads, unaligned_reads


def get_auto_score_threshold(scoring_scheme, std_devs_over_mean):
    """
    This function determines a good low score threshold for the alignments. To do this it examines
    the distribution of scores acquired by aligning random sequences.
    """
    # If the scoring scheme is a typical one, don't actually do the random alignments now - just
    # use precomputed values (made with a lot of iterations so they should be pretty good).
    scoring_scheme_str = str(scoring_scheme)
    if scoring_scheme_str == '1,0,0,0':
        mean, std_dev = 50.225667, 2.467919
    elif scoring_scheme_str == '0,-1,-1,-1':
        mean, std_dev = 49.024927, 2.724548
    elif scoring_scheme_str == '1,-1,-1,-1':
        mean, std_dev = 51.741783, 2.183467
    elif scoring_scheme_str == '5,-4,-8,-6':   # GraphMap
        mean, std_dev = 42.707636, 2.435548
    elif scoring_scheme_str == '5,-6,-10,0':   # BLASR
        mean, std_dev = 58.65047, 0.853201
    elif scoring_scheme_str == '2,-5,-2,-1':   # BWA-MEM
        mean, std_dev = 72.712148, 0.95266
    elif scoring_scheme_str == '1,-3,-5,-2':   # CUSHAW2 / blastn-short
        mean, std_dev = 46.257408, 2.162765
    elif scoring_scheme_str == '5,-11,-2,-4':  # proovread
        mean, std_dev = 73.221967, 1.363692
    elif scoring_scheme_str == '3,-6,-5,-2':   # Unicycler-align
        mean, std_dev = 61.656918, 1.314624
    elif scoring_scheme_str == '2,-3,-5,-2':   # blastn / dc-megablast
        mean, std_dev = 47.453862, 1.985947
    elif scoring_scheme_str == '1,-2,0,0':     # megablast
        mean, std_dev = 81.720641, 0.77204
    elif scoring_scheme_str == '0,-6,-5,-3':   # Bowtie2 end-to-end
        mean, std_dev = 62.647055, 1.738603
    elif scoring_scheme_str == '2,-6,-5,-3':   # Bowtie2 local
        mean, std_dev = 59.713806, 1.641191
    elif scoring_scheme_str == '1,-4,-6,-1':   # BWA
        mean, std_dev = 60.328393, 1.176776

    # If scheme doesn't match any of the above, then we have to actually do the random alignments.
    else:
        mean, std_dev = get_random_sequence_alignment_mean_and_std_dev(100, 25000, scoring_scheme)

    threshold = mean + (std_devs_over_mean * std_dev)

    # Keep the threshold bounded to sane levels.
    threshold = min(threshold, 95.0)
    threshold = max(threshold, 50.0)

    return threshold, mean, std_dev


def update_expected_slope(read, low_score_threshold):
    """
    This function updates the EXPECTED_SLOPE and ALIGNMENTS_CONTRIBUTING_TO_EXPECTED_SLOPE global
    variables using a read, but only if the read's alignment looks good.
    """
    if len(read.alignments) == 1 and read.alignments[0].read_start_pos == 0 and \
            read.alignments[0].read_end_gap == 0 and \
            read.alignments[0].scaled_score > low_score_threshold:
        global EXPECTED_SLOPE
        global TOTAL_REF_LENGTH
        global TOTAL_READ_LENGTH
        TOTAL_REF_LENGTH += read.alignments[0].get_aligned_ref_length()
        TOTAL_READ_LENGTH += read.alignments[0].get_aligned_read_length()
        EXPECTED_SLOPE = TOTAL_REF_LENGTH / TOTAL_READ_LENGTH
