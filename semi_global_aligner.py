#!/usr/bin/env python3
'''
Bikepath - a sensitive semi-global long read aligner

This is a script to align error-prone long reads (e.g. PacBio or Nanopore) to one or more
references in a semi-global manner. Semi-global alignment allows for unpenalised end gaps, but the
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
'''

import subprocess
import sys
import os
import argparse
import time
from multiprocessing.dummy import Pool as ThreadPool
from multiprocessing import cpu_count
import threading

SCRIPT_DIR = os.path.dirname(os.path.realpath(__file__))
sys.path.append(os.path.join(SCRIPT_DIR, 'lib'))
from misc import int_to_str, float_to_str, check_file_exists, quit_with_error, check_graphmap, \
                 get_mean_and_st_dev, print_progress_line, reverse_complement, \
                 print_section_header
from cpp_function_wrappers import semi_global_alignment, new_kmer_positions, add_kmer_positions, \
                                  delete_all_kmer_positions, start_extension_alignment, \
                                  end_extension_alignment, \
                                  get_random_sequence_alignment_mean_and_std_dev
from read_ref import Read, Reference, load_references, load_long_reads
from alignment import Alignment, AlignmentScoringScheme

# Used to ensure that multiple threads writing to the same SAM file don't write at the same time.
sam_write_lock = threading.Lock()

'''
VERBOSITY controls how much the script prints to the screen.
0 = nothing is printed
1 = a relatively simple output is printed
2 = a more thorough output is printed, including details on each Seqan alignment
3 = even more output is printed, including stuff from the C++ code
4 = tons of stuff is printed, including all k-mer positions in each Seqan alignment
'''
VERBOSITY = 0

'''
EXPECTED_SLOPE is the anticipated reference to read ratio. It is used by the C++ Seqan code to
rotate the common k-mer rectangles when looking for alignment lines. It is a global because it will
be constantly updated as reads are aligned.
TOTAL_REF_LENGTH and TOTAL_READ_LENGTH are the totals used to calculate EXPECTED_SLOPE. They start
at 10000 (not the more literal value of 0) so EXPECTED_SLOPE isn't too prone to fluctuation at the
start.
'''
EXPECTED_SLOPE = 1.0
TOTAL_REF_LENGTH = 10000
TOTAL_READ_LENGTH = 10000


def main():
    '''
    If this script is run on its own, execution starts here.
    '''
    full_command = ' '.join(sys.argv)
    args = get_arguments()
    check_file_exists(args.ref)
    check_file_exists(args.reads)
    if not args.no_graphmap:
        check_graphmap(args.graphmap_path)

    references = load_references(args.ref, VERBOSITY)
    read_dict, read_names = load_long_reads(args.reads, VERBOSITY)
    scoring_scheme = AlignmentScoringScheme(args.scores)

    semi_global_align_long_reads(references, args.ref, read_dict, read_names, args.reads,
                                 args.temp_dir, args.graphmap_path, args.threads, scoring_scheme,
                                 args.low_score, not args.no_graphmap, args.keep_bad, args.kmer,
                                 args.min_len, args.sam, full_command, args.allowed_overlap,
                                 VERBOSITY)
    sys.exit(0)

def get_arguments():
    '''
    Specifies the command line arguments required by the script.
    '''
    parser = argparse.ArgumentParser(description='Semi-global long read aligner',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('--ref', type=str, required=True, default=argparse.SUPPRESS,
                        help='FASTA file containing one or more reference sequences')
    parser.add_argument('--reads', type=str, required=True, default=argparse.SUPPRESS,
                        help='FASTQ file of long reads')
    parser.add_argument('--sam', type=str, required=True, default=argparse.SUPPRESS,
                        help='SAM file of resulting alignments')
    
    add_aligning_arguments(parser, False)

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

def add_aligning_arguments(parser, hide_help):
    '''
    Adds the aligning-specific arguments to the parser.
    '''
    temp_dir_help = argparse.SUPPRESS if hide_help else 'Temp directory for working files ' + \
                                                        '("PID" will be replaced with the ' + \
                                                        'process ID)'
    parser.add_argument('--temp_dir', type=str, required=False, default='align_temp_PID',
                        help=temp_dir_help)

    no_graphmap_help = argparse.SUPPRESS if hide_help else 'Do not use GraphMap as a ' + \
                                                           'first-pass aligner (default: ' + \
                                                           'GraphMap is used)'
    parser.add_argument('--no_graphmap', action='store_true', default=argparse.SUPPRESS,
                        help=no_graphmap_help)

    graphmap_path_help = argparse.SUPPRESS if hide_help else 'Path to the GraphMap executable'
    parser.add_argument('--graphmap_path', type=str, required=False, default='graphmap',
                        help=graphmap_path_help)

    scores_help = argparse.SUPPRESS if hide_help else 'Comma-delimited string of alignment ' + \
                                                      'scores: match, mismatch, gap open, gap ' + \
                                                      'extend'
    parser.add_argument('--scores', type=str, required=False, default='3,-6,-5,-2',
                        help=scores_help)

    low_score_help = argparse.SUPPRESS if hide_help else 'Score threshold - alignments below ' + \
                                                         'this are considered poor (default: ' + \
                                                         'set threshold automatically)'
    parser.add_argument('--low_score', type=float, required=False, default=argparse.SUPPRESS,
                        help=low_score_help)

    min_len_help = argparse.SUPPRESS if hide_help else 'Minimum alignment length (bp) - ' + \
                                                       'exclude alignments shorter than this ' + \
                                                       'length'
    parser.add_argument('--min_len', type=float, required=False, default=100,
                        help=min_len_help)

    keep_bad_help = argparse.SUPPRESS if hide_help else 'Include alignments in the results ' + \
                                                        'even if they are below the low score ' + \
                                                        'threshold (default: low-scoring ' + \
                                                        'alignments are discarded)'
    parser.add_argument('--keep_bad', action='store_true', default=argparse.SUPPRESS,
                        help=keep_bad_help)

    allowed_overlap_help = argparse.SUPPRESS if hide_help else 'Allow this much overlap ' + \
                                                               'between alignments in a ' + \
                                                               'single read'
    parser.add_argument('--allowed_overlap', type=int, required=False, default=100,
                        help=allowed_overlap_help)

    kmer_help = argparse.SUPPRESS if hide_help else 'K-mer size used for seeding alignments'
    parser.add_argument('--kmer', type=int, required=False, default=7,
                        help=kmer_help)

def fix_up_arguments(args):
    '''
    Repairs issues with the arguments, like not existing. We don't use None/False as a default
    in add_argument because it makes the help text look weird.
    '''
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
                                 low_score_threshold, use_graphmap, keep_bad, kmer_size,
                                 min_align_length, sam_filename, full_command, allowed_overlap,
                                 verbosity=None):
    '''
    This function does the primary work of this module: aligning long reads to references in an
    end-gap-free, semi-global manner. It returns a list of Read objects which contain their
    alignments.
    If seqan_all is True, then every Alignment object will be refined by using Seqan.
    If seqan_all is False, then only the overlap alignments and a small set of long contained
    alignments will be run through Seqan.
    '''
    if verbosity:
        global VERBOSITY
        VERBOSITY = verbosity
        
    # If the user supplied a low score threshold, we use that. Otherwise, we'll use the median
    # score minus three times the MAD.
    print_section_header('Determining low-score threshold', VERBOSITY)
    if low_score_threshold:
        if VERBOSITY > 0:
            print('Using user-supplied low score threshold: ' +
                  float_to_str(low_score_threshold, 2) + '\n')
    else:
        if VERBOSITY > 0:
            print('Automatically choosing a low score threshold using random alignments.\n')
        std_devs_over_mean = 5
        low_score_threshold, rand_mean, rand_std_dev = get_auto_score_threshold(scoring_scheme,
                                                                                std_devs_over_mean)
        if VERBOSITY > 0:
            print('Random alignment mean score: ' + float_to_str(rand_mean, 2))
            print('         standard deviation: ' + float_to_str(rand_std_dev, 2, rand_mean))
            print()
            print('Low score threshold = ' + float_to_str(rand_mean, 2) + ' + ' + \
                  str(std_devs_over_mean) + ' x ' +  float_to_str(rand_std_dev, 2) + ' = ' + \
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
        sam_file.write('ID:' + 'ALIGNER_NAME')
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
                                                  scoring_scheme)
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
            print()

        # Use Seqan to extend the GraphMap alignments so they are fully semi-global. In this
        # process, some alignments will be discarded (those too far from being semi-global).
        semi_global_graphmap_alignments = extend_to_semi_global(graphmap_alignments, scoring_scheme)
        if VERBOSITY > 3 and semi_global_graphmap_alignments:
            print_section_header('GraphMap alignments after extension', VERBOSITY)
            for alignment in semi_global_graphmap_alignments:
                print(alignment)
                if VERBOSITY > 4:
                    print(alignment.cigar)
            print()

        # Gather some statistics about the alignments.
        percent_ids = [x.percent_identity for x in graphmap_alignments]
        scores = [x.scaled_score for x in graphmap_alignments]
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

    if VERBOSITY > 0:
        if VERBOSITY < 3:
            print_section_header('Aligning reads', VERBOSITY)
        num_realignments = len(reads_to_align)
        max_v = len(read_dict)
        if use_graphmap:
            print('Reads completed by GraphMap:', int_to_str(len(completed_reads), max_v))
            print('Reads to be realigned:      ', int_to_str(num_realignments, max_v))
            print()
    if VERBOSITY == 1:
        print_progress_line(0, num_realignments, prefix='Read: ')
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
                                     sam_filename, allowed_overlap)
            completed_count += 1
            if VERBOSITY == 1:
                print_progress_line(completed_count, num_realignments, prefix='Read: ')
            if VERBOSITY > 1:
                print(output, end='')

    # If multi-threaded, use a thread pool.
    else:
        pool = ThreadPool(threads)
        arg_list = []
        for read in reads_to_align:
            arg_list.append((read, reference_dict, scoring_scheme, kmer_positions_ptr,
                             low_score_threshold, keep_bad, kmer_size, min_align_length,
                             sam_filename, allowed_overlap))
        for output in pool.imap(seqan_alignment_one_arg, arg_list, 1):
            completed_count += 1
            if VERBOSITY == 1:
                print_progress_line(completed_count, num_realignments, prefix='Read: ')
            if VERBOSITY > 1:
                print(output, end='')

    # We're done with the C++ KmerPositions object, so delete it now.
    delete_all_kmer_positions(kmer_positions_ptr)
    
    if VERBOSITY == 1:
        print()

    # Output a summary of the reads' alignments.
    fully_aligned, partially_aligned, unaligned = group_reads_by_fraction_aligned(read_dict)
    ref_bases_aligned = 0
    for read in read_dict.values():
        ref_bases_aligned += read.get_reference_bases_aligned()
    print_section_header('Read alignment summary', VERBOSITY)
    if VERBOSITY > 0:
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

    return read_dict

def print_graphmap_summary_table(graphmap_alignments, percent_id_mean, percent_id_std_dev,
                                 score_mean, score_std_dev):
    '''
    Prints a small table showing some details about the GraphMap alignments.
    '''
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
    print()

def extend_to_semi_global(alignments, scoring_scheme):
    '''
    This function returns truly semi-global alignments made from the input alignments.
    '''
    if VERBOSITY > 3 and alignments:
        print_section_header('Extending alignments', VERBOSITY)

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
    '''
    This function runs GraphMap for the given inputs and produces a SAM file at the given location.
    '''
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
    # progress lines. This helps when piping the output to file (otherwise the output can be
    # excessively large).
    process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    line = ''
    last_progress = -1.0
    read_progress_started = False
    read_progress_finished = False
    while process.poll() is None:
        graphmap_output = process.stderr.read(1)
        if VERBOSITY > 0:
            line += graphmap_output
            if line.endswith('\n') or line.endswith('\r'):
                if line.strip():
                    if 'CPU time' in line:
                        read_progress_started = True
                        trimmed_line = line.strip().split('] ')[2].split(', l')[0]
                        progress = float(trimmed_line.split('(')[1].split(')')[0][:-1])
                        progress_rounded_down = float(int(progress))
                        if progress_rounded_down > last_progress:
                            print('\r' + trimmed_line, end='')
                            last_progress = progress_rounded_down
                    elif VERBOSITY > 1:
                        if read_progress_started and not read_progress_finished:
                            print()
                            read_progress_finished = True
                        print(line, end='')
                line = ''
    if VERBOSITY == 1:
        print('\n')
    if VERBOSITY > 1:
        print()

    # Clean up.
    if os.path.isfile(fasta + '.gmidx'):
        os.remove(fasta + '.gmidx')
    if os.path.isfile(fasta + '.gmidxsec'):
        os.remove(fasta + '.gmidxsec')

    if not os.path.isfile(sam_file):
        quit_with_error('GraphMap failure')

def get_graphmap_version(graphmap_path):
    '''
    Returns the version of GraphMap.
    '''
    command = [graphmap_path, '-h']
    process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out, err = process.communicate()
    allout = out + err
    if 'Version: v' not in allout:
        command = [graphmap_path, 'align', '-h']
        process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        out, err = process.communicate()
        allout = out + err
    if 'Version: v' not in allout:
        return 0.0
    version_i = allout.find('Version: v')
    version = allout[version_i + 10:]
    version = version.split()[0]
    version = '.'.join(version.split('.')[0:2])
    return float(version)

def load_sam_alignments(sam_filename, read_dict, reference_dict, scoring_scheme):
    '''
    This function returns a list of Alignment objects from the given SAM file.
    '''
    sam_alignments = []
    sam_file = open(sam_filename, 'r')
    for line in sam_file:
        line = line.strip()
        if line and not line.startswith('@') and line.split('\t', 3)[2] != '*':
            sam_alignments.append(Alignment(sam_line=line, read_dict=read_dict,
                                            reference_dict=reference_dict,
                                            scoring_scheme=scoring_scheme))
    sam_file.close()
    return sam_alignments

def seqan_alignment_one_arg(all_args):
    '''
    This is just a one-argument version of seqan_alignment to make it easier to use that function
    in a thread pool.
    '''
    read, reference_dict, scoring_scheme, kmer_positions_ptr, low_score_threshold, keep_bad, \
        kmer_size, min_align_length, sam_filename, allowed_overlap = all_args
    return seqan_alignment(read, reference_dict, scoring_scheme, kmer_positions_ptr,
                           low_score_threshold, keep_bad, kmer_size, min_align_length,
                           sam_filename, allowed_overlap)

def seqan_alignment(read, reference_dict, scoring_scheme, kmer_positions_ptr, low_score_threshold,
                    keep_bad, kmer_size, min_align_length, sam_filename, allowed_overlap):
    '''
    Aligns a single read against all reference sequences using Seqan.
    '''
    start_time = time.time()
    output = ''
    if VERBOSITY > 1:
        output += str(read) + '\n'
    if VERBOSITY > 2:
        output += '-' * len(str(read)) + '\n'

    starting_graphmap_alignments = len(read.alignments)
    if VERBOSITY > 2:
        output += 'Graphmap alignments:\n'
        if starting_graphmap_alignments:
            for alignment in read.alignments:
                output += '  ' + str(alignment) + '\n'
        else:
            output += '  None\n'

    # print(read, EXPECTED_SLOPE, flush=True) # TEMP

    results = semi_global_alignment(read.name, read.sequence, VERBOSITY,
                                    EXPECTED_SLOPE, kmer_positions_ptr,
                                    scoring_scheme.match, scoring_scheme.mismatch,
                                    scoring_scheme.gap_open, scoring_scheme.gap_extend,
                                    low_score_threshold, keep_bad, kmer_size).split(';')
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
                    output += '  ' + str(alignment) + '\n'
                    if VERBOSITY > 3:
                        output += alignment.cigar + '\n'

    read.remove_conflicting_alignments(allowed_overlap)
    if not keep_bad:
        read.remove_low_score_alignments(low_score_threshold)
    read.remove_short_alignments(min_align_length)

    if VERBOSITY > 2:
        output += 'Final alignments:\n'
    if VERBOSITY > 1:
        if read.alignments:
            for alignment in read.alignments:
                output += '  ' + str(alignment) + '\n'
        else:
            output += '  None\n'
        output += '\n'

    # Write alignments to SAM.
    if sam_filename and read.alignments:
        sam_write_lock.acquire()
        sam_file = open(sam_filename, 'a')
        for alignment in read.alignments:
            sam_file.write(alignment.get_sam_line())
        sam_file.close()
        sam_write_lock.release()

    update_expected_slope(read, low_score_threshold)
    return output

def group_reads_by_fraction_aligned(read_dict):
    '''
    Groups reads into three lists:
      1) Fully aligned
      2) Partially aligned
      3) Unaligned
    '''
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
    '''
    This function determines a good low score threshold for the alignments. To do this it examines
    the distribution of scores acquired by aligning random sequences.
    '''
    # TO DO: make the random alignments run in separate threads to be a bit faster
    mean, std_dev = get_random_sequence_alignment_mean_and_std_dev(100, 10000, scoring_scheme)
    threshold = mean + (std_devs_over_mean * std_dev)

    # Keep the threshold bounded to sane levels.
    threshold = min(threshold, 95.0)
    threshold = max(threshold, 50.0)
    return threshold, mean, std_dev

def update_expected_slope(read, low_score_threshold):
    '''
    This function updates the EXPECTED_SLOPE and ALIGNMENTS_CONTRIBUTING_TO_EXPECTED_SLOPE global
    variables using a read, but only if the read's alignment looks good.
    '''
    if len(read.alignments) == 1 and \
       read.alignments[0].read_start_pos == 0 and read.alignments[0].read_end_gap == 0 and \
       read.alignments[0].scaled_score > low_score_threshold:
        global EXPECTED_SLOPE
        global TOTAL_REF_LENGTH
        global TOTAL_READ_LENGTH
        TOTAL_REF_LENGTH += read.alignments[0].get_aligned_ref_length()
        TOTAL_READ_LENGTH += read.alignments[0].get_aligned_read_length()
        EXPECTED_SLOPE = TOTAL_REF_LENGTH / TOTAL_READ_LENGTH

if __name__ == '__main__':
    main()
