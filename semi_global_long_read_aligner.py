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
import ctypes

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
    long_reads = semi_global_align_long_reads(args.ref, args.reads, args.sam_raw,
                                              args.sam, args.temp_dir, args.path,
                                              True, args.threads)
    write_reference_errors_to_table(args.ref, long_reads, args.table, True)
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
    parser.add_argument('--reads', type=str, required=True, default=argparse.SUPPRESS,
                        help='FASTQ file of long reads')
    parser.add_argument('--sam', type=str, required=True, default=argparse.SUPPRESS,
                        help='SAM file of alignments after QC filtering')
    parser.add_argument('--sam_raw', type=str, required=False,
                        help='SAM file of unfiltered alignments')
    parser.add_argument('--table', type=str, required=False,
                        help='Path and/or prefix for table files summarising reference errors')
    parser.add_argument('--temp_dir', type=str, required=False, default='align_temp',
                        help='Temp directory for working files')
    parser.add_argument('--path', type=str, required=False, default='graphmap',
                        help='Path to the GraphMap executable')
    parser.add_argument('--threads', type=int, required=False, default=8,
                        help='Number of threads used by GraphMap')
    return parser.parse_args()

def semi_global_align_long_reads(ref_fasta, long_reads_fastq, sam_raw, sam_filtered, temp_dir,
                                 graphmap_path, print_summary, threads):
    '''
    This function does the primary work of this module: aligning long reads to references in an
    end-gap-free, semi-global manner. It returns a list of LongRead objects which contain their
    alignments.
    '''
    long_reads = load_long_reads(long_reads_fastq)

    if not sam_raw:
        temp_sam_raw = True
        sam_raw = os.path.join(temp_dir, 'alignments.sam')
    else:
        temp_sam_raw = False

    seq_by_seq_graphmap_alignment(ref_fasta, long_reads_fastq, sam_raw, graphmap_path, temp_dir,
                                  print_summary, threads)
    alignments = load_alignments(sam_raw, ref_fasta)
    if print_summary:
        max_v = max(100, len(alignments))
        print()
        print('Alignment summary')
        print('-----------------')
        print('Total raw GraphMap alignments:      ', int_to_str(len(alignments), max_v))

    # Give the alignments to their corresponding reads.
    for alignment in alignments:
        long_reads[alignment.read_name].alignments.append(alignment)

    # Filter the alignments based on conflicting read position.
    filtered_alignments = []
    for read in long_reads.itervalues():
        read.remove_conflicting_alignments()
        filtered_alignments += read.alignments
    if print_summary:
        print('Alignments after conflict filtering:', int_to_str(len(filtered_alignments), max_v))

    # Filter the alignments based on identity.
    mean_id, std_dev_id = get_mean_and_st_dev_identity(filtered_alignments, True)
    if mean_id == 0.0 or std_dev_id == -1:
        low_id_cutoff = 75.0
        if print_summary:
            print('Not enough alignments to automatically set a low identity cutoff. Using 75%.')
    else:
        low_id_cutoff = mean_id - (3.0 * std_dev_id)
        if print_summary:
            print('Complete alignment identity mean:   ', float_to_str(mean_id, max_v) + '%')
            print('              standard deviation:   ', float_to_str(std_dev_id, max_v) + '%')
            print('Low identity cutoff:                ', float_to_str(low_id_cutoff, max_v) + '%')
    filtered_alignments = []
    for read in long_reads.itervalues():
        read.remove_low_id_alignments(low_id_cutoff)
        filtered_alignments += read.alignments

    if print_summary:
        print('Alignments after identity filtering:', int_to_str(len(filtered_alignments), max_v))
        print()
        print('Read summary')
        print('------------')
        max_v = len(long_reads)
        print('Total read count:       ', int_to_str(len(long_reads), max_v))
        fully_aligned_count = 0
        partially_aligned_count = 0
        unaligned_count = 0
        for read in long_reads.itervalues():
            fraction_aligned = read.get_fraction_aligned()
            if fraction_aligned == 1.0:
                fully_aligned_count += 1
            elif fraction_aligned == 0.0:
                unaligned_count += 1
            else:
                partially_aligned_count += 1
        print('Fully aligned reads:    ', int_to_str(fully_aligned_count, max_v))
        print('Partially aligned reads:', int_to_str(partially_aligned_count, max_v))
        print('Unaligned reads:        ', int_to_str(unaligned_count, max_v))
        print()

    # FUTURE POSSIBILITY: FOR ANY READS WHICH ARE LACKING MAPPED REGIONS, TRY AGAIN WITH A MORE
    # SENSITIVE SEARCH.

    write_sam_file(filtered_alignments, sam_filtered)
    if temp_sam_raw:
        os.remove(sam_raw)

    return long_reads

def write_reference_errors_to_table(ref_fasta, long_reads, table_prefix, print_summary):
    '''
    Writes a table file summarising the alignment errors in terms of reference sequence position.
    Works in a brute force manner - could be made more efficient later if necessary.
    '''
    # If we are not making table files or printing a summary, then quit now because there's nothing
    # else to do.
    if not table_prefix and not print_summary:
        return

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

    if print_summary:
        max_v = max(100, sum([len(x.alignments) for x in long_reads.itervalues()]))
        print('Alignment summaries per reference')
        print('---------------------------------')

    ref_headers_and_seqs = load_fasta(ref_fasta)
    for header, seq in ref_headers_and_seqs:
        nice_header = get_nice_header(header)
        if table_prefix:
            header_for_filename = clean_str_for_filename(nice_header)
            if table_prefix.endswith('/'):
                table_filename = table_prefix + header_for_filename + '.txt'
            else:
                table_filename = table_prefix + '_' + header_for_filename + '.txt'
            table = open(table_filename, 'w')
            table.write('\t'.join(['base', 'read depth', 'mismatches', 'insertions',
                                   'deletions']) + '\n')
        seq_len = len(seq)
        depths = [0] * seq_len
        mismatches = [0] * seq_len
        insertions = [0] * seq_len
        deletions = [0] * seq_len
        alignments = []
        for read in long_reads.itervalues():
            for alignment in read.alignments:
                if alignment.ref_name == nice_header:
                    alignments.append(alignment)
                    for pos in xrange(alignment.ref_start_pos, alignment.ref_end_pos):
                        depths[pos] += 1
                    for pos in alignment.ref_mismatch_positions:
                        mismatches[pos] += 1
                    for pos in alignment.ref_insertion_positions:
                        insertions[pos] += 1
                    for pos in alignment.ref_deletion_positions:
                        deletions[pos] += 1
        if table_prefix:
            for i in xrange(seq_len):
                table.write('\t'.join([str(i+1), str(depths[i]), str(mismatches[i]),
                                       str(insertions[i]), str(deletions[i])]) + '\n')
        if print_summary:
            mismatch_rates = []
            insertion_rates = []
            deletion_rates = []
            for i in xrange(seq_len):
                depth = depths[i]
                if depth > 0.0:
                    mismatch_rates.append(mismatches[i] / depth)
                    insertion_rates.append(insertions[i] / depth)
                    deletion_rates.append(deletions[i] / depth)
            mean_depth = sum(depths) / seq_len
            mean_mismatch_rate = 100.0 * sum(mismatch_rates) / seq_len
            mean_insertion_rate = 100.0 * sum(insertion_rates) / seq_len
            mean_deletion_rate = 100.0 * sum(deletion_rates) / seq_len
            mean_id, std_dev_id = get_mean_and_st_dev_identity(alignments, False)
            contained_alignment_count = 0
            overlapping_alignment_count = 0
            for alignment in alignments:
                if alignment.is_whole_read():
                    contained_alignment_count += 1
                else:
                    overlapping_alignment_count += 1
            print(get_nice_header_and_len(header, seq))
            if alignments:
                print('  Total alignments:      ', int_to_str(len(alignments), max_v))
                print('  Contained alignments:  ', int_to_str(contained_alignment_count, max_v))
                print('  Overlapping alignments:', int_to_str(overlapping_alignment_count, max_v))
                print('  Mean read depth:       ', float_to_str(mean_depth, max_v))
                print('  Mismatch rate:         ', float_to_str(mean_mismatch_rate, max_v) + '%')
                print('  Insertion rate:        ', float_to_str(mean_insertion_rate, max_v) + '%')
                print('  Deletion rate:         ', float_to_str(mean_deletion_rate, max_v) + '%')
                print('  Mean identity:         ', float_to_str(mean_id, max_v) + '%')
                if std_dev_id != -1:
                    print('  Identity std dev:      ', float_to_str(std_dev_id, max_v) + '%')
            else:
                print('  Filtered alignments:   ', int_to_str(0, max_v))
            print()
        if table_prefix:
            table.close()

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

def seq_by_seq_graphmap_alignment(ref_fasta, long_reads_fastq, sam_file, graphmap_path,
                                  working_dir, print_summary, threads):
    '''
    This function runs GraphMap separately for each individual sequence in the reference.
    Resulting alignments are collected in a single SAM file.
    '''
    if print_summary:
        print()
        print('Raw GraphMap alignments per reference')
        print('-------------------------------------')

    final_sam = open(sam_file, 'w')
    ref_headers_and_seqs = load_fasta(ref_fasta)
    if print_summary and ref_headers_and_seqs:
        longest_header = max([len(get_nice_header_and_len(header, seq)) for \
                                   header, seq in ref_headers_and_seqs])
        read_count = line_count(long_reads_fastq) / 4
    for header, seq in ref_headers_and_seqs:
        nice_header = get_nice_header(header)
        one_seq_fasta = os.path.join(working_dir, nice_header + '.fasta')
        save_to_fasta(nice_header, seq, one_seq_fasta)
        sam_filename = os.path.join(working_dir, nice_header + '.sam')
        run_graphmap(one_seq_fasta, long_reads_fastq, sam_filename, graphmap_path, threads)

        # Copy the segment's SAM alignments to the final SAM file.
        alignment_count = 0
        one_seq_sam = open(sam_filename, 'r')
        for line in one_seq_sam:
            if not line.startswith('@') and line.split('\t', 3)[2] != '*':
                final_sam.write(line)
                alignment_count += 1
        one_seq_sam.close()

        if print_summary:
            nice_header = get_nice_header_and_len(header, seq, pad_length=longest_header) + ':'
            print(nice_header, int_to_str(alignment_count, read_count))
            sys.stdout.flush()

        # Clean up
        os.remove(one_seq_fasta)
        os.remove(sam_filename)

    final_sam.close()

def run_graphmap(fasta, long_reads_fastq, sam_file, graphmap_path, threads):
    '''
    This function runs GraphMap for the given inputs and produces a SAM file at the given location.
    '''
    command = [graphmap_path, '-r', fasta, '-d', long_reads_fastq, '-o',
               sam_file, '-Z', '-F', '1.0', '-t', str(threads), '-a', 'anchorgotoh']
    process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    _, _ = process.communicate()

    # Clean up.
    os.remove(fasta + '.gmidx')
    os.remove(fasta + '.gmidxsec')

def load_alignments(sam_filename, ref_fasta):
    '''
    This function returns a list of Alignment objects from the given SAM file.
    '''
    references = {get_nice_header(header): seq for header, seq in load_fasta(ref_fasta)}
    alignments = []
    sam_file = open(sam_filename, 'r')
    for line in sam_file:
        if not line.startswith('@') and line.split('\t', 3)[2] != '*':
            alignments.append(Alignment.from_sam(line.strip(), references))
    return alignments

def seqan_alignment_one_read_all_refs(ref_headers_and_seqs, read, print_summary, threads, k_size,
                                      band_size, allowed_length_discrepancy):
    '''
    Aligns a single read against all reference sequences using Seqan.
    '''
    # Prepare the c functions.
    c_lib_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'seqan_align.so')
    c_lib = ctypes.CDLL(c_lib_path)
    c_lib.semiGlobalAlign.argtypes = [ctypes.c_char_p, ctypes.c_char_p, ctypes.c_int, ctypes.c_int,
                                      ctypes.c_int, ctypes.c_int, ctypes.c_double]
    c_lib.semiGlobalAlign.restype = ctypes.c_void_p
    c_lib.free_c_string.argtypes = [ctypes.c_void_p]
    c_lib.free_c_string.restype = None

    alignments = []
    for header, seq in ref_headers_and_seqs:
        nice_header = get_nice_header(header)
        forward_alignment = run_one_seqan_alignment(nice_header, seq, read, False, k_size,
                                                    band_size, allowed_length_discrepancy, c_lib)
        reverse_alignment = run_one_seqan_alignment(nice_header, seq, read, True, k_size,
                                                    band_size, allowed_length_discrepancy, c_lib)
        if forward_alignment:
            alignments.append(forward_alignment)
        if reverse_alignment:
            alignments.append(reverse_alignment)






def run_one_seqan_alignment(ref_name, ref_seq, read, rev_comp, k_size, band_size,
                            allowed_length_discrepancy, c_lib):
    '''
    Runs a single alignment using Seqan. Returns either an Alignment object (if successful) or None
    (if not).
    '''
    if rev_comp:
        read_seq = reverse_complement(read.sequence)
    else:
        read_seq = read.sequence

    ptr = c_lib.semiGlobalAlign(read_seq, ref_seq, len(read_seq), len(ref_seq), k_size, band_size,
                                allowed_length_discrepancy)
    alignment_result = ctypes.cast(ptr, ctypes.c_char_p).value
    c_lib.free_c_string(ptr)

    alignment_result_parts = alignment_result.split(',')
    if len(alignment_result_parts) < 13:
        return None

    cigar = alignment_result_parts[0]
    read_start = int(alignment_result_parts[1])
    read_end = int(alignment_result_parts[2])
    ref_start = int(alignment_result_parts[3])
    ref_end = int(alignment_result_parts[4])
    alignment_length = int(alignment_result_parts[5])
    match_count = int(alignment_result_parts[6])
    mismatch_count = int(alignment_result_parts[7])
    ref_mismatch_positions = alignment_result_parts[8].split(";")
    insertion_count = int(alignment_result_parts[9])
    ref_insertion_positions = alignment_result_parts[10].split(";")
    deletion_count = int(alignment_result_parts[11])
    ref_deletion_positions = alignment_result_parts[12].split(";")
    edit_distance = int(alignment_result_parts[13])
    percent_identity = float(alignment_result_parts[14])
    milliseconds = int(alignment_result_parts[15])

    return Alignment(False, read.name, read_seq, ref_name, ref_seq, 0, rev_comp, cigar,
                     read.qualities, read_start, read_end, ref_start, ref_end,
                     match_count, mismatch_count, insertion_count, deletion_count,
                     alignment_length, edit_distance, percent_identity,
                     ref_mismatch_positions, ref_insertion_positions, ref_deletion_positions)

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


def simplify_ranges(ranges):
    '''
    Collapses overlapping ranges together.
    '''
    fixed_ranges = []
    for int_range in ranges:
        if int_range[0] > int_range[1]:
            fixed_ranges.append((int_range[1], int_range[0]))
        elif int_range[0] < int_range[1]:
            fixed_ranges.append(int_range)
    starts_ends = [(x[0], 1) for x in fixed_ranges]
    starts_ends += [(x[1], -1) for x in fixed_ranges]
    starts_ends.sort(key=lambda x: x[0])
    current_sum = 0
    cumulative_sum = []
    for start_end in starts_ends:
        current_sum += start_end[1]
        cumulative_sum.append((start_end[0], current_sum))
    prev_depth = 0
    start = 0
    combined = []
    for pos, depth in cumulative_sum:
        if prev_depth == 0:
            start = pos
        elif depth == 0:
            combined.append((start, pos))
        prev_depth = depth
    return combined

def range_is_contained(test_range, other_ranges):
    '''
    Returns True if test_range is entirely contained within any range in other_ranges.
    '''
    start, end = test_range
    for other_range in other_ranges:
        if other_range[0] <= start and other_range[1] >= end:
            return True
    return False


def write_sam_file(alignments, sam_filename):
    '''
    Writes the given alignments to a SAM file.
    '''
    sam_file = open(sam_filename, 'w')
    for alignment in alignments:
        sam_file.write(alignment.get_sam_line())
        sam_file.write('\n')
    sam_file.close()


def load_fasta(filename): # type: (str) -> list[tuple[str, str]]
    '''
    Returns the names and sequences for the given fasta file.
    '''
    fasta_seqs = []
    fasta_file = open(filename, 'r')
    name = ''
    sequence = ''
    for line in fasta_file:
        line = line.strip()
        if not line:
            continue
        if line[0] == '>': # Header line = start of new contig
            if name:
                fasta_seqs.append((name.split()[0], sequence))
                name = ''
                sequence = ''
            name = line[1:]
        else:
            sequence += line
    if name:
        fasta_seqs.append((name.split()[0], sequence))
    return fasta_seqs

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

def get_nice_header_and_len(header, seq, pad_length=0):
    '''
    Add the length in base pairs to the nice header. If there is a pad length, it will add spaces
    in between the header and length so things can line up nicely.
    '''
    part_1 = get_nice_header(header) + ' '
    part_2 = '(' + '{:,}'.format(len(seq)) + ' bp)'
    if len(part_1) + len(part_2) < pad_length:
        spaces = ' ' * (pad_length - len(part_1) - len(part_2))
    else:
        spaces = ''
    return part_1 + spaces + part_2

def save_to_fasta(header, sequence, filename):
    '''
    Saves the header/sequence to FASTA file.
    '''
    fasta = open(filename, 'w')
    fasta.write('>' + header + '\n')
    fasta.write(add_line_breaks_to_sequence(sequence, 60))
    fasta.close()

def save_reads_to_fastq(reads, fastq_filename):
    '''
    Writes the given reads to a FASTQ file.
    '''
    fastq = open(fastq_filename, 'w')
    for read in reads:
        fastq.write(read.get_fastq())

def add_line_breaks_to_sequence(sequence, length):
    '''
    Wraps sequences to the defined length.  All resulting sequences end in a line break.
    '''
    seq_with_breaks = ''
    while len(sequence) > length:
        seq_with_breaks += sequence[:length] + '\n'
        sequence = sequence[length:]
    if len(sequence) > 0:
        seq_with_breaks += sequence
        seq_with_breaks += '\n'
    return seq_with_breaks

def get_mean_and_st_dev_identity(alignments, limit_to_safe_alignments):
    '''
    This function returns the mean and standard deviation for the identities of the given
    alignments. If limit_to_safe_alignments is True, it only considers alignments that are
    entirely contained within contigs without having been extended.
    If there are 0 alignments, it returns 0 for both mean and std dev.
    If there is 1 alignment, it returns the real mean (the identity of that alignment) and a std
    dev of -1 (because we need two to get a std dev).
    '''
    identities = []
    for alignment in alignments:
        if limit_to_safe_alignments and \
           (alignment.extended or not alignment.is_whole_read()):
            continue
        identities.append(alignment.percent_identity)
    num = len(identities)
    if num == 0:
        return 0.0, 0.0
    mean = sum(identities) / num
    if num == 1:
        return mean, -1
    sum_squares = sum((x - mean) ** 2 for x in identities)
    st_dev = (sum_squares / (num - 1)) ** 0.5
    return mean, st_dev

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

def float_to_str(num, max_num=0):
    '''
    Converts a number to a string. Will add left padding based on the max value to ensure numbers
    align well.
    '''
    num_str = '%.1f' % num
    after_decimal = num_str.split('.')[1]
    num_str = int_to_str(int(num)) + '.' + after_decimal
    if max_num > 0:
        max_str = float_to_str(max_num)
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

def line_count(filename):
    '''
    Counts the lines in the given file.
    '''
    i = 0
    with open(filename) as file_to_count:
        for i, _ in enumerate(file_to_count):
            pass
    return i + 1

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



class LongRead(object):
    '''
    This class holds a long read, e.g. from PacBio or Oxford Nanopore.
    '''
    def __init__(self, name, sequence, qualities):
        self.name = name
        self.sequence = sequence
        self.qualities = qualities
        self.alignments = []

    def remove_conflicting_alignments(self):
        '''
        This function removes alignments from the read which are likely to be spurious. It sorts
        alignments by identity and works through them from highest identity to lowest identity,
        only keeping alignments that cover new parts of the read.
        It also uses an identity threshold to remove alignments with very poor identity.
        '''
        self.alignments = sorted(self.alignments, reverse=True,
                                 key=lambda x: (x.percent_identity, random.random()))
        kept_alignments = []
        read_ranges = []
        for alignment in self.alignments:
            read_range = alignment.read_start_end_positive_strand()
            if not range_is_contained(read_range, read_ranges):
                read_ranges.append(read_range)
                read_ranges = simplify_ranges(read_ranges)
                kept_alignments.append(alignment)
        self.alignments = kept_alignments

    def remove_low_id_alignments(self, id_threshold):
        '''
        This function removes alignments with identity below the cutoff.
        '''
        self.alignments = [x for x in self.alignments if x.percent_identity >= id_threshold]

    def get_fastq(self):
        '''
        Returns a string for the read in FASTQ format. It contains four lines and ends in a line
        break.
        '''
        return '@' + self.name + '\n' + \
               self.sequence + '\n' + \
               '+' + self.name + '\n' + \
               self.qualities + '\n'

    def get_descriptive_string(self):
        '''
        Returns a multi-line string that describes the read and its alignments.
        '''
        header = self.name + ' (' + str(len(self.sequence)) + ' bp)'
        line = '-' * len(header)
        description = header + '\n' + line + '\n'
        if not self.alignments:
            description += 'no alignments'
        else:
            description += '%.2f' % (100.0 * self.get_fraction_aligned()) + '% aligned\n'
            description += '\n'.join([str(x) for x in self.alignments])
        return description + '\n\n'

    def get_fraction_aligned(self):
        '''
        This function returns the fraction of the read which is covered by any of the read's
        alignments.
        '''
        read_ranges = [x.read_start_end_positive_strand() for x in self.alignments]
        read_ranges = simplify_ranges(read_ranges)
        aligned_length = sum([x[1] - x[0] for x in read_ranges])
        return aligned_length / len(self.sequence)



class Alignment(object):
    '''
    This class describes an alignment between a long read and a contig.
    '''
    def __init__(self, made_from_sam, read_name, full_read_sequence, ref_name, full_ref_sequence,
                 mapping_quality, rev_comp, cigar, read_qualities,
                 read_start_pos, read_end_pos, ref_start_pos, ref_end_pos,
                 match_count, mismatch_count, insertion_count, deletion_count,
                 alignment_length, edit_distance, percent_identity,
                 ref_mismatch_positions, ref_insertion_positions, ref_deletion_positions):
        # Save the many arguments to the class.
        self.read_name = read_name
        self.full_read_sequence = full_read_sequence
        self.ref_name = ref_name
        self.full_ref_sequence = full_ref_sequence
        self.mapping_quality = mapping_quality
        self.rev_comp = rev_comp
        self.cigar = cigar
        self.read_qualities = read_qualities
        self.read_start_pos = read_start_pos
        self.read_end_pos = read_end_pos
        self.ref_start_pos = ref_start_pos
        self.ref_end_pos = ref_end_pos
        self.match_count = match_count
        self.mismatch_count = mismatch_count
        self.insertion_count = insertion_count
        self.deletion_count = deletion_count
        self.alignment_length = alignment_length
        self.edit_distance = edit_distance
        self.percent_identity = percent_identity
        self.ref_mismatch_positions = ref_mismatch_positions
        self.ref_insertion_positions = ref_insertion_positions
        self.ref_deletion_positions = ref_deletion_positions

        self.cigar_parts = re.findall(r'\d+\w', self.cigar)
        self.read_length = len(self.full_read_sequence)
        self.ref_length = len(self.full_ref_sequence)

        # If the alignment was done by Seqan, then most of the other details were already done in
        # the C++ code. If it was done by GraphMap, then made_from_sam is True and various other
        # steps are necessary in this constructor.

        if made_from_sam:
            # Use the CIGAR to determine read/reference positions.
            self.ref_end_pos = self.ref_start_pos
            for cigar_part in self.cigar_parts:
                self.ref_end_pos += get_ref_shift_from_cigar_part(cigar_part)
            self.read_start_pos = self.get_start_soft_clips()
            self.read_end_pos = self.read_length - self.get_end_soft_clips()

        self.reference_end_gap = self.ref_length - self.ref_end_pos
        self.read_end_gap = self.get_end_soft_clips()

        if made_from_sam:
            # Extend the alignment so it is fully semi-global, reaching the end of the sequence.
            self.extended = False
            self.extend_alignment()

        # Get the aligned parts of the read and reference sequences.
        self.aligned_read_seq = self.full_read_sequence[self.read_start_pos:self.read_end_pos]
        self.aligned_ref_seq = self.full_ref_sequence[self.ref_start_pos:self.ref_end_pos]

        if made_from_sam:
            # Count matches, mismatches, insertions and deletions.
            # Insertions and deletions are counted per base. E.g. 5M3I4M has 3 insertions, not 1.
            self.tally_up_alignment()

    @classmethod
    def from_sam(cls, sam_line, references):
        '''
        Construct an Alignment using a SAM line.
        '''
        sam_parts = sam_line.split('\t')
        read_name = sam_parts[0].split('/')[0]
        sam_flag = int(sam_parts[1])
        rev_comp = bool(sam_flag & 0x10)
        ref_name = sam_parts[2]
        ref_start_pos = int(sam_parts[3]) - 1
        mapping_quality = int(sam_parts[4])
        cigar = sam_parts[5]
        full_read_sequence = sam_parts[9]
        read_qualities = sam_parts[10]
        return cls(True, read_name, full_read_sequence, ref_name, references[ref_name],
                   mapping_quality, rev_comp, cigar, read_qualities,
                   None, None, ref_start_pos, None,
                   None, None, None, None,
                   None, None, None,
                   None, None, None)

    def __repr__(self):
        if self.rev_comp:
            strand = '-'
        else:
            strand = '+'
        read_start, read_end = self.read_start_end_positive_strand()
        return self.read_name + ' (' + str(read_start) + '-' + str(read_end) + \
               ', strand: ' + strand + '), ' + self.ref_name + ' (' + \
               str(self.ref_start_pos) + '-' + str(self.ref_end_pos) + '), ' + \
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
        return self.ref_end_pos - self.ref_start_pos

    def get_read_to_ref_ratio(self):
        '''
        Returns the length ratio between the aligned parts of the read and reference.
        '''
        return self.get_alignment_length_read() / self.get_alignment_length_reference()

    def read_start_end_positive_strand(self):
        '''
        This function returns the read start/end coordinates for the positive strand of the read.
        For alignments on the positive strand, this is just the normal start/end. But for
        alignments on the negative strand, the coordinates are flipped to the other side.
        '''
        if not self.rev_comp:
            return self.read_start_pos, self.read_end_pos
        else:
            start = self.read_length - self.read_end_pos
            end = self.read_length - self.read_start_pos
            return start, end

    def extend_alignment(self):
        '''
        This function extends the alignment as much as possible in both directions so the alignment
        only terminates when it reaches the end of either the read or the reference.
        It does not actually perform the alignment - it just counts each alignment as a match. This
        means that very long extensions will probably result in terrible alignments, but that's
        okay because we'll filter alignments by quality later.
        '''
        missing_bases_at_start = min(self.read_start_pos, self.ref_start_pos)
        missing_bases_at_end = min(self.read_end_gap, self.reference_end_gap)

        if missing_bases_at_start:
            self.extended = True

            # Adjust the start of the reference.
            self.ref_start_pos -= missing_bases_at_start

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
            self.extended = True

            # Adjust the end of the reference.
            self.ref_end_pos += missing_bases_at_end
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

    def tally_up_alignment(self):
        '''
        Counts the matches, mismatches, indels and deletions. Also calculates the percent identity,
        which it does like BLAST: matches / alignment positions.
        '''
        # Reset any existing tallies.
        self.match_count = 0
        self.mismatch_count = 0
        self.insertion_count = 0
        self.deletion_count = 0
        self.percent_identity = 0.0
        self.ref_mismatch_positions = []
        self.ref_insertion_positions = []
        self.ref_deletion_positions = []

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
                self.insertion_count += cigar_count
                self.ref_insertion_positions += [ref_i + self.ref_start_pos] * cigar_count
                read_i += cigar_count
            elif cigar_type == 'D':
                self.deletion_count += cigar_count
                for i in xrange(cigar_count):
                    self.ref_deletion_positions.append(ref_i + self.ref_start_pos + i)
                ref_i += cigar_count
            else: # match/mismatch
                for _ in xrange(cigar_count):
                    read_base = self.aligned_read_seq[read_i]
                    ref_base = self.aligned_ref_seq[ref_i]
                    if read_base == ref_base:
                        self.match_count += 1
                    else:
                        self.mismatch_count += 1
                        self.ref_mismatch_positions.append(ref_i + self.ref_start_pos)
                    read_i += 1
                    ref_i += 1
            align_i += cigar_count
        self.percent_identity = 100.0 * self.match_count / align_i

    def get_sam_line(self):
        '''
        Returns a SAM alignment line.
        '''
        edit_distance = self.mismatch_count + self.insertion_count + self.deletion_count
        sam_flag = 0 #TEMP
        return '\t'.join([self.read_name, str(sam_flag), self.ref_name,
                          str(self.ref_start_pos + 1), str(self.mapping_quality), self.cigar,
                          '*', '0', '0', self.full_read_sequence, self.read_qualities,
                          'NM:i:' + str(edit_distance)])

    def is_whole_read(self):
        '''
        Returns True if the alignment covers the entirety of the read.
        '''
        return self.read_start_pos == 0 and self.read_end_gap == 0



if __name__ == '__main__':
    main()
