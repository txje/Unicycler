#!/usr/bin/env python3
"""
Tool for comparing the assemblies made by Unicycler, SPAdes, hybridSPAdes, ABySS, Cerulean and
npScarf.

Can use ART and PBSIM to generate synthetic reads and then compare assemblies back to the
reference using QUAST. Or it can take real reads and a reference for those reads.

Author: Ryan Wick
email: rrwick@gmail.com
"""

import random
import os
import subprocess
import argparse
import sys
import time
import shutil
import gzip
import datetime


def main():
    set_up_env_var()
    args = get_args()

    if args.reference:
        ref_name = get_reference_name_from_filename(args.reference)
        scaled_ref_length = get_scaled_ref_length(args)
    else:
        ref_name = args.genome_name
        scaled_ref_length = args.genome_size

    # Make a directory for this ref, if necessary.
    ref_dir = os.path.abspath(ref_name)
    if not os.path.exists(ref_dir):
        os.makedirs(ref_dir)
    os.chdir(ref_dir)

    if args.real_short_1 and args.real_short_2 and args.real_long:
        real_reads(args, scaled_ref_length, ref_dir)

    elif args.short_depth and args.long_accs and args.long_lens and args.model_qc:
        simulated_reads(args, scaled_ref_length, ref_dir)

    else:
        quit_with_error('Options not compatible with either real or simulated reads')


def get_args():
    """
    Specifies the command line arguments required by the script.
    """
    parser = argparse.ArgumentParser(description='Assembly tester')

    # Used for both real and simulated reads.
    parser.add_argument('--reference', type=str,
                        help='The reference genome to shred and reassemble')

    # Used for real reads.
    parser.add_argument('--real_short_1', type=str,
                        help='Real short reads, forward half of pair')
    parser.add_argument('--real_short_2', type=str,
                        help='Real short reads, second half of pair')
    parser.add_argument('--real_long', type=str,
                        help='Real long reads')

    # Used for real reads (no reference).
    parser.add_argument('--genome_name', type=str,
                        help='The name of the genome being assembled')
    parser.add_argument('--genome_size', type=int,
                        help='The size of the whole genome, scaled for relative depths')

    # Used for simulated reads.
    parser.add_argument('--short_depth', type=float, default=50.0,
                        help='Base read depth for fake short reads')
    parser.add_argument('--long_accs', type=str, default='75',
                        help='Mean accuracies for long reads (comma delimited)')
    parser.add_argument('--long_lens', type=str, default='10000',
                        help='Mean lengths for long reads (comma delimited)')
    parser.add_argument('--model_qc', type=str, default='model_qc_clr',
                        help='Model QC file for pbsim')
    parser.add_argument('--rotation_count', type=int, default=20, required=False,
                        help='The number of times to run read simulators with random start '
                             'positions')

    # Used for both real and simulated reads.
    parser.add_argument('--threads', type=int, required=False, default=8,
                        help='Number of CPU threads')
    parser.add_argument('--long_depths', type=str, default='0.25,0.5,1,2,4,8,16',
                        help='Maximum read depth for long reads')
    parser.add_argument('--iterations', type=int, default=10,
                        help='Number of replicate tests')
    parser.add_argument('--no_cerulean', action='store_true',
                        help='Skips Cerulean and ABySS assemblies')
    parser.add_argument('--no_npscarf', action='store_true',
                        help='Skips npScarf assemblies')
    parser.add_argument('--no_spades', action='store_true',
                        help='Skips SPAdes and hybridSPAdes assemblies')
    parser.add_argument('--no_unicycler', action='store_true',
                        help='Skips Unicycler assemblies')

    args = parser.parse_args()
    if args.reference:
        args.reference = os.path.abspath(args.reference)
    if args.real_short_1:
        args.real_short_1 = os.path.abspath(args.real_short_1)
    if args.real_short_2:
        args.real_short_2 = os.path.abspath(args.real_short_2)
    if args.real_long:
        args.real_long = os.path.abspath(args.real_long)
    if args.model_qc:
        args.model_qc = os.path.abspath(args.model_qc)
    return args


def real_reads(args, scaled_ref_length, ref_dir):
    """
    Runs tests using real reads and comparing to a reference.
    """
    using_ref = args.reference is not None

    if using_ref:
        print_with_timestamp('Running in real read mode (with reference): ' + args.real_short_1
                             + ', ' + args.real_short_2 + ', ' + args.real_long)
    else:
        print_with_timestamp('Running in real read mode (with reference): ' + args.real_short_1
                             + ', ' + args.real_short_2 + ', ' + args.real_long)

    long_reads = load_fastq(args.real_long, '')

    if using_ref:
        # Create and move into a directory for the specific long read set being used.
        long_read_set_name = get_reference_name_from_filename(args.real_long)
        long_read_set_dir = os.path.abspath(long_read_set_name)
        if not os.path.exists(long_read_set_dir):
            os.makedirs(long_read_set_dir)
        os.chdir(long_read_set_dir)
        ref_dir = long_read_set_dir
        quast_results, simple_quast_results = create_quast_results_tables(args)
    else:
        quast_results, simple_quast_results = None, None

    # Create and move into a directory for the short-only and all-long assemblies.
    dir_name, dir_num = get_next_available_set_number(ref_dir, 1)
    iter_dir = os.path.join(ref_dir, dir_name)
    print()
    print_with_timestamp('Changing to new directory: ' + str(iter_dir))
    os.chdir(iter_dir)

    abyss_dir, spades_no_long_dir, unicycler_no_long_dir, \
        abyss_time, spades_no_long_time, unicycler_no_long_time = \
        assemble_short_reads_only(args, args.real_short_1, args.real_short_2, quast_results,
                                  simple_quast_results)

    for i in range(args.iterations):

        # Create and move into a directory for this iteration.
        dir_name, dir_num = get_next_available_set_number(ref_dir, 1)
        iter_dir = os.path.join(ref_dir, dir_name)
        print()
        print_with_timestamp('Changing to new directory: ' + str(iter_dir))
        os.chdir(iter_dir)

        depths = [float(x) for x in args.long_depths.split(',')]
        for subsampled_depth in depths:

            subsampled_reads, subsampled_filename, subsampled_long_read_depth = \
                subsample_long_reads(args, long_reads, subsampled_depth, scaled_ref_length, False)

            assemble_long_reads(args, args.real_short_1, args.real_short_2, subsampled_reads,
                                subsampled_filename, len(subsampled_reads),
                                subsampled_long_read_depth, quast_results, simple_quast_results,
                                unicycler_no_long_dir, spades_no_long_dir, abyss_dir,
                                unicycler_no_long_time, spades_no_long_time, abyss_time)


def simulated_reads(args, scaled_ref_length, ref_dir):
    """
    Runs tests using fake reads.
    """
    quast_results, simple_quast_results = create_quast_results_tables(args)

    # Set up each read accuracy-length combination.
    accuracies = [float(x) for x in args.long_accs.split(',')]
    lengths = [int(x) for x in args.long_lens.split(',')]
    accuracies_and_lengths = []
    for accuracy in accuracies:
        accuracies_and_lengths += [(accuracy, length) for length in lengths]

    acc_len_str = ', '.join(str(x) + '% + ' + str(y) + ' bp' for x, y in accuracies_and_lengths)
    print_with_timestamp('Running in simulated read mode with these accuracy and length '
                         'combinations: ' + acc_len_str)

    dir_num = 1
    for i in range(args.iterations):

        # Create and move into a directory for this iteration.
        dir_name, dir_num = get_next_available_set_number(ref_dir, dir_num)
        iter_dir = os.path.join(ref_dir, dir_name)
        print()
        print_with_timestamp('Changing to new directory: ' + str(iter_dir))
        os.chdir(iter_dir)

        # Make new short reads every time through this loop.
        short_1, short_2 = make_fake_short_reads(args)

        abyss_dir, spades_no_long_dir, unicycler_no_long_dir, \
            abyss_time, spades_no_long_time, unicycler_no_long_time = \
            assemble_short_reads_only(args, short_1, short_2, quast_results, simple_quast_results)

        # Shuffle the acc/len order so if repeated runs are cut short I'm not always getting the
        # first combination.
        random.shuffle(accuracies_and_lengths)

        for accuracy, length in accuracies_and_lengths:
            args.long_acc, args.long_len = accuracy, length
            depths = [float(x) for x in args.long_depths.split(',')]
            for depth in depths:
                long_filename, long_reads = make_fake_long_reads(args, depth)
                long_read_count = len(long_reads)
                long_read_depth = get_long_read_depth(long_reads, scaled_ref_length)

                assemble_long_reads(args, short_1, short_2, long_reads, long_filename,
                                    long_read_count, long_read_depth, quast_results,
                                    simple_quast_results, unicycler_no_long_dir,
                                    spades_no_long_dir, abyss_dir, unicycler_no_long_time,
                                    spades_no_long_time, abyss_time)


def assemble_short_reads_only(args, short_1, short_2, quast_results, simple_quast_results):
    """
    Run ABySS, SPAdes and Unicycler on short read data alone.
    """
    args.long_acc, args.long_len = 0.0, 0
    if not args.no_cerulean:
        abyss_dir, abyss_time = run_abyss(short_1, short_2, args, quast_results,
                                          simple_quast_results)
    else:
        abyss_dir = None
        abyss_time = 0.0

    if not args.no_spades:
        spades_no_long_dir, spades_time = run_regular_spades(short_1, short_2, args, quast_results,
                                                             simple_quast_results)
    else:
        spades_no_long_dir = None
        spades_time = 0.0

    if not args.no_unicycler:
        unicycler_no_long_dir, unicycler_time = \
            run_unicycler(args, short_1, short_2, None, 0, 0.0, quast_results,
                          simple_quast_results, 'conservative')
        run_unicycler(args, short_1, short_2, None, 0, 0.0, quast_results,
                      simple_quast_results, 'normal', unicycler_no_long_dir,
                      replacement_time=unicycler_time)
        run_unicycler(args, short_1, short_2, None, 0, 0.0, quast_results,
                      simple_quast_results, 'bold', unicycler_no_long_dir,
                      replacement_time=unicycler_time)
    else:
        unicycler_no_long_dir = None
        unicycler_time = 0.0

    return abyss_dir, spades_no_long_dir, unicycler_no_long_dir, abyss_time, spades_time, \
        unicycler_time


def assemble_long_reads(args, short_1, short_2, long_reads, long_filename, long_read_count,
                        long_read_depth, quast_results, simple_quast_results,
                        unicycler_no_long_dir, spades_no_long_dir, abyss_dir,
                        unicycler_no_long_time, spades_no_long_time, abyss_time):
    # Run Unicycler on the full set of long reads - will be the source of alignments for
    # subsampled Unicycler runs.
    if not args.no_unicycler:
        unicycler_all_long_dir, unicycler_time = \
            run_unicycler(args, short_1, short_2, long_filename, long_read_count,
                          long_read_depth, quast_results, simple_quast_results, 'conservative',
                          unicycler_no_long_dir, extra_time=unicycler_no_long_time)
        run_unicycler(args, short_1, short_2, long_filename, long_read_count,
                      long_read_depth, quast_results, simple_quast_results, 'normal',
                      unicycler_no_long_dir, unicycler_all_long_dir,
                      replacement_time=unicycler_time)
        run_unicycler(args, short_1, short_2, long_filename, long_read_count,
                      long_read_depth, quast_results, simple_quast_results, 'bold',
                      unicycler_no_long_dir, unicycler_all_long_dir,
                      replacement_time=unicycler_time)

    else:
        unicycler_all_long_dir = None

    # Run hybridSPAdes on the full set of long reads.
    if not args.no_spades:
        run_hybrid_spades(short_1, short_2, long_filename, long_read_count, long_read_depth,
                          args, quast_results, simple_quast_results, spades_no_long_dir)

    # Run npScarf on the full set of long reads - will be the source of alignments for
    # subsampled npScarf runs.
    if not args.no_npscarf:
        np_scarf_all_long_dir = run_np_scarf(long_filename, long_reads, long_read_count,
                                             long_read_depth, args, quast_results,
                                             simple_quast_results, spades_no_long_dir,
                                             extra_time=spades_no_long_time)
    else:
        np_scarf_all_long_dir = None

    # Run Cerulean on the full set of long reads - will be the source of alignments for
    # subsampled Cerulean runs.
    if not args.no_cerulean:
        cerulean_all_long_dir = run_cerulean(long_filename, long_reads, long_read_count,
                                             long_read_depth, args, quast_results,
                                             simple_quast_results, abyss_dir, extra_time=abyss_time)
    else:
        cerulean_all_long_dir = None

    return unicycler_all_long_dir, np_scarf_all_long_dir, cerulean_all_long_dir


def subsample_long_reads(args, long_reads, subsampled_depth, scaled_ref_length, include_acc_len):
    read_indices = list(range(len(long_reads)))
    subsampled_reads = []
    random.shuffle(read_indices)
    for i in read_indices:
        subsampled_reads.append(long_reads[i])
        subsampled_long_read_depth = get_long_read_depth(subsampled_reads, scaled_ref_length)
        if subsampled_long_read_depth == subsampled_depth:
            break
        elif subsampled_long_read_depth > subsampled_depth:
            last_read = subsampled_reads.pop()
            under_shot_depth = get_long_read_depth(subsampled_reads, scaled_ref_length)
            depth_from_last_read = subsampled_long_read_depth - under_shot_depth
            required_depth = subsampled_depth - under_shot_depth
            fraction_last_read_needed = required_depth / depth_from_last_read
            last_read_length = len(last_read[1])
            trimmed_last_read_length = int(round(last_read_length * fraction_last_read_needed))
            trimmed_last_read = (last_read[0],
                                 last_read[1][:trimmed_last_read_length],
                                 last_read[2],
                                 last_read[3][:trimmed_last_read_length])
            subsampled_reads.append(trimmed_last_read)
            break

    subsampled_long_read_depth = get_long_read_depth(subsampled_reads, scaled_ref_length)
    if include_acc_len:
        subsampled_filename = 'long_subsampled_' + str(args.long_acc) + '_' + \
                              str(args.long_len) + '_' + str(subsampled_depth)
    else:
        subsampled_filename = 'long_subsampled_' + str(subsampled_depth)
    subsampled_filename = subsampled_filename.replace('.', '_')
    file_index = 1
    full_subsampled_filename = ''
    while True:
        if file_index == 1:
            full_subsampled_filename = os.path.abspath(subsampled_filename + '.fastq')
        else:
            full_subsampled_filename = os.path.abspath(subsampled_filename + '_' +
                                                       str(file_index) + '.fastq')
        if not os.path.isfile(full_subsampled_filename + '.gz'):
            break
        else:
            file_index += 1

    print()
    print_with_timestamp('Subsampling to ' + str(subsampled_depth) + 'x')
    save_long_reads_to_fastq(subsampled_reads, full_subsampled_filename)
    try:
        subprocess.check_output(['gzip', full_subsampled_filename], stderr=subprocess.STDOUT)
    except subprocess.CalledProcessError as e:
        quit_with_error('gzip encountered an error:\n' + e.output.decode())
    full_subsampled_filename += '.gz'

    return subsampled_reads, full_subsampled_filename, subsampled_long_read_depth


class AssemblyError(Exception):
    pass


def make_fake_short_reads(args):
    """
    Runs ART to generate fake Illumina reads. Runs ART separate for each sequence in the reference
    file (to control relative depth) and at multiple sequence rotations (to ensure circular
    assembly).
    """
    read_filename_1 = os.path.abspath('short_1.fastq')
    read_filename_2 = os.path.abspath('short_2.fastq')

    if os.path.isfile(read_filename_1):
        os.remove(read_filename_1)
    if os.path.isfile(read_filename_2):
        os.remove(read_filename_2)
    if os.path.isfile(read_filename_1 + '.gz'):
        os.remove(read_filename_1 + '.gz')
    if os.path.isfile(read_filename_2 + '.gz'):
        os.remove(read_filename_2 + '.gz')

    print_with_timestamp('Generating synthetic short reads')

    references = load_fasta(args.reference)
    relative_depths = get_relative_depths(args.reference)

    # This will hold all simulated short reads. Each read is a list of 8 strings: the first four are
    # for the first read in the pair, the second four are for the second.
    short_read_pairs = []

    read_prefix = 1  # Used to prevent duplicate read names.
    for i, ref in enumerate(references):

        short_depth = relative_depths[i] * args.short_depth

        ref_seq = ref[1]
        circular = ref[3]

        if circular:
            short_depth_per_rotation = short_depth / args.rotation_count

            for j in range(args.rotation_count):

                # Randomly rotate the sequence.
                random_start = random.randint(0, len(ref_seq) - 1)
                rotated = ref_seq[random_start:] + ref_seq[:random_start]

                # Save the rotated sequence to FASTA.
                temp_fasta_filename = 'temp_rotated.fasta'
                temp_fasta = open(temp_fasta_filename, 'w')
                temp_fasta.write('>' + ref[0] + '\n')
                temp_fasta.write(rotated + '\n')
                temp_fasta.close()

                short_read_pairs += run_art(temp_fasta_filename, short_depth_per_rotation,
                                            str(read_prefix))
                os.remove(temp_fasta_filename)
                read_prefix += 1

        else:  # linear
            temp_fasta_filename = 'temp.fasta'
            temp_fasta = open(temp_fasta_filename, 'w')
            temp_fasta.write('>' + ref[0] + '\n')
            temp_fasta.write(ref_seq + '\n')
            temp_fasta.close()

            short_read_pairs += run_art(temp_fasta_filename, short_depth, str(read_prefix))
            os.remove(temp_fasta_filename)
            read_prefix += 1

    random.shuffle(short_read_pairs)
    reads_1 = open(read_filename_1, 'w')
    reads_2 = open(read_filename_2, 'w')
    for read_pair in short_read_pairs:
        reads_1.write(read_pair[0])
        reads_1.write('\n')
        reads_1.write(read_pair[1])
        reads_1.write('\n')
        reads_1.write(read_pair[2])
        reads_1.write('\n')
        reads_1.write(read_pair[3])
        reads_1.write('\n')
        reads_2.write(read_pair[4])
        reads_2.write('\n')
        reads_2.write(read_pair[5])
        reads_2.write('\n')
        reads_2.write(read_pair[6])
        reads_2.write('\n')
        reads_2.write(read_pair[7])
        reads_2.write('\n')
    reads_1.close()
    reads_2.close()

    # Gzip the reads to save space.
    try:
        subprocess.check_output(['gzip', read_filename_1], stderr=subprocess.STDOUT)
        subprocess.check_output(['gzip', read_filename_2], stderr=subprocess.STDOUT)
    except subprocess.CalledProcessError as e:
        quit_with_error('gzip encountered an error:\n' + e.output.decode())
    read_filename_1 += '.gz'
    read_filename_2 += '.gz'
    print_with_timestamp(read_filename_1)
    print_with_timestamp(read_filename_2)

    return read_filename_1, read_filename_2


def make_fake_long_reads(args, long_depth):
    """
    Runs pbsim to generate fake long reads. Runs pbsim separate for each sequence in the reference
    file (to control relative depth) and at multiple sequence rotations (to ensure circular
    assembly).
    """
    long_filename = 'long_' + '{0:g}'.format(args.long_acc) + '_' + str(args.long_len) + '_' + \
        '{0:g}'.format(long_depth) + 'x'
    long_filename = long_filename.replace('.', '_')
    long_filename = os.path.abspath(long_filename + '.fastq')

    print()
    print_with_timestamp('Generating synthetic long reads: ' + str(args.long_acc) +
                         '% accuracy, ' + str(args.long_len) + ' bp, ' +
                         '{0:g}'.format(long_depth) + 'x')

    references = load_fasta(args.reference)
    relative_depths = get_relative_depths(args.reference)

    # This will hold all simulated long reads. Each read is a list of 4 strings: one for each
    # line of the FASTQ.
    long_reads = []

    read_prefix = 1  # Used to prevent duplicate read names.
    for i, ref in enumerate(references):
        ref_seq = ref[1]
        ref_len = len(ref_seq)
        ref_long_depth = relative_depths[i] * long_depth

        ref_seq = ref[1]
        circular = ref[3]

        if circular:
            rotation_count = args.rotation_count
            long_depth_per_rotation = ref_long_depth / rotation_count

            while long_depth_per_rotation * ref_len < args.long_len and rotation_count > 1:
                rotation_count = max(1, rotation_count // 2)
                long_depth_per_rotation = ref_long_depth / rotation_count

            for j in range(rotation_count):

                # Randomly rotate the sequence.
                random_start = random.randint(0, ref_len - 1)
                rotated = ref_seq[random_start:] + ref_seq[:random_start]

                # If the sequence is very short compared to the read length, then we duplicate the
                # sequence. This is because PBSIM does weird stuff when you, for example, ask for
                # 25 kb reads from an 8 kb sequence.
                copy_adjusted_depth = long_depth_per_rotation
                copy_adjusted_sequence = rotated
                while len(copy_adjusted_sequence) < 5 * args.long_len:
                    copy_adjusted_sequence += copy_adjusted_sequence
                    copy_adjusted_depth /= 2.0

                # Save the rotated sequence to FASTA.
                temp_fasta_filename = 'temp_rotated.fasta'
                temp_fasta = open(temp_fasta_filename, 'w')
                temp_fasta.write('>' + ref[0] + '\n')
                temp_fasta.write(copy_adjusted_sequence + '\n')
                temp_fasta.close()

                long_reads += run_pbsim(temp_fasta_filename, copy_adjusted_depth, args,
                                        str(read_prefix), ref_len)
                os.remove(temp_fasta_filename)
                read_prefix += 1

        else:  # linear
            temp_fasta_filename = 'temp.fasta'
            temp_fasta = open(temp_fasta_filename, 'w')
            temp_fasta.write('>' + ref[0] + '\n')
            temp_fasta.write(ref_seq + '\n')
            temp_fasta.close()

            long_reads += run_pbsim(temp_fasta_filename, ref_long_depth, args, str(read_prefix),
                                    ref_len)
            os.remove(temp_fasta_filename)
            read_prefix += 1

    random.shuffle(long_reads)
    save_long_reads_to_fastq(long_reads, long_filename)

    # Gzip the reads to save space.
    try:
        subprocess.check_output(['gzip', long_filename], stderr=subprocess.STDOUT)
    except subprocess.CalledProcessError as e:
        quit_with_error('gzip encountered an error:\n' + e.output.decode())
    long_filename += '.gz'

    print_with_timestamp(long_filename)
    return long_filename, long_reads


def save_long_reads_to_fastq(long_reads, fastq):
    fastq_file = open(fastq, 'w')
    for read in long_reads:
        fastq_file.write(read[0] + '\n')
        fastq_file.write(read[1] + '\n')
        fastq_file.write(read[2] + '\n')
        fastq_file.write(read[3] + '\n')
    fastq_file.close()


def get_scaled_ref_length(args):
    references = load_fasta(args.reference)
    relative_depths = get_relative_depths(args.reference)
    if len(references) != len(relative_depths):
        quit_with_error('you must provide exactly one relative depth for each reference sequence')
    total_length = 0
    for i, ref in enumerate(references):
        total_length += relative_depths[i] * len(ref[1])
    return total_length


def get_long_read_depth(long_reads, scaled_ref_length):
    total_read_length = sum(len(x[1]) for x in long_reads)
    return total_read_length / scaled_ref_length


def create_subsampled_sam(full_sam_file, subsampled_sam_file, subsampled_long_reads):
    subsampled_read_names = set()
    for read in subsampled_long_reads:
        subsampled_read_names.add(read[0][1:])
    with open(full_sam_file, 'rt') as full_sam:
        with open(subsampled_sam_file, 'w') as subsampled_sam:
            for line in full_sam:
                if line.startswith('@'):
                    subsampled_sam.write(line)
                else:
                    read_name = line.split('\t', 1)[0]
                    if read_name in subsampled_read_names:
                        subsampled_sam.write(line)


def create_subsampled_m4(full_m4_file, subsampled_m4_file, subsampled_long_reads):
    subsampled_read_names = set()
    for read in subsampled_long_reads:
        subsampled_read_names.add(read[0][1:])
    with open(full_m4_file, 'rt') as full_m4:
        with open(subsampled_m4_file, 'w') as subsampled_m4:
            for line in full_m4:
                read_name = line.split()[0]
                if read_name in subsampled_read_names:
                    subsampled_m4.write(line)


def run_art(input_fasta, depth, read_prefix):
    """
    Runs ART with some fixed settings: 125 bp paired-end reads, 400 bp fragments, HiSeq 2500.
    Returns reads as list of list of strings.
    """
    art_command = ['art_illumina',
                   '--seqSys', 'HS25',
                   '--in', input_fasta,
                   '--len', '125',
                   '--mflen', '400',
                   '--sdev', '60',
                   '--fcov', str(depth),
                   '--out', 'art_output']
    try:
        subprocess.check_output(art_command, stderr=subprocess.STDOUT)
    except subprocess.CalledProcessError as e:
        quit_with_error('ART encountered an error:\n' + e.output.decode())

    output_fastq_1_filename = 'art_output1.fq'
    output_fastq_2_filename = 'art_output2.fq'
    try:
        with open(output_fastq_1_filename, 'rt') as f:
            fastq_1_lines = f.read().splitlines()
        with open(output_fastq_2_filename, 'rt') as f:
            fastq_2_lines = f.read().splitlines()
        pair_count = int(len(fastq_1_lines) / 4)
    except FileNotFoundError:
        pair_count = 0
        quit_with_error('Could not find ART output read files')

    os.remove('art_output1.fq')
    os.remove('art_output2.fq')
    os.remove('art_output1.aln')
    os.remove('art_output2.aln')

    read_pairs = []
    i = 0
    for _ in range(pair_count):
        name_1 = '@' + read_prefix + '_' + fastq_1_lines[i].split('-')[1]
        name_2 = '@' + read_prefix + '_' + fastq_2_lines[i].split('-')[1]
        seq_1 = fastq_1_lines[i + 1]
        seq_2 = fastq_2_lines[i + 1]
        qual_1 = fastq_1_lines[i + 3]
        qual_2 = fastq_2_lines[i + 3]

        read_pairs.append((name_1, seq_1, '+', qual_1, name_2, seq_2, '+', qual_2))
        i += 4

    return read_pairs


def run_pbsim(input_fasta, depth, args, read_prefix, seq_length):

    # Reduce the max and mean sequence lengths when the reference sequence is short.
    max_length = min(100000, seq_length)
    mean_length = min(args.long_len, seq_length)

    pbsim_command = ['pbsim',
                     '--depth', str(depth),
                     '--length-min', '100',
                     '--length-max', str(max_length),
                     '--accuracy-min', '0.5',
                     '--accuracy-max', '1.0',
                     '--model_qc', args.model_qc,
                     '--length-mean', str(mean_length),
                     '--length-sd', str(mean_length / 10),
                     '--accuracy-mean', str(args.long_acc / 100.0),
                     '--accuracy-sd', '0.02',
                     '--difference-ratio', '10:40:30',
                     '--seed', str(random.randint(0, 1000000)),
                     input_fasta]
    try:
        subprocess.check_output(pbsim_command, stderr=subprocess.STDOUT)
    except subprocess.CalledProcessError as e:
        quit_with_error('pbsim encountered an error:\n' + e.output.decode())

    reads = load_fastq('sd_0001.fastq', read_prefix)
    os.remove('sd_0001.fastq')
    os.remove('sd_0001.maf')
    os.remove('sd_0001.ref')

    return reads


def get_compression_type(filename):
    """
    Attempts to guess the compression (if any) on a file using the first few bytes.
    http://stackoverflow.com/questions/13044562
    """
    magic_dict = {'gz': (b'\x1f', b'\x8b', b'\x08'),
                  'bz2': (b'\x42', b'\x5a', b'\x68'),
                  'zip': (b'\x50', b'\x4b', b'\x03', b'\x04')}
    max_len = max(len(x) for x in magic_dict)

    unknown_file = open(filename, 'rb')
    file_start = unknown_file.read(max_len)
    unknown_file.close()
    compression_type = 'plain'
    for filetype, magic_bytes in magic_dict.items():
        if file_start.startswith(magic_bytes):
            compression_type = filetype
    if compression_type == 'bz2':
        quit_with_error('cannot use bzip2 format - use gzip instead')
    if compression_type == 'zip':
        quit_with_error('cannot use zip format - use gzip instead')
    return compression_type


def load_fastq(fastq_filename, read_prefix):
    if get_compression_type(fastq_filename) == 'gz':
        open_func = gzip.open
    else:  # plain text
        open_func = open
    reads = []
    with open_func(fastq_filename, 'rt') as fastq:
        for line in fastq:
            name = line.strip()[1:].split()[0]
            if read_prefix:
                name = read_prefix + '_' + name
            name = '@' + name
            sequence = next(fastq).strip()
            _ = next(fastq)
            qualities = next(fastq).strip()
            reads.append((name, sequence, '+', qualities))
    return reads


def run_abyss(short_1, short_2, args, all_quast_results, simple_quast_results):
    """
    Runs ABySS on the short reads to compare with SPAdes.
    """
    run_name, abyss_dir = get_run_name_and_run_dir_name('ABySS', args.reference, 0.0, args)
    abyss_contigs = os.path.join(abyss_dir, 'run-contigs.fa')
    abyss_scaffolds = os.path.join(abyss_dir, 'run-scaffolds.fa')
    starting_path = os.path.abspath('.')

    abyss_start_time = time.time()
    try:
        print_with_timestamp('Running ' + run_name)
        if not os.path.exists(abyss_dir):
            os.makedirs(abyss_dir)
        shutil.copyfile(short_1, os.path.join(abyss_dir, 'reads_1.fastq.gz'))
        shutil.copyfile(short_2, os.path.join(abyss_dir, 'reads_2.fastq.gz'))
        os.chdir(abyss_dir)
        try:
            subprocess.check_output(['gunzip', 'reads_1.fastq.gz'], stderr=subprocess.STDOUT)
            subprocess.check_output(['gunzip', 'reads_2.fastq.gz'], stderr=subprocess.STDOUT)
        except subprocess.CalledProcessError as e:
            quit_with_error('gzip encountered an error:\n' + e.output.decode())

        abyss_command = "abyss-pe k=64 j=" + str(args.threads) + \
                        " in='reads_1.fastq reads_2.fastq' name=run"
        print_with_timestamp(abyss_command)
        try:
            abyss_out = subprocess.check_output(abyss_command, stderr=subprocess.STDOUT, shell=True)
            with open(os.path.join('abyss.out'), 'wb') as f:
                f.write(abyss_out)
        except subprocess.CalledProcessError as e:
            print_with_timestamp('ABySS encountered an error:\n' + e.output.decode())
            raise AssemblyError

    except AssemblyError:
        abyss_scaffolds = None
        abyss_contigs = None

    abyss_time = time.time() - abyss_start_time
    os.chdir(starting_path)

    run_quast(abyss_scaffolds, args, all_quast_results, simple_quast_results, 'ABySS scaffolds',
              0, 0.0, abyss_time, abyss_time, run_name, abyss_dir)
    run_quast(abyss_contigs, args, all_quast_results, simple_quast_results, 'ABySS contigs',
              0, 0.0, abyss_time, abyss_time, run_name, abyss_dir)

    return os.path.abspath(abyss_dir), abyss_time


def run_regular_spades(short_1, short_2, args, all_quast_results, simple_quast_results):
    """
    Runs SPAdes with only short reads (i.e. not hybridSPAdes).
    """
    run_name, spades_dir = get_run_name_and_run_dir_name('SPAdes', args.reference, 0.0, args)
    spades_scaffolds = os.path.join(spades_dir, 'scaffolds.fasta')
    spades_contigs = os.path.join(spades_dir, 'contigs.fasta')
    spades_before_rr = os.path.join(spades_dir, 'before_rr.fasta')

    spades_start_time = time.time()
    try:
        print_with_timestamp('Running ' + run_name)
        if not os.path.exists(spades_dir):
            os.makedirs(spades_dir)
        spades_command = ['spades.py',
                          '-1', short_1,
                          '-2', short_2,
                          '--careful',
                          '--threads', str(args.threads),
                          '-o', spades_dir]
        print_with_timestamp(' '.join(spades_command))
        try:
            spades_out = subprocess.check_output(spades_command, stderr=subprocess.STDOUT)
            with open(os.path.join(spades_dir, 'spades.out'), 'wb') as f:
                f.write(spades_out)
        except subprocess.CalledProcessError as e:
            print_with_timestamp('SPAdes encountered an error:\n' + e.output.decode())
            raise AssemblyError

    except AssemblyError:
        spades_scaffolds = None
        spades_contigs = None
        spades_before_rr = None

    spades_time = time.time() - spades_start_time

    run_quast(spades_scaffolds, args, all_quast_results, simple_quast_results, 'SPAdes scaffolds',
              0, 0.0, spades_time, spades_time, run_name, spades_dir)
    run_quast(spades_contigs, args, all_quast_results, simple_quast_results, 'SPAdes contigs',
              0, 0.0, spades_time, spades_time, run_name, spades_dir)
    run_quast(spades_before_rr, args, all_quast_results, simple_quast_results, 'SPAdes before_rr',
              0, 0.0, spades_time, spades_time, run_name, spades_dir)

    clean_up_spades_dir(spades_dir)
    return os.path.abspath(spades_dir), spades_time


def run_hybrid_spades(short_1, short_2, long, long_count, long_depth, args, all_quast_results,
                      simple_quast_results, spades_no_long_dir):
    """
    Runs hybridSPAdes using the --nanopore option.
    """
    run_name, spades_dir = get_run_name_and_run_dir_name('SPAdes', args.reference, long_depth, args)
    spades_assembly = os.path.join(spades_dir, 'scaffolds.fasta')

    # If we can find corrected reads from a previous SPAdes run, use those to save time.
    corrected_1, corrected_2, corrected_u = '', '', ''
    corrected_dir = os.path.join(spades_no_long_dir, 'corrected')
    for item in os.listdir(corrected_dir):
        if item.endswith('.fastq.gz'):
            if 'short_1' in item:
                corrected_1 = os.path.join(corrected_dir, item)
            elif 'short_2' in item:
                corrected_2 = os.path.join(corrected_dir, item)
            elif 'unpaired' in item:
                corrected_u = os.path.join(corrected_dir, item)

    spades_start_time = time.time()
    try:
        print_with_timestamp('Running ' + run_name)
        if not os.path.exists(spades_dir):
            os.makedirs(spades_dir)
        spades_command = ['spades.py']
        if corrected_1 and corrected_2 and corrected_u:
            spades_command += ['-1', corrected_1,
                               '-2', corrected_2,
                               '-s', corrected_u,
                               '--only-assembler']
        else:
            spades_command += ['-1', short_1,
                               '-2', short_2]
        if 'PACBIO' in long.upper() or args.long_acc >= 80.0:
            spades_command += ['--pacbio', long]
        else:
            spades_command += ['--nanopore', long]
        spades_command += ['--careful',
                           '--threads', str(args.threads),
                           '-o', spades_dir]
        print_with_timestamp(' '.join(spades_command))
        try:
            spades_out = subprocess.check_output(spades_command, stderr=subprocess.STDOUT)
            with open(os.path.join(spades_dir, 'spades.out'), 'wb') as f:
                f.write(spades_out)
        except subprocess.CalledProcessError as e:
            print_with_timestamp('SPAdes encountered an error:\n' + e.output.decode())
            raise AssemblyError

    except AssemblyError:
        spades_assembly = None

    spades_time = time.time() - spades_start_time
    run_quast(spades_assembly, args, all_quast_results, simple_quast_results, 'SPAdes hybrid',
              long_count, long_depth, spades_time, spades_time, run_name, spades_dir)
    clean_up_spades_dir(spades_dir)
    return os.path.abspath(spades_dir)


def run_np_scarf(long_read_file, long_reads, long_count, long_depth, args, all_quast_results,
                 simple_quast_results, spades_dir, np_scarf_all_long_dir=None, extra_time=0.0):
    """
    Runs npScarf using the SPAdes assembly results:

    jsa.seq.sort -r -n --input scaffolds.fasta --output np_scarf.fasta
    bwa index np_scarf.fasta
    bwa mem -t 8 -k 11 -W 20 -r 10 -A 1 -B 1 -O 1 -E 1 -L 0 -a -Y np_scarf.fasta long_reads.fastq
        > alignments.sam
    jsa.np.gapcloser -b alignments.sam -seq np_scarf.fasta
    """
    run_name, np_scarf_dir = get_run_name_and_run_dir_name('npScarf', args.reference, long_depth,
                                                           args)
    np_scarf_assembly = os.path.join(np_scarf_dir, 'out.fin.fasta')

    np_scarf_start_time = time.time()
    try:
        print_with_timestamp('Running ' + run_name)
        if not os.path.exists(np_scarf_dir):
            os.makedirs(np_scarf_dir)

        np_scarf_fasta = os.path.join(np_scarf_dir, 'np_scarf.fasta')
        jsa_seq_sort_command = ['jsa.seq.sort',
                                '-r',
                                '-n',
                                '--input', os.path.join(spades_dir, 'contigs.fasta'),
                                '--output', np_scarf_fasta]
        print_with_timestamp(' '.join(jsa_seq_sort_command))
        try:
            subprocess.check_output(jsa_seq_sort_command, stderr=subprocess.STDOUT)
        except subprocess.CalledProcessError as e:
            print_with_timestamp('jsa.seq.sort encountered an error:\n' + e.output.decode())
            raise AssemblyError

        # If we are running npScarf with all long reads, then we run BWA Mem to get the SAM file.
        alignments_file = os.path.join(np_scarf_dir, 'alignments.sam')
        if not np_scarf_all_long_dir:
            bwa_index_command = ['bwa', 'index',
                                 np_scarf_fasta]
            print_with_timestamp(' '.join(bwa_index_command))
            try:
                subprocess.check_output(bwa_index_command, stderr=subprocess.STDOUT)
            except subprocess.CalledProcessError as e:
                print_with_timestamp('BWA index encountered an error:\n' + e.output.decode())
                raise AssemblyError

            alignments_sam = open(alignments_file, 'w')
            dev_null = open(os.devnull, 'w')
            bwa_mem_command = ['bwa', 'mem',
                               '-t', str(args.threads),
                               '-k', '11',
                               '-W', '20',
                               '-r', '10',
                               '-A', '1', '-B', '1', '-O', '1', '-E', '1', '-L', '0',
                               '-a', '-Y',
                               np_scarf_fasta,
                               long_read_file]
            print_with_timestamp(' '.join(bwa_mem_command))
            bwa_mem = subprocess.Popen(bwa_mem_command, stdout=alignments_sam, stderr=dev_null)
            try:
                bwa_mem.communicate()
                bwa_mem.wait()
            except subprocess.CalledProcessError:
                print_with_timestamp('BWA Mem encountered an error')
                raise AssemblyError
            dev_null.close()
            alignments_sam.close()

        # If we are running npScarf with a subset, there's no need to do the aligning again. Instead
        # we can just make a SAM by sampling from the previously done full set.
        else:
            create_subsampled_sam(os.path.join(np_scarf_all_long_dir, 'alignments.sam'),
                                  alignments_file, long_reads)

        jsa_np_gapcloser_command = ['jsa.np.gapcloser',
                                    '-b', alignments_file,
                                    '-seq', np_scarf_fasta]
        print_with_timestamp(' '.join(jsa_np_gapcloser_command))
        try:
            np_scarf_out = subprocess.check_output(jsa_np_gapcloser_command,
                                                   stderr=subprocess.STDOUT)
            with open(os.path.join(np_scarf_dir, 'np_scarf.out'), 'wb') as f:
                f.write(np_scarf_out)

        except subprocess.CalledProcessError as e:
            print_with_timestamp('jsa.np.gapcloser encountered an error:\n' + e.output.decode())
            raise AssemblyError

        shutil.move('out.fin.fasta', np_scarf_assembly)
        os.remove('out.fin.japsa')

    except AssemblyError:
        np_scarf_assembly = None

    np_scarf_time = time.time() - np_scarf_start_time
    full_time = np_scarf_time + extra_time
    run_quast(np_scarf_assembly, args, all_quast_results, simple_quast_results, 'npScarf',
              long_count, long_depth, np_scarf_time, full_time, run_name, np_scarf_dir)

    clean_up_np_scarf_dir(np_scarf_dir)
    return os.path.abspath(np_scarf_dir)


def run_cerulean(long_read_file, long_reads, long_count, long_depth, args, all_quast_results,
                 simple_quast_results, abyss_dir, cerulean_all_long_dir=None, extra_time=0.0):
    run_name, cerulean_dir = get_run_name_and_run_dir_name('Cerulean', args.reference, long_depth,
                                                           args)
    starting_path = os.path.abspath('.')
    cerulean_assembly = os.path.join(cerulean_dir, 'cerulean_final.fasta')

    print_with_timestamp('Running ' + run_name)
    if not os.path.exists(cerulean_dir):
        os.makedirs(cerulean_dir)

    cerulean_start_time = time.time()
    shutil.copy(os.path.join(abyss_dir, 'run-contigs.fa'),
                os.path.join(cerulean_dir, 'run-contigs.fa'))
    shutil.copy(os.path.join(abyss_dir, 'run-contigs.dot'),
                os.path.join(cerulean_dir, 'run-contigs.dot'))

    shutil.copyfile(long_read_file, os.path.join(cerulean_dir, 'long_reads.fastq.gz'))
    os.chdir(cerulean_dir)
    try:
        subprocess.check_output(['gunzip', 'long_reads.fastq.gz'], stderr=subprocess.STDOUT)
    except subprocess.CalledProcessError as e:
        quit_with_error('gzip encountered an error:\n' + e.output.decode())

    try:
        try:
            jelly_path = os.environ['JELLYPATH']
            fakequals_path = os.path.join(jelly_path, 'bin', 'fakeQuals.py')
            jelly_py_path = os.path.join(jelly_path, 'bin', 'Jelly.py')
        except KeyError:
            print_with_timestamp('JELLYPATH not set')
            raise AssemblyError
        try:
            cerulean_path = os.path.join(os.environ['CERULEANPATH'], 'Cerulean.py')
        except KeyError:
            print_with_timestamp('CERULEANPATH not set')
            raise AssemblyError
        try:
            sawriter_command = ['sawriter', 'run-contigs.fa']
            print_with_timestamp(' '.join(sawriter_command))
            subprocess.check_output(sawriter_command, stderr=subprocess.STDOUT)
        except subprocess.CalledProcessError as e:
            print_with_timestamp('sawriter encountered an error:\n' + e.output.decode())
            raise AssemblyError
        fastq_to_fasta('long_reads.fastq', 'long_reads.fasta')
        os.remove('long_reads.fastq')
        if cerulean_all_long_dir:
            create_subsampled_m4(os.path.join(cerulean_all_long_dir,
                                              'run_pacbio_contigs_mapping.fasta.m4'),
                                 os.path.join(cerulean_dir, 'run_pacbio_contigs_mapping.fasta.m4'),
                                 long_reads)
        else:
            try:
                blasr_command = ['blasr',
                                 'long_reads.fasta',
                                 'run-contigs.fa',
                                 '-minMatch', '10',
                                 '-minPctIdentity', '70',
                                 '-bestn', '30',
                                 '-nCandidates', '30',
                                 '-maxScore', '500',
                                 '-nproc', str(args.threads),
                                 '-noSplitSubreads',
                                 '-out', 'long_reads_contigs_mapping.fasta.m4']
                print_with_timestamp(' '.join(blasr_command))
                subprocess.check_output(blasr_command, stderr=subprocess.STDOUT, timeout=3600)
            except subprocess.CalledProcessError as e:
                print_with_timestamp('blasr encountered an error:\n' + e.output.decode())
                raise AssemblyError
            os.rename('long_reads_contigs_mapping.fasta.m4', 'run_pacbio_contigs_mapping.fasta.m4')
        try:
            cerulean_command = ['python', cerulean_path,
                                '--dataname', 'run',
                                '--basedir', '.',
                                '--nproc', str(args.threads)]
            print_with_timestamp(' '.join(cerulean_command))
            subprocess.check_output(cerulean_command, stderr=subprocess.STDOUT)
        except subprocess.CalledProcessError as e:
            print_with_timestamp('Cerulean.py encountered an error:\n' + e.output.decode())
            raise AssemblyError
        try:
            fakequals_command = ['python', fakequals_path, 'run_cerulean.fasta',
                                 'run_cerulean.qual']
            print_with_timestamp(' '.join(fakequals_command))
            subprocess.check_output(fakequals_command, stderr=subprocess.STDOUT)
        except subprocess.CalledProcessError as e:
            print_with_timestamp('fakeQuals.py encountered an error:\n' + e.output.decode())
            raise AssemblyError
        try:
            fakequals_command = ['python', fakequals_path, 'long_reads.fasta', 'long_reads.qual']
            print_with_timestamp(' '.join(fakequals_command))
            subprocess.check_output(fakequals_command, stderr=subprocess.STDOUT)
        except subprocess.CalledProcessError as e:
            print_with_timestamp('fakeQuals.py encountered an error:\n' + e.output.decode())
            raise AssemblyError
        if not os.path.exists('PBJelly'):
            os.makedirs('PBJelly')
        this_dir = os.path.abspath('.')
        template_protocol = open(os.path.join(jelly_path, 'docs', 'jellyExample', 'Protocol.xml'),
                                 'rt')
        this_protocol = open('Protocol.xml', 'wt')
        for line in template_protocol:
            line = line.replace('<reference>/__PATH__/_TO_/jellyExample/data/reference/'
                                'lambda.fasta',
                                '<reference>' + os.path.join(this_dir, 'run_cerulean.fasta'))
            line = line.replace('<outputDir>/__PATH__/_TO_/jellyExample/',
                                '<outputDir>' + os.path.join(this_dir, 'PBJelly') + '/')
            line = line.replace('-nproc 4',
                                '-nproc ' + str(args.threads))
            line = line.replace('baseDir="/__PATH__/_TO_/jellyExample/data/reads/"',
                                'baseDir="' + this_dir + '/"')
            line = line.replace('filtered_subreads.fastq', 'long_reads.fasta')
            this_protocol.write(line)
        this_protocol.close()
        template_protocol.close()
        try:
            setup_command = ['python', jelly_py_path, 'setup', 'Protocol.xml']
            mapping_command = ['python', jelly_py_path, 'mapping', 'Protocol.xml']
            support_command = ['python', jelly_py_path, 'support', 'Protocol.xml']
            extraction_command = ['python', jelly_py_path, 'extraction', 'Protocol.xml']
            assembly_command = ['python', jelly_py_path, 'assembly', 'Protocol.xml']
            output_command = ['python', jelly_py_path, 'output', 'Protocol.xml']
            print_with_timestamp(' '.join(setup_command))
            subprocess.check_output(setup_command, stderr=subprocess.STDOUT)
            print_with_timestamp(' '.join(mapping_command))
            subprocess.check_output(mapping_command, stderr=subprocess.STDOUT)
            print_with_timestamp(' '.join(support_command))
            subprocess.check_output(support_command, stderr=subprocess.STDOUT)
            print_with_timestamp(' '.join(extraction_command))
            subprocess.check_output(extraction_command, stderr=subprocess.STDOUT)
            print_with_timestamp(' '.join(assembly_command))
            subprocess.check_output(assembly_command, stderr=subprocess.STDOUT)
            print_with_timestamp(' '.join(output_command))
            subprocess.check_output(output_command, stderr=subprocess.STDOUT)
        except subprocess.CalledProcessError as e:
            print_with_timestamp('Jelly.py encountered an error:\n' + e.output.decode())
            raise AssemblyError
        shutil.copy(os.path.join(cerulean_dir, 'PBJelly', 'jelly.out.fasta'),
                    cerulean_assembly)

    except AssemblyError:
        cerulean_assembly = None
    else:
        clean_up_cerulean_dir(cerulean_dir)

    cerulean_time = time.time() - cerulean_start_time
    full_time = cerulean_time + extra_time
    os.chdir(starting_path)
    run_quast(cerulean_assembly, args, all_quast_results, simple_quast_results, 'Cerulean',
              long_count, long_depth, cerulean_time, full_time, run_name, cerulean_dir)

    return os.path.abspath(cerulean_dir)


def run_unicycler(args, short_1, short_2, long_read_filename, long_read_count,
                  long_read_depth, all_quast_results, simple_quast_results, bridging_mode,
                  unicycler_no_long_dir=None, unicycler_all_long_dir=None, replacement_time=None,
                  extra_time=0.0):
    run_name, unicycler_dir = get_run_name_and_run_dir_name('Unicycler_' + bridging_mode,
                                                            args.reference, long_read_depth, args)
    unicycler_assembly = os.path.join(unicycler_dir, 'assembly.fasta')

    if not os.path.exists(unicycler_dir):
        os.makedirs(unicycler_dir)

    # If an unbridged graph is already available, copy it over to save time (we won't have to run
    # SPAdes).
    if unicycler_no_long_dir:
        shutil.copyfile(os.path.join(unicycler_no_long_dir, '001_unbridged_graph.gfa'),
                        os.path.join(unicycler_dir, '001_unbridged_graph.gfa'))

    # If read alignments are already available, copy them over instead of doing the alignment again.
    if unicycler_all_long_dir:
        read_align_dir = os.path.join(unicycler_dir, 'read_alignment_temp')
        if not os.path.exists(read_align_dir):
            os.makedirs(read_align_dir)
            shutil.copyfile(os.path.join(unicycler_all_long_dir, 'read_alignment_temp',
                                         'long_read_alignments.sam'),
                            os.path.join(read_align_dir, 'long_read_alignments.sam'))

    unicycler_start_time = time.time()
    print_with_timestamp('Running ' + run_name)
    if not os.path.exists(unicycler_dir):
        os.makedirs(unicycler_dir)
    unicycler_command = ['unicycler',
                         '--short1', short_1,
                         '--short2', short_2]
    if long_read_filename:
        unicycler_command += ['--long', long_read_filename]
    else:
        unicycler_command += ['--no_long']
    unicycler_command += ['--out', unicycler_dir,
                          '--mode', bridging_mode,
                          '--keep_temp', '1',
                          '--threads', str(args.threads),
                          '--verbosity', '2',
                          '--no_rotate']
    print_with_timestamp(' '.join(unicycler_command))
    try:
        unicycler_out = subprocess.check_output(unicycler_command, stderr=subprocess.STDOUT)
        with open(os.path.join(unicycler_dir, 'unicycler.out'), 'wb') as f:
            f.write(unicycler_out)
    except subprocess.CalledProcessError as e:
        quit_with_error('Unicycler encountered an error:\n' + e.output.decode())

    unicycler_time = time.time() - unicycler_start_time
    if replacement_time:
        full_time = replacement_time
    else:
        full_time = unicycler_time
    full_time += extra_time

    run_quast(unicycler_assembly, args, all_quast_results, simple_quast_results,
              'Unicycler (' + bridging_mode + ')', long_read_count, long_read_depth, unicycler_time,
              full_time, run_name, unicycler_dir)

    clean_up_unicycler_dir(unicycler_dir)
    return os.path.abspath(unicycler_dir), full_time


# def run_canu(long, long_count, long_depth, args, all_quast_results, simple_quast_results):
#     """
#     Runs Canu.
#     """
#     run_name, canu_dir = get_run_name_and_run_dir_name('Canu', args.reference, long_depth, args)
#     canu_assembly = os.path.join(canu_dir, 'canu_out.fasta')
#
#     canu_start_time = time.time()
#     try:
#         print_with_timestamp('Running ' + run_name)
#         if not os.path.exists(canu_dir):
#             os.makedirs(canu_dir)
#         canu_command = ['canu',
#                         '-p', 'canu_out',
#                         '-d', canu_dir,
#                         '-nanopore-raw', long]
#         print_with_timestamp(' '.join(canu_command))
#         try:
#             canu_out = subprocess.check_output(canu_command, stderr=subprocess.STDOUT)
#             with open(os.path.join(canu_dir, 'canu.out'), 'wb') as f:
#                 f.write(canu_out)
#         except subprocess.CalledProcessError as e:
#             print_with_timestamp('Canu encountered an error:\n' + e.output.decode())
#             raise AssemblyError
#
#     except AssemblyError:
#         spades_assembly = None
#
#     canu_time = time.time() - canu_start_time
#     run_quast(spades_assembly, args, all_quast_results, simple_quast_results, 'Canu',
#               long_count, long_depth, canu_time, run_name, canu_dir)
#     clean_up_spades_dir(canu_dir)
#     return os.path.abspath(canu_dir)


def create_quast_results_tables(args):
    quast_results_filename = 'quast_results.tsv'
    simple_quast_results_filename = 'quast_results_simple.tsv'

    if not os.path.isfile(simple_quast_results_filename):
        simple_quast_results = open(simple_quast_results_filename, 'w')
        simple_quast_results.write("Reference name\t"
                                   "Assembler\t")
        if args.real_long:
            simple_quast_results.write("Full long read filename\t")
        else:
            simple_quast_results.write("Synthetic long read accuracy (%)\t"
                                       "Synthetic long read mean length (bp)\t")
        simple_quast_results.write("Long read depth\t"
                                   "Run time (seconds)\t"
                                   "Full run time (seconds)\t"
                                   "Reference pieces\t"
                                   "# contigs\t"
                                   "NGA50\t"
                                   "Completeness (%)\t"
                                   "Complete\t"
                                   "# misassemblies\t"
                                   "# local misassemblies\t"
                                   "# mismatches and indels per 100 kbp\n")
        simple_quast_results.close()

    if not os.path.isfile(quast_results_filename):
        quast_results = open(quast_results_filename, 'w')
        quast_results.write("Reference name\t"
                            "Reference size (bp)\t"
                            "Reference pieces\t"
                            "Assembler\t")
        if args.real_long:
            quast_results.write("Full long read filename\t")
        else:
            quast_results.write("Synthetic long read accuracy (%)\t"
                                "Synthetic long read mean length (bp)\t")
        quast_results.write("Long read depth\t"
                            "Long read count\t"
                            "Run time (seconds)\t"
                            "# contigs (>= 0 bp)\t"
                            "# contigs (>= 1000 bp)\t"
                            "# contigs (>= 5000 bp)\t"
                            "# contigs (>= 10000 bp)\t"
                            "# contigs (>= 25000 bp)\t"
                            "# contigs (>= 50000 bp)\t"
                            "Total length (>= 0 bp)\t"
                            "Total length (>= 1000 bp)\t"
                            "Total length (>= 5000 bp)\t"
                            "Total length (>= 10000 bp)\t"
                            "Total length (>= 25000 bp)\t"
                            "Total length (>= 50000 bp)\t"
                            "# contigs\t"
                            "Largest contig\t"
                            "Total length\t"
                            "Reference length\t"
                            "GC (%)\t"
                            "Reference GC (%)\t"
                            "N50\t"
                            "NG50\t"
                            "N75\t"
                            "NG75\t"
                            "L50\t"
                            "LG50\t"
                            "L75\t"
                            "LG75\t"
                            "# misassemblies\t"
                            "# misassembled contigs\t"
                            "Misassembled contigs length\t"
                            "# local misassemblies\t"
                            "# unaligned contigs\t"
                            "Unaligned length\t"
                            "Genome fraction (%)\t"
                            "Duplication ratio\t"
                            "# N's per 100 kbp\t"
                            "# mismatches per 100 kbp\t"
                            "# indels per 100 kbp\t"
                            "Largest alignment\t"
                            "Reference lengths\t"
                            "Largest alignment per reference\t"
                            "Largest alignment per reference sum\t"
                            "Percent complete per reference\t"
                            "Percent complete total\t"
                            "NA50\t"
                            "NGA50\t"
                            "NA75\t"
                            "NGA75\t"
                            "LA50\t"
                            "LGA50\t"
                            "LA75\t"
                            "LGA75\n")
        quast_results.close()
    return os.path.abspath(quast_results_filename), os.path.abspath(simple_quast_results_filename)


def run_quast(assembly, args, all_quast_results, simple_quast_results, assembler_name,
              long_read_count, long_read_depth, run_time, full_time, run_name, run_dir_name):
    if all_quast_results is None or simple_quast_results is None:
        return

    reference_name = get_reference_name_from_filename(args.reference)
    long_read_acc = float_to_str(args.long_acc, 1) if long_read_count > 0 else ''
    long_read_len = str(args.long_len) if long_read_count > 0 else ''

    ref_length, ref_count = get_fasta_length_and_seq_count(args.reference)
    quast_line = [reference_name, str(ref_length), str(ref_count), assembler_name]
    if args.real_long:
        quast_line.append(args.real_long)
    else:
        quast_line += [long_read_acc, long_read_len]
    quast_line += [float_to_str(long_read_depth, 5), str(long_read_count), str(run_time)]

    simple_quast_line = [reference_name, assembler_name]
    if args.real_long:
        simple_quast_line.append(args.real_long)
    else:
        simple_quast_line += [long_read_acc, long_read_len]
    simple_quast_line += [float_to_str(long_read_depth, 5), str(run_time), str(full_time)]

    if assembly is None:
        print_with_timestamp('Skipping QUAST for ' + run_name)
        quast_line += [''] * 51
        simple_quast_line += [''] * 3

    else:
        print_with_timestamp('Running QUAST for ' + run_name)
        quast_dir = os.path.join(os.path.dirname(os.path.normpath(run_dir_name)), 'quast_results',
                                 os.path.basename(os.path.normpath(run_dir_name)))
        this_quast_results = os.path.join(quast_dir, 'transposed_report.tsv')

        quast_command = ['quast.py',
                         assembly,
                         '-R', args.reference,
                         '-o', quast_dir,
                         '-l', '"' + run_name.replace(',', '') + '"',
                         '--threads', str(args.threads)]
        print_with_timestamp(' '.join(quast_command))
        try:
            subprocess.check_output(quast_command, stderr=subprocess.STDOUT)
        except subprocess.CalledProcessError as e:
            print_with_timestamp('QUAST encountered an error:\n' + e.output.decode())
            quast_line += [''] * 51
            simple_quast_line += [''] * 3
        else:
            with open(this_quast_results, 'rt') as results:
                headers = results.readline().strip().split('\t')
                results = results.readline().strip().split('\t')
                quast_line += results[1:]

                try:
                    contig_count = int(get_quast_result(headers, results, '# contigs (>= 0 bp)'))
                except ValueError:
                    contig_count = 0
                try:
                    percent_complete = float(get_quast_result(headers, results,
                                                              'Percent complete total'))
                except ValueError:
                    percent_complete = 0.0
                try:
                    misassemblies = int(get_quast_result(headers, results, '# misassemblies'))
                except ValueError:
                    misassemblies = 0

                complete_per_ref = get_quast_result(headers, results,
                                                    'Percent complete per reference')
                if not complete_per_ref:
                    complete = 'N'
                else:
                    percents = [float(x.split(':')[1]) for x in complete_per_ref.split(' ')]
                    if all(x >= 99.0 for x in percents) and ref_count == contig_count and \
                            misassemblies == 0:
                        complete = 'Y'
                    else:
                        complete = 'N'

                simple_quast_line.append(str(ref_count))
                simple_quast_line.append(str(contig_count))
                simple_quast_line.append(get_quast_result(headers, results, 'NGA50'))
                simple_quast_line.append(str(percent_complete))
                simple_quast_line.append(complete)
                simple_quast_line.append(str(misassemblies))
                simple_quast_line.append(get_quast_result(headers, results,
                                                          '# local misassemblies'))
                try:
                    mismatches = float(get_quast_result(headers, results,
                                                        '# mismatches per 100 kbp'))
                    indels = float(get_quast_result(headers, results, '# indels per 100 kbp'))
                    simple_quast_line.append(str(mismatches + indels))
                except ValueError:
                    simple_quast_line.append('')

    all_results_line = '\t'.join(quast_line) + '\n'
    with open(all_quast_results, 'at') as all_results:
        all_results.write(all_results_line)

    simple_results_line = '\t'.join(simple_quast_line) + '\n'
    with open(simple_quast_results, 'at') as simple_results:
        simple_results.write(simple_results_line)


def get_quast_result(headers, results, column_name):
    try:
        return results[headers.index(column_name)]
    except ValueError:
        return ''


def load_fasta(filename):
    """
    Returns a list of tuples (header, seq) for each record in the fasta file.
    """
    fasta_seqs = []
    fasta_file = open(filename, 'rt')
    name = ''
    sequence = ''
    for line in fasta_file:
        line = line.strip()
        if not line:
            continue
        if line[0] == '>':  # Header line = start of new contig
            if name:
                name_parts = name.split()
                seq_name = name_parts[0]
                relative_depth = name_parts[1]
                if len(name_parts) > 2 and name_parts[2] == 'linear':
                    circular = False
                else:
                    circular = True
                fasta_seqs.append((seq_name, sequence, relative_depth, circular))
                sequence = ''
            name = line[1:]
        else:
            sequence += line
    if name:
        name_parts = name.split()
        seq_name = name_parts[0]
        relative_depth = name_parts[1]
        if len(name_parts) > 2 and name_parts[2] == 'linear':
            circular = False
        else:
            circular = True
        fasta_seqs.append((seq_name, sequence, relative_depth, circular))
    fasta_file.close()
    return fasta_seqs


def get_fasta_length_and_seq_count(filename):
    records = load_fasta(filename)
    seq_count = len(records)
    length = sum(len(x[1]) for x in records)
    return length, seq_count


def quit_with_error(message):
    print_with_timestamp('Error: ' + message)
    sys.exit(1)


def float_to_str(num, decimals):
    format_str = '%.' + str(decimals) + 'f'
    return format_str % num


def get_reference_name_from_filename(reference):
    return reference.split('/')[-1].split('.')[0]


def get_relative_depths(reference):
    references = load_fasta(reference)
    longest_len = 0
    longest_depth = 0.0
    relative_depths = []
    for ref in references:
        try:
            depth = float(ref[2])
        except ValueError:
            quit_with_error('Reference sequences must indicate depth in each header line')
        relative_depths.append(depth)
        length = len(ref[1])
        if length > longest_len:
            longest_len = length
            longest_depth = depth
    return [x / longest_depth for x in relative_depths]


def get_run_name_and_run_dir_name(assembler_name, reference, long_read_depth, args):
    ref_name = get_reference_name_from_filename(reference)
    read_accuracy = '{0:g}'.format(args.long_acc)
    read_length = '{0:g}'.format(args.long_len)
    read_depth = '{:0>5.2f}'.format(long_read_depth)
    run_name = ref_name + ', ' + assembler_name
    run_dir_name = assembler_name
    if long_read_depth > 0.0 and not args.real_long:
        run_name += ', ' + read_accuracy + '%, ' + read_length + ' bp'
        run_dir_name += '_' + read_accuracy + '_' + read_length
    run_name += ', ' + read_depth + 'x'
    run_dir_name += '_' + read_depth + 'x'
    run_dir_name = run_dir_name.replace(',', '_').replace(' ', '_')

    full_run_name = run_name
    full_run_dir_name = os.path.abspath(run_dir_name)
    index = 1
    while os.path.isdir(full_run_dir_name):
        index += 1
        full_run_name = run_name + '_' + str(index)
        full_run_dir_name = os.path.abspath(run_dir_name) + '_' + str(index)

    return full_run_name, full_run_dir_name


def clean_up_spades_dir(spades_dir):
    """
    Deletes everything in spades_dir except for assembly_graph.fastg, scaffolds.paths,
    scaffolds.fasta and spades.out.
    """
    for item in os.listdir(spades_dir):
        path = os.path.join(spades_dir, item)
        if 'assembly_graph.fastg' in path or 'before_rr.fasta' in path or \
                'scaffolds.paths' in path or 'scaffolds.fasta' in path or \
                'contigs.paths' in path or 'contigs.fasta' in path or \
                'spades.out' in path or 'corrected' in path:
            continue
        if os.path.isfile(path):
            os.remove(path)
        elif os.path.isdir(path):
            shutil.rmtree(path)


def clean_up_cerulean_dir(cerulean_dir):
    """
    Deletes everything in cerulean_dir except for cerulean_final.fasta and m4 files.
    """
    for item in os.listdir(cerulean_dir):
        path = os.path.join(cerulean_dir, item)
        if 'cerulean_final.fasta' in path or '.m4' in path:
            continue
        if os.path.isfile(path):
            os.remove(path)
        elif os.path.isdir(path):
            shutil.rmtree(path)


def clean_up_unicycler_dir(unicycler_dir):
    """
    Deletes everything in unicycler_dir except for assembly.gfa, assembly.fasta,
    read_alignment_temp, 001_unbridged_graph.gfa and unicycler.out.
    """
    for item in os.listdir(unicycler_dir):
        path = os.path.join(unicycler_dir, item)
        if 'gfa' in path or 'assembly.fasta' in path or 'read_alignment_temp' in path or \
                '001_unbridged_graph.gfa' in path or 'unicycler.out' in path:
            continue
        if os.path.isfile(path):
            os.remove(path)
        elif os.path.isdir(path):
            shutil.rmtree(path)
    alignment_dir = os.path.join(unicycler_dir, 'read_alignment_temp')
    if os.path.isdir(alignment_dir):
        for item in os.listdir(alignment_dir):
            path = os.path.join(alignment_dir, item)
            if 'long_read_alignments.sam' in path:
                continue
            if os.path.isfile(path):
                os.remove(path)
            elif os.path.isdir(path):
                shutil.rmtree(path)


def clean_up_np_scarf_dir(np_scarf_dir):
    """
    Deletes everything in np_scarf_dir except for out.fin.fasta, alignments.sam and np_scarf.out.
    """
    for item in os.listdir(np_scarf_dir):
        path = os.path.join(np_scarf_dir, item)
        if 'out.fin.fasta' in path or 'alignments.sam' in path or 'np_scarf.out' in path:
            continue
        if os.path.isfile(path):
            os.remove(path)
        elif os.path.isdir(path):
            shutil.rmtree(path)


def check_for_existing_assembly(assembly_path):
    """
    Checks to see if the assembly exists, and if so, returns its directory.
    """
    if os.path.isfile(assembly_path):
        assembly_path_parts = assembly_path.split('/')
        assembly_dir = '/'.join(assembly_path_parts[:-1])
        return assembly_dir
    else:
        return ''


def get_next_available_set_number(starting_path, num):
    while True:
        dir_name = str(num).rjust(6, '0')
        full_dir_path = os.path.join(starting_path, dir_name)
        if not os.path.exists(full_dir_path):
            os.makedirs(full_dir_path)
            return full_dir_path, num
        num += 1


def fastq_to_fasta(fastq_filename, fasta_filename):
    with open(fastq_filename, 'rt') as fastq:
        with open(fasta_filename, 'wt') as fasta:
            for line in fastq:
                fasta.write('>' + line[1:])
                fasta.write(next(fastq))
                next(fastq)
                next(fastq)


# def run_command_through_time(command, time_path):
#     command_with_time = [time_path, '-v'] + command
#     output = subprocess.check_output(command_with_time, stderr=subprocess.STDOUT)
#     output_lines = output.splitlines()
#     max_mem = 0
#     for line in output_lines:
#         if b'Maximum resident set size' in line:
#             max_mem = int(line.split(b': ')[1])
#     return output, max_mem
#
#
# def run_shell_command_through_time(command, time_path):
#     command_with_time = time_path + ' -v ' + command
#     output = subprocess.check_output(command_with_time, stderr=subprocess.STDOUT, shell=True)
#     output_lines = output.splitlines()
#     max_mem = 0
#     for line in output_lines:
#         if b'Maximum resident set size' in line:
#             max_mem = int(line.split(b': ')[1])
#     return output, max_mem


def set_up_env_var():
    """
    Set up environment variables for Cerulean and PBJelly.
    """
    if os.path.isdir('/Users/Ryan/Applications/PBSuite_15.8.24'):
        local_pb_path = '/Users/Ryan/Applications/PBSuite_15.8.24'
    else:
        local_pb_path = '/vlsci/SG0006/rwick/PBSuite_15.8.24'
    if os.path.isdir('/Users/Ryan/Applications/Cerulean/src'):
        local_cerulean_path = '/Users/Ryan/Applications/Cerulean/src'
    else:
        local_cerulean_path = '/vlsci/SG0006/rwick/Cerulean/src'
    if 'SWEETPATH' not in os.environ:
        os.environ['SWEETPATH'] = local_pb_path
    if 'PYTHONPATH' not in os.environ:
        os.environ['PYTHONPATH'] = local_pb_path
    elif 'PBSuite' not in os.environ['PYTHONPATH']:
        os.environ['PYTHONPATH'] += ':' + local_pb_path
    if 'PBSuite' not in os.environ['PATH']:
        os.environ['PATH'] += ':' + os.environ['SWEETPATH'] + '/bin/'
    if 'JELLYPATH' not in os.environ:
        os.environ['JELLYPATH'] = local_pb_path
    if 'CERULEANPATH' not in os.environ:
        os.environ['CERULEANPATH'] = local_cerulean_path


def print_with_timestamp(string):
    print('{:%Y-%m-%d %H:%M:%S}'.format(datetime.datetime.now()) + '  ' + string, flush=True)


if __name__ == '__main__':
    main()
