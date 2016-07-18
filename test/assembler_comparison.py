#!/usr/bin/env python3
"""
Tool for comparing the assemblies made by Unicycler, SPAdes, hybridSPAdes and npScarf.

Usage:
assembler_comparison.py --reference path/to/reference.fasta

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


def main():
    args = get_args()
    starting_path = os.path.abspath('.')

    scaled_ref_length = get_scaled_ref_length(args)
    ref_name = get_reference_name_from_filename(args.reference)
    short_1, short_2 = make_fake_short_reads(args, starting_path, ref_name)
    quast_results, simple_quast_results = create_quast_results_tables()

    # Run SPAdes and Unicycler on short read data alone.
    short_read_only_dir = os.path.abspath('000000_only_short_reads')
    if not os.path.exists(short_read_only_dir):
        os.makedirs(short_read_only_dir)
    os.chdir(short_read_only_dir)
    abyss_dir = run_abyss(short_1, short_2, args, quast_results, simple_quast_results)
    spades_no_long_dir = run_regular_spades(short_1, short_2, args, quast_results,
                                            simple_quast_results)
    unicycler_no_long_dir = run_unicycler_no_long(short_1, short_2, args, quast_results,
                                                  simple_quast_results)
    # Run Cerulean with short reads only?
    os.chdir(starting_path)

    # The program runs indefinitely, always running more tests until the user kills it.
    dir_num = 1
    while True:
        dir_name, dir_num = get_next_available_set_number(starting_path, dir_num)
        new_path = os.path.join(starting_path, dir_name)
        os.chdir(new_path)

        long_filename, long_reads = make_fake_long_reads(args)
        long_read_count = len(long_reads)

        # Run Unicycler on the full set of long reads - will be the source of alignments for
        # subsampled Unicycler runs.
        long_read_depth = get_long_read_depth(long_reads, scaled_ref_length)
        unicycler_all_long_dir = run_unicycler_all_long(short_1, short_2, long_filename, args,
                                                        long_read_count, long_read_depth,
                                                        quast_results, simple_quast_results,
                                                        unicycler_no_long_dir)

        # Run hybridSPAdes on the full set of long reads.
        run_hybrid_spades(short_1, short_2, long_filename, long_read_count, long_read_depth, args,
                          quast_results, simple_quast_results, spades_no_long_dir)

        # Run npScarf on the full set of long reads - will be the source of alignments for
        # subsampled npScarf runs.
        np_scarf_all_long_dir = run_np_scarf(long_filename, long_reads, long_read_count,
                                             long_read_depth, args, quast_results,
                                             simple_quast_results, spades_no_long_dir)

        # Randomly subsample this read set 10 different times.
        for _ in range(10):

            # Choose a random number of reads to include in the subsample.
            subsampled_count = random.randint(1, long_read_count - 1)

            # Subsample the long reads.
            subsampled_reads = random.sample(long_reads, subsampled_count)
            subsampled_long_read_depth = get_long_read_depth(subsampled_reads, scaled_ref_length)
            subsampled_filename = ref_name + '_long_subsampled_' + str(subsampled_count) + '.fastq'
            print('\nSubsampling to', subsampled_count, 'reads')
            save_long_reads_to_fastq(subsampled_reads, subsampled_filename)
            try:
                subprocess.check_output(['gzip', subsampled_filename], stderr=subprocess.STDOUT)
            except subprocess.CalledProcessError as e:
                quit_with_error('gzip encountered an error:\n' + e.output.decode())
            subsampled_filename += '.gz'

            # Run Unicycler on the subsampled long reads. This should be relatively fast, because we
            # can skip the short read assembly and the alignment. The polishing step still takes a
            # while, though.
            run_unicycler_subsampled_long(short_1, short_2, subsampled_reads, subsampled_filename,
                                          args, subsampled_count, subsampled_long_read_depth,
                                          quast_results, simple_quast_results,
                                          unicycler_all_long_dir)

            # Run hybridSPAdes on the subsampled long reads.
            run_hybrid_spades(short_1, short_2, subsampled_filename, subsampled_count,
                              subsampled_long_read_depth, args, quast_results, simple_quast_results,
                              spades_no_long_dir)

            # Run npScarf on the subsampled long reads. This is very fast because we can skip the
            # alignment step.
            run_np_scarf(subsampled_filename, subsampled_reads, subsampled_count,
                         subsampled_long_read_depth, args, quast_results, simple_quast_results,
                         spades_no_long_dir, np_scarf_all_long_dir=np_scarf_all_long_dir)

            # Run Cerulean?

        os.chdir(starting_path)


def get_args():
    """
    Specifies the command line arguments required by the script.
    """
    parser = argparse.ArgumentParser(description='Assembly tester')
    parser.add_argument('--reference', type=str, required=True,
                        help='The reference genome to shred and reassemble')
    parser.add_argument('--short_depth', type=float, default=50.0,
                        help='Base read depth for fake short reads')
    parser.add_argument('--long_depth', type=float, default=20.0,
                        help='Base read depth for fake long reads')
    parser.add_argument('--long_acc', type=float, default=80.0,
                        help='Mean accuracy for long reads')
    parser.add_argument('--long_len', type=int, default=10000,
                        help='Mean length for long reads')
    parser.add_argument('--model_qc', type=str, default='model_qc_clr',
                        help='Model QC file for pbsim')
    parser.add_argument('--rotation_count', type=int, default=100, required=False,
                        help='The number of times to run read simulators with random start '
                             'positions')
    parser.add_argument('--threads', type=int, required=False, default=8,
                        help='Number of CPU threads')

    args = parser.parse_args()
    args.reference = os.path.abspath(args.reference)
    args.model_qc = os.path.abspath(args.model_qc)
    return args


class AssemblyError(Exception):
    pass


def make_fake_short_reads(args, current_path, ref_name):
    """
    Runs ART to generate fake Illumina reads. Runs ART separate for each sequence in the reference
    file (to control relative depth) and at multiple sequence rotations (to ensure circular
    assembly).
    """
    read_filename_1 = os.path.abspath(ref_name + '_short_1.fastq')
    read_filename_2 = os.path.abspath(ref_name + '_short_2.fastq')

    if os.path.isfile(read_filename_1 + '.gz') and os.path.isfile(read_filename_2 + '.gz'):
        read_filename_1 += '.gz'
        read_filename_2 += '.gz'
        print('\nReads already exist:', read_filename_1, read_filename_2)
        return read_filename_1, read_filename_2

    print('\nGenerating synthetic short reads', flush=True)

    temp_dir = os.path.join(current_path, 'temp_' + ref_name)
    if not os.path.exists(temp_dir):
        os.makedirs(temp_dir)
    os.chdir(temp_dir)

    references = load_fasta(args.reference)
    relative_depths = get_relative_depths(args.reference)

    # This will hold all simulated short reads. Each read is a list of 8 strings: the first four are
    # for the first read in the pair, the second four are for the second.
    short_read_pairs = []

    read_prefix = 1  # Used to prevent duplicate read names.
    for i, ref in enumerate(references):

        short_depth = relative_depths[i] * args.short_depth
        short_depth_per_rotation = short_depth / args.rotation_count

        ref_seq = ref[1]
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

    os.chdir(current_path)

    random.shuffle(short_read_pairs)
    reads_1 = open(read_filename_1, 'w')
    reads_2 = open(read_filename_2, 'w')
    for read_pair in short_read_pairs:
        reads_1.write(read_pair[0] + '\n')
        reads_1.write(read_pair[1] + '\n')
        reads_1.write(read_pair[2] + '\n')
        reads_1.write(read_pair[3] + '\n')
        reads_2.write(read_pair[4] + '\n')
        reads_2.write(read_pair[5] + '\n')
        reads_2.write(read_pair[6] + '\n')
        reads_2.write(read_pair[7] + '\n')
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
    print(read_filename_1)
    print(read_filename_2)

    shutil.rmtree(temp_dir)
    return read_filename_1, read_filename_2


def make_fake_long_reads(args):
    """
    Runs pbsim to generate fake long reads. Runs pbsim separate for each sequence in the reference
    file (to control relative depth) and at multiple sequence rotations (to ensure circular
    assembly).
    """
    ref_name = get_reference_name_from_filename(args.reference)
    long_filename = os.path.abspath(ref_name + '_long.fastq')
    print('\nGenerating synthetic long reads', flush=True)

    references = load_fasta(args.reference)
    relative_depths = get_relative_depths(args.reference)

    # This will hold all simulated long reads. Each read is a list of 4 strings: one for each
    # line of the FASTQ.
    long_reads = []

    read_prefix = 1  # Used to prevent duplicate read names.
    for i, ref in enumerate(references):

        long_depth = relative_depths[i] * args.long_depth
        long_depth_per_rotation = long_depth / args.rotation_count

        ref_seq = ref[1]
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

            long_reads += run_pbsim(temp_fasta_filename, long_depth_per_rotation, args,
                                    str(read_prefix))

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

    print(long_filename)
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
                    read_name = line.split('\t')[0]
                    if read_name in subsampled_read_names:
                        subsampled_sam.write(line)


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


def run_pbsim(input_fasta, depth, args, read_prefix):
    pbsim_command = ['pbsim',
                     '--depth', str(depth),
                     '--length-min', '100',
                     '--length-max', '100000',
                     '--accuracy-min', '0.5',
                     '--accuracy-max', '1.0',
                     '--model_qc', args.model_qc,
                     '--length-mean', str(args.long_len),
                     '--length-sd', str(args.long_len / 10),
                     '--accuracy-mean', str(args.long_acc / 100.0),
                     '--accuracy-sd', str(0.02),
                     input_fasta]
    try:
        subprocess.check_output(pbsim_command, stderr=subprocess.STDOUT)
    except subprocess.CalledProcessError as e:
        quit_with_error('pbsim encountered an error:\n' + e.output.decode())

    reads = []
    with open('sd_0001.fastq', 'rt') as fastq:
        for line in fastq:
            name = '@' + read_prefix + '_' + line.strip().split('_')[1]
            sequence = next(fastq).strip()
            _ = next(fastq)
            qualities = next(fastq).strip()
            reads.append((name, sequence, '+', qualities))

    os.remove('sd_0001.fastq')
    os.remove('sd_0001.maf')
    os.remove('sd_0001.ref')

    return reads


def run_abyss(short_1, short_2, args, all_quast_results, simple_quast_results):
    """
    Runs ABySS on the short reads to compare with SPAdes.
    """
    run_name, abyss_dir = get_run_name_and_run_dir_name('ABySS', args.reference, 0.0, 0)
    abyss_contigs = os.path.join(abyss_dir, 'run-contigs.fa')
    abyss_scaffolds = os.path.join(abyss_dir, 'run-scaffolds.fa')
    existing_dir = check_for_existing_assembly(abyss_scaffolds)
    if existing_dir:
        return os.path.abspath(existing_dir)
    starting_path = os.path.abspath('.')

    abyss_start_time = time.time()
    try:
        print('\nRunning', run_name, flush=True)
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
        print(abyss_command)
        try:
            abyss_out = subprocess.check_output(abyss_command, stderr=subprocess.STDOUT, shell=True)
            with open(os.path.join(abyss_dir, 'abyss.out'), 'wb') as f:
                f.write(abyss_out)
        except subprocess.CalledProcessError as e:
            print('ABySS encountered an error:\n' + e.output.decode())
            raise AssemblyError

    except AssemblyError:
        abyss_scaffolds = None
        abyss_contigs = None

    abyss_time = time.time() - abyss_start_time
    os.chdir(starting_path)

    run_quast(abyss_scaffolds, args, all_quast_results, simple_quast_results, 'ABySS scaffolds',
              0, 0.0, abyss_time)
    run_quast(abyss_contigs, args, all_quast_results, simple_quast_results, 'ABySS scaffolds',
              0, 0.0, abyss_time)

    return os.path.abspath(abyss_dir)


def run_regular_spades(short_1, short_2, args, all_quast_results, simple_quast_results):
    """
    Runs SPAdes with only short reads (i.e. not hybridSPAdes).
    """
    run_name, spades_dir = get_run_name_and_run_dir_name('SPAdes', args.reference, 0.0, 0)
    spades_scaffolds = os.path.join(spades_dir, 'scaffolds.fasta')
    spades_contigs = os.path.join(spades_dir, 'contigs.fasta')
    spades_before_rr = os.path.join(spades_dir, 'before_rr.fasta')
    existing_dir = check_for_existing_assembly(spades_scaffolds)
    if existing_dir:
        return os.path.abspath(existing_dir)

    spades_start_time = time.time()
    try:
        print('\nRunning', run_name, flush=True)
        if not os.path.exists(spades_dir):
            os.makedirs(spades_dir)
        spades_command = ['spades.py',
                          '-1', short_1,
                          '-2', short_2,
                          '--careful',
                          '--threads', str(args.threads),
                          '-o', spades_dir]
        print(' '.join(spades_command))
        try:
            spades_out = subprocess.check_output(spades_command, stderr=subprocess.STDOUT)
            with open(os.path.join(spades_dir, 'spades.out'), 'wb') as f:
                f.write(spades_out)
        except subprocess.CalledProcessError as e:
            print('SPAdes encountered an error:\n' + e.output.decode())
            raise AssemblyError

    except AssemblyError:
        spades_scaffolds = None
        spades_contigs = None
        spades_before_rr = None

    spades_time = time.time() - spades_start_time

    run_quast(spades_scaffolds, args, all_quast_results, simple_quast_results, 'SPAdes scaffolds',
              0, 0.0, spades_time)
    run_quast(spades_contigs, args, all_quast_results, simple_quast_results, 'SPAdes contigs',
              0, 0.0, spades_time)
    run_quast(spades_before_rr, args, all_quast_results, simple_quast_results, 'SPAdes before_rr',
              0, 0.0, spades_time)

    clean_up_spades_dir(spades_dir)
    return os.path.abspath(spades_dir)


def run_hybrid_spades(short_1, short_2, long, long_count, long_depth, args, all_quast_results,
                      simple_quast_results, spades_no_long_dir):
    """
    Runs hybridSPAdes using the --nanopore option.
    """
    run_name, spades_dir = get_run_name_and_run_dir_name('SPAdes', args.reference, long_depth,
                                                         long_count)
    spades_assembly = os.path.join(spades_dir, 'scaffolds.fasta')
    existing_dir = check_for_existing_assembly(spades_assembly)
    if existing_dir:
        return os.path.abspath(existing_dir)

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
        print('\nRunning', run_name, flush=True)
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
        spades_command += ['--nanopore', long,
                           '--careful',
                           '--threads', str(args.threads),
                           '-o', spades_dir]
        print(' '.join(spades_command))
        try:
            spades_out = subprocess.check_output(spades_command, stderr=subprocess.STDOUT)
            with open(os.path.join(spades_dir, 'spades.out'), 'wb') as f:
                f.write(spades_out)
        except subprocess.CalledProcessError as e:
            print('SPAdes encountered an error:\n' + e.output.decode())
            raise AssemblyError

    except AssemblyError:
        spades_assembly = None

    spades_time = time.time() - spades_start_time
    run_quast(spades_assembly, args, all_quast_results, simple_quast_results, 'SPAdes hybrid',
              long_count, long_depth, spades_time)
    clean_up_spades_dir(spades_dir)
    return os.path.abspath(spades_dir)


def run_np_scarf(long_read_file, long_reads, long_count, long_depth, args, all_quast_results,
                 simple_quast_results, spades_dir, np_scarf_all_long_dir=None):
    """
    Runs npScarf using the SPAdes assembly results:

    jsa.seq.sort -r -n --input scaffolds.fasta --output np_scarf.fasta
    bwa index np_scarf.fasta
    bwa mem -t 8 -k 11 -W 20 -r 10 -A 1 -B 1 -O 1 -E 1 -L 0 -a -Y np_scarf.fasta long_reads.fastq
        > alignments.sam
    jsa.np.gapcloser -b alignments.sam -seq np_scarf.fasta
    """
    run_name, np_scarf_dir = get_run_name_and_run_dir_name('npScarf', args.reference, long_depth,
                                                           long_count)
    np_scarf_assembly = os.path.join(np_scarf_dir, 'out.fin.fasta')
    existing_dir = check_for_existing_assembly(np_scarf_assembly)
    if existing_dir:
        return os.path.abspath(existing_dir)

    np_scarf_start_time = time.time()

    try:
        print('\nRunning', run_name, flush=True)
        if not os.path.exists(np_scarf_dir):
            os.makedirs(np_scarf_dir)

        np_scarf_fasta = os.path.join(np_scarf_dir, 'np_scarf.fasta')
        jsa_seq_sort_command = ['jsa.seq.sort',
                                '-r',
                                '-n',
                                '--input', os.path.join(spades_dir, 'contigs.fasta'),
                                '--output', np_scarf_fasta]
        print(' '.join(jsa_seq_sort_command))
        try:
            subprocess.check_output(jsa_seq_sort_command, stderr=subprocess.STDOUT)
        except subprocess.CalledProcessError as e:
            print('jsa.seq.sort encountered an error:\n' + e.output.decode(), flush=True)
            raise AssemblyError

        # If we are running npScarf with all long reads, then we run BWA Mem to get the SAM file.
        alignments_file = os.path.join(np_scarf_dir, 'alignments.sam')
        if not np_scarf_all_long_dir:
            bwa_index_command = ['bwa', 'index',
                                 np_scarf_fasta]
            print(' '.join(bwa_index_command))
            try:
                subprocess.check_output(bwa_index_command, stderr=subprocess.STDOUT)
            except subprocess.CalledProcessError as e:
                print('BWA index encountered an error:\n' + e.output.decode())
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
            print(' '.join(bwa_mem_command))
            bwa_mem = subprocess.Popen(bwa_mem_command, stdout=alignments_sam, stderr=dev_null)
            try:
                bwa_mem.communicate()
                bwa_mem.wait()
            except subprocess.CalledProcessError:
                print('BWA Mem encountered an error')
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
        print(' '.join(jsa_np_gapcloser_command))
        try:
            np_scarf_out = subprocess.check_output(jsa_np_gapcloser_command,
                                                   stderr=subprocess.STDOUT)
            with open(os.path.join(np_scarf_dir, 'np_scarf.out'), 'wb') as f:
                f.write(np_scarf_out)

        except subprocess.CalledProcessError as e:
            print('jsa.np.gapcloser encountered an error:\n' + e.output.decode())
            raise AssemblyError

        shutil.move('out.fin.fasta', np_scarf_assembly)
        os.remove('out.fin.japsa')

    except AssemblyError:
        np_scarf_assembly = None

    np_scarf_time = time.time() - np_scarf_start_time
    run_quast(np_scarf_assembly, args, all_quast_results, simple_quast_results, 'npScarf',
              long_count, long_depth, np_scarf_time)

    clean_up_np_scarf_dir(np_scarf_dir)
    return os.path.abspath(np_scarf_dir)


def run_unicycler_no_long(short_1, short_2, args, all_quast_results, simple_quast_results):
    run_name, unicycler_dir = get_run_name_and_run_dir_name('Unicycler', args.reference, 0.0, 0)
    unicycler_assembly = os.path.join(unicycler_dir, 'assembly.fasta')
    existing_dir = check_for_existing_assembly(unicycler_assembly)
    if existing_dir:
        return os.path.abspath(existing_dir)

    unicycler_start_time = time.time()
    print('\nRunning', run_name, flush=True)
    if not os.path.exists(unicycler_dir):
        os.makedirs(unicycler_dir)
    unicycler_command = ['unicycler',
                         '--short1', short_1,
                         '--short2', short_2,
                         '--no_long',
                         '--out', unicycler_dir,
                         '--keep_temp', '0',
                         '--threads', str(args.threads),
                         '--verbosity', '2']
    print(' '.join(unicycler_command))
    try:
        unicycler_out = subprocess.check_output(unicycler_command, stderr=subprocess.STDOUT)
        with open(os.path.join(unicycler_dir, 'unicycler.out'), 'wb') as f:
            f.write(unicycler_out)
    except subprocess.CalledProcessError as e:
        quit_with_error('Unicycler encountered an error:\n' + e.output.decode())
    unicycler_time = time.time() - unicycler_start_time
    run_quast(unicycler_assembly, args, all_quast_results, simple_quast_results, 'Unicycler',
              0, 0.0, unicycler_time)

    clean_up_unicycler_dir(unicycler_dir)
    return os.path.abspath(unicycler_dir)


def run_unicycler_all_long(short_1, short_2, long, args, long_read_count, long_read_depth,
                           all_quast_results, simple_quast_results, unicycler_no_long_dir):
    run_name, unicycler_dir = get_run_name_and_run_dir_name('Unicycler', args.reference,
                                                            long_read_depth, long_read_count)
    unicycler_assembly = os.path.join(unicycler_dir, 'assembly.fasta')
    existing_dir = check_for_existing_assembly(unicycler_assembly)
    if existing_dir:
        return os.path.abspath(existing_dir)

    if not os.path.exists(unicycler_dir):
        os.makedirs(unicycler_dir)

    # Copy over the unbridged graph, to save time.
    shutil.copyfile(os.path.join(unicycler_no_long_dir, '001_unbridged_graph.gfa'),
                    os.path.join(unicycler_dir, '001_unbridged_graph.gfa'))

    unicycler_start_time = time.time()
    print('\nRunning', run_name, flush=True)
    unicycler_command = ['unicycler',
                         '--short1', short_1,
                         '--short2', short_2,
                         '--long', long,
                         '--out', unicycler_dir,
                         '--keep_temp', '1',
                         '--threads', str(args.threads),
                         '--verbosity', '2']
    print(' '.join(unicycler_command))
    try:
        unicycler_out = subprocess.check_output(unicycler_command, stderr=subprocess.STDOUT)
        with open(os.path.join(unicycler_dir, 'unicycler.out'), 'wb') as f:
            f.write(unicycler_out)
    except subprocess.CalledProcessError as e:
        quit_with_error('Unicycler encountered an error:\n' + e.output.decode())
    unicycler_time = time.time() - unicycler_start_time
    run_quast(unicycler_assembly, args, all_quast_results, simple_quast_results, 'Unicycler',
              long_read_count, long_read_depth, unicycler_time)

    clean_up_unicycler_dir(unicycler_dir)
    return os.path.abspath(unicycler_dir)


def run_unicycler_subsampled_long(short_1, short_2, subsampled_reads, subsampled_filename,
                                  args, subsampled_count, subsampled_depth,
                                  all_quast_results, simple_quast_results, unicycler_all_long_dir):
    run_name, unicycler_dir = get_run_name_and_run_dir_name('Unicycler', args.reference,
                                                            subsampled_depth, subsampled_count)
    unicycler_assembly = os.path.join(unicycler_dir, 'assembly.fasta')
    existing_dir = check_for_existing_assembly(unicycler_assembly)
    if existing_dir:
        return os.path.abspath(existing_dir)

    if not os.path.exists(unicycler_dir):
        os.makedirs(unicycler_dir)

    # Copy over the unbridged graph, to save time.
    shutil.copyfile(os.path.join(unicycler_all_long_dir, '001_unbridged_graph.gfa'),
                    os.path.join(unicycler_dir, '001_unbridged_graph.gfa'))

    # Instead of making Unicycler realign the subsampled reads, we just copy the appropriate
    # alignments over from the full read set Unicycler run.
    read_align_dir = os.path.join(unicycler_dir, 'read_alignment_temp')
    if not os.path.exists(read_align_dir):
        os.makedirs(read_align_dir)
    create_subsampled_sam(os.path.join(unicycler_all_long_dir, 'read_alignment_temp',
                                       'long_read_alignments.sam'),
                          os.path.join(read_align_dir, 'long_read_alignments.sam'),
                          subsampled_reads)

    unicycler_start_time = time.time()
    print('\nRunning', run_name, flush=True)
    if not os.path.exists(unicycler_dir):
        os.makedirs(unicycler_dir)
    unicycler_command = ['unicycler',
                         '--short1', short_1,
                         '--short2', short_2,
                         '--long', subsampled_filename,
                         '--out', unicycler_dir,
                         '--keep_temp', '0',
                         '--threads', str(args.threads),
                         '--verbosity', '2']
    print(' '.join(unicycler_command))
    try:
        unicycler_out = subprocess.check_output(unicycler_command, stderr=subprocess.STDOUT)
        with open(os.path.join(unicycler_dir, 'unicycler.out'), 'wb') as f:
            f.write(unicycler_out)
    except subprocess.CalledProcessError as e:
        quit_with_error('Unicycler encountered an error:\n' + e.output.decode())
    unicycler_time = time.time() - unicycler_start_time
    run_quast(unicycler_assembly, args, all_quast_results, simple_quast_results, 'Unicycler',
              subsampled_count, subsampled_depth, unicycler_time)

    clean_up_unicycler_dir(unicycler_dir)
    return os.path.abspath(unicycler_dir)


def create_quast_results_tables():
    quast_results_filename = 'quast_results.tsv'
    simple_quast_results_filename = 'quast_results_simple.tsv'

    if not os.path.isfile(simple_quast_results_filename):
        simple_quast_results = open(simple_quast_results_filename, 'w')
        simple_quast_results.write("Reference name\t"
                                   "Assembler\t"
                                   "Long read depth\t"
                                   "Long read accuracy (%)\t"
                                   "Long read mean size (bp)\t"
                                   "Completeness (%)\t"
                                   "# misassemblies\t"
                                   "# mismatches and indels per 100 kbp\n")
        simple_quast_results.close()

    if not os.path.isfile(quast_results_filename):
        quast_results = open(quast_results_filename, 'w')
        quast_results.write("Reference name\t"
                            "Reference size (bp)\t"
                            "Reference pieces\t"
                            "Assembler\t"
                            "Long read count\t"
                            "Long read depth\t"
                            "Long read accuracy (%)\t"
                            "Long read mean size (bp)\t"
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
              long_read_count, long_read_depth, run_time):
    reference_name = get_reference_name_from_filename(args.reference)
    run_name, run_dir_name = get_run_name_and_run_dir_name(assembler_name, args.reference,
                                                           long_read_depth, long_read_count)

    ref_length, ref_count = get_fasta_length_and_seq_count(args.reference)
    quast_line = [reference_name, str(ref_length), str(ref_count), assembler_name,
                  str(long_read_count), float_to_str(long_read_depth, 5),
                  float_to_str(args.long_acc, 1), str(args.long_len),
                  str(run_time)]
    simple_quast_line = [reference_name, assembler_name, float_to_str(long_read_depth, 5),
                         float_to_str(args.long_acc, 1), str(args.long_len)]

    if assembly is None:
        print('\nSkipping QUAST for', run_name, flush=True)
        quast_line += [''] * 51
        simple_quast_line += [''] * 3

    else:
        print('\nRunning QUAST for', run_name, flush=True)
        quast_dir = os.path.join('quast_results', run_dir_name)
        this_quast_results = os.path.join(quast_dir, 'transposed_report.tsv')

        quast_command = ['quast.py',
                         assembly,
                         '-R', args.reference,
                         '-o', quast_dir,
                         '-l', '"' + run_name.replace(',', '') + '"',
                         '--threads', str(args.threads)]
        print(' '.join(quast_command))
        try:
            subprocess.check_output(quast_command, stderr=subprocess.STDOUT)
        except subprocess.CalledProcessError as e:
            print('QUAST encountered an error:\n' + e.output.decode(), flush=True)
            quast_line += [''] * 51
            simple_quast_line += [''] * 3
        else:
            with open(this_quast_results, 'rt') as results:
                results.readline()  # header line
                result_parts = results.readline().split('\t')
                quast_line += result_parts[1:]
                simple_quast_line.append(result_parts[43])
                simple_quast_line.append(result_parts[27])
                simple_quast_line.append(str(float(result_parts[36]) + float(result_parts[37])) +
                                         '\n')

    with open(all_quast_results, 'at') as all_results:
        all_results.write('\t'.join(quast_line))

    with open(simple_quast_results, 'at') as simple_results:
        simple_results.write('\t'.join(simple_quast_line))


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
                fasta_seqs.append((name.split()[0], sequence, name.split()[-1]))
                sequence = ''
            name = line[1:]
        else:
            sequence += line
    if name:
        fasta_seqs.append((name.split()[0], sequence, name.split()[-1]))
    fasta_file.close()
    return fasta_seqs


def get_fasta_length_and_seq_count(filename):
    records = load_fasta(filename)
    seq_count = len(records)
    length = sum(len(x[1]) for x in records)
    return length, seq_count


def quit_with_error(message):
    print('\nError:', message, file=sys.stderr, flush=True)
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


def get_run_name_and_run_dir_name(assembler_name, reference, long_read_depth, long_read_count):
    ref_name = get_reference_name_from_filename(reference)
    run_name = ref_name + ', ' + assembler_name + ', ' + str(long_read_count) + ' long reads, ' + \
        float_to_str(long_read_depth, 2) + 'x'
    run_dir_name = run_name.replace(', ', '_').replace(' ', '_').replace('.', '_')
    return run_name, run_dir_name


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


def clean_up_unicycler_dir(unicycler_dir):
    """
    Deletes everything in unicycler_dir except for assembly.gfa, assembly.fasta,
    read_alignment_temp, 001_unbridged_graph.gfa and unicycler.out.
    """
    for item in os.listdir(unicycler_dir):
        path = os.path.join(unicycler_dir, item)
        if 'assembly.gfa' in path or 'assembly.fasta' in path or 'read_alignment_temp' in path or \
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
    Checks to see if the assembly exists, and if so, returns its directory. The tricky part is
    that it doesn't care if the depth part matches exactly, but the assembler name, sample name
    and read count do have to match exactly.
    """
    assembly_path_parts = assembly_path.split('/')
    assembly_name = assembly_path_parts[-1]
    assembly_dir = '/'.join(assembly_path_parts[:-1])
    assembly_dir_first_part = assembly_dir.split('_long_reads')[0]
    for item in os.listdir('.'):
        if os.path.isdir(item) and item.startswith(assembly_dir_first_part):
            existing_assembly = os.path.join(item, assembly_name)
            if os.path.isfile(existing_assembly):
                print('\nAssembly already exists:', existing_assembly)
                return item
            else:
                shutil.rmtree(item)
                return ''
    return ''


def get_next_available_set_number(starting_path, num):
    while True:
        dir_name = str(num).rjust(6, '0')
        if not os.path.exists(os.path.join(starting_path, dir_name)):
            os.makedirs(dir_name)
            return dir_name, num
        num += 1


if __name__ == '__main__':
    main()
