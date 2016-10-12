#!/usr/bin/env python3

import argparse
import os
import sys
import subprocess
import collections
import textwrap


def main():
    args = get_arguments()
    check_files(args)
    set_env_vars(args.pb_path)
    clean_up()

    if args.bax:
        print_round_header('Converting bax.h5 reads to BAM format')
        make_reads_bam(args.bax)
        args.bam = 'subreads.bam'

    round_num = 1
    files = [f for f in os.listdir('.') if os.path.isfile(f)]
    while any(f.startswith('%03d' % round_num) for f in files):
        round_num += 1

    latest_assembly = args.fasta
    done = False
    while not done:
        print_round_header('Round ' + str(round_num))
        done, latest_assembly = polish_assembly(latest_assembly, round_num,
                                                args.min_align_length, args.bam, args.threads,
                                                args.min_ref_length, args.max_homopolymer)
        round_num += 1
    print_finished(latest_assembly)


def get_arguments():
    parser = argparse.ArgumentParser(description='PacBio polishing loop - runs arrow algorithm '
                                                 'until results no longer change')
    parser.add_argument('--bax', nargs='+', type=str,
                        help='PacBio raw bax.h5 read files')
    parser.add_argument('--bam', type=str,
                        help='PacBio BAM of subreads (either --bam or --bax is required)')
    parser.add_argument('--fasta', type=str, required=True,
                        help='Assembly to polish')
    parser.add_argument('--threads', type=int, required=True,
                        help='CPU threads to use in alignment and consensus')
    parser.add_argument('--pb_path', type=str, default='/home/UNIMELB/wickr-u/pb_tools',
                        help='Path to PacBio tools (should contain bin and lib directories)')
    parser.add_argument('--min_align_length', type=int, default=1000,
                        help='Minimum BLASR alignment length')
    parser.add_argument('--min_ref_length', type=int, default=10000,
                        help='Minimum reference size for PacBio polishing')
    parser.add_argument('--max_homopolymer', type=int, default=3,
                        help='Changes to the length of homopolymers larger than this will be '
                             'ignored')

    args = parser.parse_args()
    if not args.bax and not args.bam:
        parser.error('either --bax or --bam is required')
    if args.bax and args.bam:
        parser.error('only one of the following can be used: --bax or --bam')
    return args


def check_files(args):
    if args.bam and not os.path.isfile(args.bam + '.pbi'):
        sys.exit('Error: missing ' + args.bam + '.pbi')
    if not os.path.isdir(args.pb_path):
        sys.exit('Error: path does not exist: ' + args.pb_path)
    if not os.path.isdir(args.pb_path):
        sys.exit('Error: path does not exist: ' + args.pb_path)
    args.bax2bam = os.path.join(args.pb_path, 'bin', 'bax2bam')
    args.pbalign = os.path.join(args.pb_path, 'bin', 'pbalign')
    args.arrow = os.path.join(args.pb_path, 'bin', 'arrow')


def clean_up():
    files_to_delete = []
    for f in ['pbalign_alignments.bam', 'pbalign_alignments.bam.pbi', 'pbalign_alignments.bam.bai']:
        if os.path.isfile(f):
            files_to_delete.append(f)
    fai_files = [f for f in os.listdir('.') if os.path.isfile(f) and f.endswith('.fai')]
    for f in fai_files:
        files_to_delete.append(f)
    if files_to_delete:
        print_command(['rm'] + files_to_delete)
    for f in files_to_delete:
        os.remove(f)


def set_env_vars(pb_path):
    lib_path = os.path.join(pb_path, 'lib')
    bin_path = os.path.join(pb_path, 'bin')
    python_path = os.path.join(pb_path, 'lib', 'python2.7', 'dist-packages')
    if 'LD_LIBRARY_PATH' in os.environ:
        os.environ['LD_LIBRARY_PATH'] = lib_path + ':' + os.environ['LD_LIBRARY_PATH']
    else:
        os.environ['LD_LIBRARY_PATH'] = lib_path
    if 'PATH' in os.environ:
        os.environ['PATH'] = bin_path + ':' + os.environ['PATH']
    else:
        os.environ['PATH'] = bin_path
    if 'PYTHONPATH' in os.environ:
        os.environ['PYTHONPATH'] = python_path + ':' + os.environ['PYTHONPATH']
    else:
        os.environ['PYTHONPATH'] = python_path


def make_reads_bam(bax):
    command = ['bax2bam'] + bax
    print_command(command)
    run_command_print_output(command)
    files = [f for f in os.listdir('.') if os.path.isfile(f)]
    files_to_delete = []
    for filename in files:
        if 'scrap' in filename:
            os.remove(filename)
            files_to_delete.append(filename)
        elif 'subreads' in filename and filename.endswith('.bam'):
            os.rename(filename, 'subreads.bam')
            print_command(['mv', filename, 'subreads.bam'])
        elif 'subreads' in filename and filename.endswith('.bam.pbi'):
            os.rename(filename, 'subreads.bam.pbi')
            print_command(['mv', filename, 'subreads.bam.pbi'])
    print_command(['rm'] + files_to_delete)
    files = [f for f in os.listdir('.') if os.path.isfile(f)]
    if 'subreads.bam' not in files:
        sys.exit('Error: bax2bam failed to make subreads.bam')
    if 'subreads.bam.pbi' not in files:
        sys.exit('Error: bax2bam failed to make subreads.bam.pbi')


def polish_assembly(fasta, round_num, min_align_length, reads_bam, threads, min_ref_length,
                    max_homopolymer):
    align_reads(fasta, min_align_length, reads_bam, threads)
    raw_variants_gff = '%03d' % round_num + '_raw_variants.gff'
    raw_variants = load_variants(raw_variants_gff)
    genomic_consensus(fasta, threads, raw_variants)
    clean_up()
    filtered_variants = filter_variants(fasta, raw_variants_gff, min_ref_length, max_homopolymer)
    polished_fasta = '%03d' % round_num + '_polish.fasta'
    apply_variants(fasta, filtered_variants, polished_fasta)
    print_result(raw_variants, filtered_variants, polished_fasta)
    done = len(filtered_variants) == 0
    return done, polished_fasta


def align_reads(fasta, min_align_length, reads_bam, threads):
    command = ['pbalign',
               '--nproc', str(threads),
               '--minLength', str(min_align_length),
               '--algorithmOptions="--minRawSubreadScore 800 --bestn 1"',
               reads_bam,
               fasta,
               'pbalign_alignments.bam']
    print_command(command)
    run_command_print_output(command)
    files = [f for f in os.listdir('.') if os.path.isfile(f)]
    if 'pbalign_alignments.bam' not in files:
        sys.exit('Error: pbalign failed to make pbalign_alignments.bam')
    if 'pbalign_alignments.bam.pbi' not in files:
        sys.exit('Error: pbalign failed to make pbalign_alignments.bam.pbi')
    if 'pbalign_alignments.bam.bai' not in files:
        sys.exit('Error: pbalign failed to make pbalign_alignments.bam.bai')


def genomic_consensus(fasta, threads, raw_variants_filename):
    subprocess.call(['samtools', 'faidx', fasta])
    command = ['arrow',
               'pbalign_alignments.bam',
               '-j', str(threads),
               '--noEvidenceConsensusCall', 'reference',
               '-r', fasta,
               '-o', raw_variants_filename]
    print_command(command)
    run_command_print_output(command)
    files = [f for f in os.listdir('.') if os.path.isfile(f)]
    if raw_variants_filename not in files:
        sys.exit('Error: arrow failed to make ' + raw_variants_filename)


def filter_variants(fasta, raw_variants_gff, min_ref_length, max_homopolymer):
    """
    This function produces a new GFF file without variants in homopolymer runs.
    """
    reference = load_fasta(fasta)
    raw_variants = load_variants(raw_variants_gff)
    filtered_variants = [x for x in raw_variants if
                         not x.is_homopolymer_change(reference, max_homopolymer) and
                         not x.has_small_reference(reference, min_ref_length)]
    filtered_variants_gff = raw_variants_gff.replace('_raw_', '_filtered_')
    with open(filtered_variants_gff, 'wt') as new_gff:
        with open(raw_variants_gff, 'rt') as old_gff:
            for line in old_gff:
                if line.startswith('##'):
                    new_gff.write(line)
        for variant in filtered_variants:
            new_gff.write(variant.original_line + '\n')
    return filtered_variants


def apply_variants(in_fasta, variants, out_fasta):
    """
    This function creates a new FASTA file by applying the variants to an existing FASTA file.
    """
    in_seqs = load_fasta(in_fasta)
    out_seqs = collections.OrderedDict()
    for name, seq in in_seqs.items():
        seq_variants = sorted([x for x in variants if x.ref_name == name],
                              key=lambda x: x.start_pos)
        new_seq = ''
        pos = 0
        for variant in seq_variants:
            new_seq += seq[pos:variant.start_pos]
            new_seq += variant.variant_seq
            pos = variant.end_pos
        new_seq += seq[pos:]
        out_seqs[name] = new_seq
    with open(out_fasta, 'wt') as fasta:
        for name, seq in out_seqs.items():
            fasta.write('>' + name + '\n')
            fasta.write(add_line_breaks_to_sequence(seq, 60))


def print_command(command):
    print('\033[1m' + ' '.join(command) + '\033[0m', flush=True)


def print_round_header(text):
    print()
    print('\033[1m' + '\033[36m' + '\033[4m' + text + '\033[0m', flush=True)


def print_result(raw_variants, filtered_variants, latest_assembly):
    print_result_line('Unfiltered variants:        ' + str(len(raw_variants)))
    print_result_line('Filtered variants:          ' + str(len(filtered_variants)))
    filtered_variant_positions = ', '.join([str(x.start_pos) for x in filtered_variants])
    filtered_variant_position_lines = textwrap.wrap(filtered_variant_positions, 50)
    print_result_line('Filtered variant positions: ' + filtered_variant_position_lines[0])
    for line in filtered_variant_position_lines[1:]:
        print_result_line('                            ' + line)
    print_result_line('Polished FASTA:             ' + latest_assembly)


def print_result_line(text):
    print('\033[1m' + '\033[32m' + text + '\033[0m', flush=True)


def print_finished(latest_assembly):
    print()
    result = 'No changes in ' + latest_assembly + '. All done!'
    print('\033[1m' + '\033[32m' + result + '\033[0m', flush=True)


def print_warning(warning):
    print('\033[31m' + warning + '\033[0m', flush=True)


def run_command_print_output(command):
    process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    for line in iter(process.stdout.readline, b''):
        print(line.rstrip().decode(), flush=True)


def load_fasta(filename):
    """
    Returns a dictionary (header: seq) for each record in the fasta file.
    """
    fasta_seqs = collections.OrderedDict()
    fasta_file = open(filename, 'rt')
    name = ''
    sequence = ''
    for line in fasta_file:
        line = line.strip()
        if not line:
            continue
        if line[0] == '>':  # Header line = start of new contig
            if name:
                fasta_seqs[name.split()[0]] = sequence
                sequence = ''
            name = line[1:]
        else:
            sequence += line
    if name:
        fasta_seqs[name.split()[0]] = sequence
    fasta_file.close()
    return fasta_seqs


def load_variants(gff_file):
    variants = []
    with open(gff_file, 'rt') as gff:
        for line in gff:
            line = line.strip()
            if line and not line.startswith('##'):
                variants.append(Variant(line))
    return variants


def has_multiple_bases(seq):
    seq = seq.upper()
    base_counts = [seq.count('A'), seq.count('C'), seq.count('G'), seq.count('T')]
    base_counts = [x for x in base_counts if x > 0]
    return len(base_counts) > 1


def homopolymer_size(seq, pos):
    size = 1
    starting_base = seq[pos]
    forward_pos = pos + 1
    while True:
        if forward_pos >= len(seq) or seq[forward_pos] != starting_base:
            break
        size += 1
        forward_pos += 1
    reverse_pos = pos - 1
    while True:
        if reverse_pos < 0 or seq[reverse_pos] != starting_base:
            break
        size += 1
        reverse_pos -= 1
    return size


def add_line_breaks_to_sequence(sequence, line_length):
    seq_with_breaks = ''
    pos = 0
    while pos < len(sequence):
        seq_with_breaks += sequence[pos:pos+line_length] + '\n'
        pos += line_length
    return seq_with_breaks


class Variant(object):
    def __init__(self, gff_line):
        """
        https://github.com/PacificBiosciences/GenomicConsensus/blob/master/doc/VariantsGffSpecification.rst
        """
        line_parts = gff_line.split('\t')
        attributes = {x.split('=')[0]: x.split('=')[1] for x in line_parts[8].split(';')}

        self.original_line = gff_line
        self.ref_name = line_parts[0]
        self.type = line_parts[2]
        self.ref_seq = attributes['reference'].replace('.', '')
        self.variant_seq = attributes['variantSeq'].replace('.', '')
        self.start_pos = int(line_parts[3]) - 1
        self.end_pos = self.start_pos + len(self.ref_seq)
        self.confidence = int(attributes['confidence'])
        self.coverage = int(attributes['coverage'])

    def has_small_reference(self, reference, min_ref_length):
        return len(reference[self.ref_name]) < min_ref_length

    def is_homopolymer_change(self, reference, max_homopolymer):
        if self.type != 'insertion' and self.type != 'deletion':
            return False
        if has_multiple_bases(self.ref_seq) or has_multiple_bases(self.variant_seq):
            return False
        ref_seq = reference[self.ref_name]
        return homopolymer_size(ref_seq, self.start_pos) > max_homopolymer


if __name__ == '__main__':
    main()
