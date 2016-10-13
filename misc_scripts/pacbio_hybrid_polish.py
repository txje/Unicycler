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
        done, latest_assembly = polish_assembly(latest_assembly, round_num, args.short1,
                                                args.short2, args.min_align_length, args.bam,
                                                args.threads, args.min_ref_length,
                                                args.homopolymer, args.illumina_alt)
        round_num += 1
    print_finished(latest_assembly)


def get_arguments():
    parser = argparse.ArgumentParser(description='PacBio hybrid polishing',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--short1', type=str, required=True,
                        help='Short read FASTQ - first reads in pair')
    parser.add_argument('--short2', type=str, required=True,
                        help='Short read FASTQ - second reads in pair')
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
                        help='Contigs shorter than this will not be polished')
    parser.add_argument('--homopolymer', type=int, default=3,
                        help='Changes to a homopolymer of this length or greater will be ignored')
    parser.add_argument('--illumina_alt', type=float, default=0.05,
                        help='Fraction of Illumina reads which must vary from the reference in '
                             'order for a base to be eligible for PacBio polishing')

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
    for f in ['pbalign_alignments.bam', 'pbalign_alignments.bam.pbi',  'pbalign_alignments.bam.bai',
              'illumina_alignments.bam', 'illumina_alignments.bam.bai']:
        if os.path.isfile(f):
            files_to_delete.append(f)
    for f in [f for f in os.listdir('.') if os.path.isfile(f) and f.endswith('.fai')]:
        files_to_delete.append(f)
    for f in [f for f in os.listdir('.') if os.path.isfile(f) and f.endswith('.bt2')]:
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


def polish_assembly(fasta, round_num, short1, short2, min_align_length, reads_bam, threads,
                    min_ref_length, max_homopolymer, illumina_alt):
    # Conduct PacBio alignment and consensus
    align_pacbio_reads(fasta, min_align_length, reads_bam, threads)
    raw_variants_gff = '%03d' % round_num + '_raw_variants.gff'
    genomic_consensus(fasta, threads, raw_variants_gff)
    raw_variants = load_variants(raw_variants_gff, fasta, min_ref_length, max_homopolymer)
    clean_up()

    # Conduct Illumina alignment
    align_illumina_reads(fasta, short1, short2, threads)
    for variant in raw_variants:
        variant.assess_against_illumina_alignments(fasta)

    filtered_variants = filter_variants(fasta, raw_variants_gff, min_ref_length, max_homopolymer,
                                        illumina_alt)
    polished_fasta = '%03d' % round_num + '_polish.fasta'
    apply_variants(fasta, filtered_variants, polished_fasta)
    print_result(raw_variants, filtered_variants, polished_fasta)
    done = len(filtered_variants) == 0
    return done, polished_fasta


def align_pacbio_reads(fasta, min_align_length, reads_bam, threads):
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


def align_illumina_reads(fasta, short1, short2, threads):
    bowtie2build_command = ['bowtie2-build', fasta, fasta]
    print_command(bowtie2build_command)
    subprocess.check_output(bowtie2build_command, stderr=subprocess.STDOUT)

    bowtie2_command = ['bowtie2',
                       '--end-to-end',
                       '--very-sensitive',
                       '--threads', str(threads),
                       '--no-unal',
                       '-I', '0', '-X', '2000',
                       '-x', fasta,
                       '-1', short1, '-2', short2]
    samtools_view_command = ['samtools', 'view', '-@', str(threads), '-Sb', '-']
    samtools_sort_command = ['samtools', 'sort', '-@', str(threads), '-',
                             '-o', 'illumina_alignments.bam']
    bowtie2 = subprocess.Popen(bowtie2_command, stdout=subprocess.PIPE)
    samtools_view = subprocess.Popen(samtools_view_command,
                                     stdin=bowtie2.stdout, stdout=subprocess.PIPE)
    samtools_sort = subprocess.Popen(samtools_sort_command,
                                     stdin=samtools_view.stdout, stdout=subprocess.PIPE)
    print_command(bowtie2_command + ['|'] + samtools_view_command + ['|'] + samtools_sort_command)
    bowtie2.stdout.close()
    samtools_view.stdout.close()
    samtools_sort.communicate()

    samtools_index_command = ['samtools', 'index', 'illumina_alignments.bam']
    print_command(samtools_index_command)
    run_command_print_output(samtools_index_command)


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


def filter_variants(fasta, raw_variants_gff, min_ref_length, max_homopolymer, illumina_alt):
    """
    This function produces a new GFF file without variants:
      * in homopolymer runs
      * in small references
      * TO DO: in regions where the Illumina mapping is super-solid
    """
    reference = load_fasta(fasta)
    raw_variants = load_variants(raw_variants_gff, fasta, max_homopolymer, min_ref_length)

    print(get_out_header())
    filtered_variants = []
    for variant in raw_variants:
        variant_output = variant.get_out_line() + '\t'
        if not variant.is_homopolymer_change(reference, max_homopolymer) and \
                not variant.has_small_reference(reference, min_ref_length) and \
                variant.illumina_alt_fraction >= illumina_alt:
            filtered_variants.append(variant)
            variant_output += '\033[32m' + 'PASS' + '\033[0m'
        else:
            variant_output += '\033[31m' + 'FAIL' + '\033[0m'
        print(variant_output, flush=True)

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
    if filtered_variants:
        filtered_variant_position_lines = textwrap.wrap(filtered_variant_positions, 70)
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
    try:
        process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        for line in iter(process.stdout.readline, b''):
            print(line.rstrip().decode(), flush=True)
    except subprocess.CalledProcessError as e:
        sys.exit(e.output.decode())


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


def load_variants(gff_file, fasta, min_ref_length, max_homopolymer):
    reference = load_fasta(fasta)
    variants = []
    with open(gff_file, 'rt') as gff:
        for line in gff:
            line = line.strip()
            if line and not line.startswith('##'):
                variants.append(PacBioVariant(line, reference, min_ref_length, max_homopolymer))
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


class PacBioVariant(object):
    def __init__(self, gff_line, reference, min_ref_length, max_homopolymer):
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

        ref_seq = reference[self.ref_name]
        self.has_small_reference = len(ref_seq) < min_ref_length

        if self.type != 'insertion' and self.type != 'deletion':
            self.is_homopolymer_change = False
        elif has_multiple_bases(self.ref_seq) or has_multiple_bases(self.variant_seq):
            self.is_homopolymer_change = False
        else:
            self.is_homopolymer_change = homopolymer_size(ref_seq, self.start_pos) > max_homopolymer

        self.illumina_alt_fraction = 0.0

    def assess_against_illumina_alignments(self, reference_fasta):
        ref_location = self.ref_name + ':' + str(self.start_pos - 5) + '-' + str(self.end_pos + 5)
        freebayes_command = ['freebayes',
                             '-f', reference_fasta,
                             '-p', '1',
                             '-r', ref_location,
                             '--report-monomorphic',
                             '--min-alternate-fraction', '0',
                             '--pooled-continuous',
                             '--min-alternate-count', '1',
                             '--haplotype-length', '0',
                             'illumina_alignments.bam']
        freebayes_out = subprocess.check_output(freebayes_command, stderr=subprocess.STDOUT)
        freebayes_lines = [x for x in freebayes_out.decode().split('\n')
                           if x and not x.startswith('#')]
        for line in freebayes_lines:
            line_parts = line.split('\t')
            start_pos = int(line_parts[1]) - 1
            end_pos = start_pos + len(line_parts[3]) - 1  # inclusive end
            if start_pos <= self.start_pos <= end_pos:
                ref_occurrences = int(line.split(';RO=')[1].split(';')[0])
                if ';AO=' in line:
                    alt_occurrences = int(line.split(';AO=')[1].split(';')[0])
                else:
                    alt_occurrences = 0
                alt_fraction = alt_occurrences / (ref_occurrences + alt_occurrences)
                self.illumina_alt_fraction = max(self.illumina_alt_fraction, alt_fraction)

    def get_out_line(self):
        homopolymer = 'yes' if self.is_homopolymer_change else 'no'
        small = 'yes' if self.has_small_reference else 'no'
        return '\t'.join([self.ref_name,
                          str(self.start_pos + 1),
                          self.ref_seq,
                          self.variant_seq,
                          homopolymer,
                          small,
                          "%.2f" % self.illumina_alt_fraction])


def get_out_header():
    return '\033[1m' + '\t'.join(['CONTIG',
                                  'POSITION',
                                  'REF',
                                  'ALT',
                                  'HOMO',
                                  'SMALL',
                                  'ILLUMINA',
                                  'RESULT']) + '\033[0m'


if __name__ == '__main__':
    main()
