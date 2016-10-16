#!/usr/bin/env python3
"""
This script is for polishing a PacBio genome that was finished using Unicycler.

Author: Ryan Wick
email: rrwick@gmail.com
"""

import argparse
import os
import shutil
import sys
import subprocess
import collections
import datetime
import statistics
import re
from .misc import add_line_breaks_to_sequence, load_fasta, MyHelpFormatter
from .alignment import get_ref_shift_from_cigar_part


def main():
    args = get_arguments()
    get_tool_paths(args)
    clean_up()

    if args.bax:
        print_round_header('Converting bax.h5 reads to BAM format')
        make_reads_bam(args)
        args.bam = 'subreads.bam'

    round_num = 1
    files = [f for f in os.listdir('.') if os.path.isfile(f)]
    while any(f.startswith('%03d' % round_num) for f in files):
        round_num += 1

    latest_assembly = args.fasta
    done = False
    while not done:
        print_round_header('Round ' + str(round_num))
        done, latest_assembly = polish_assembly(latest_assembly, round_num, args)
        clean_up()
        round_num += 1
    print_finished(latest_assembly)


def get_arguments():
    parser = argparse.ArgumentParser(description='PacBio hybrid polishing',
                                     formatter_class=MyHelpFormatter)
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
    parser.add_argument('--min_align_length', type=int, default=1000,
                        help='Minimum BLASR alignment length')
    parser.add_argument('--homopolymer', type=int, default=4,
                        help='Changes to a homopolymer of this length or greater will be ignored')
    parser.add_argument('--large', type=int, default=10,
                        help='Indels of this size or greater will be assess as large indels')
    parser.add_argument('--illumina_alt', type=float, default=10.0,
                        help='Percent of Illumina reads which must vary from the reference in '
                             'order for a base to be eligible for PacBio polishing')
    parser.add_argument('--pitchfork', type=str, default='',
                        help='Path to Pitchfork installation of PacBio tools (should contain bin '
                             'and lib directories)')
    parser.add_argument('--samtools', type=str, default=None,
                        help='path to samtools executable')
    parser.add_argument('--bowtie2', type=str, default=None,
                        help='path to bowtie2 executable')
    parser.add_argument('--bax2bam', type=str, default=None,
                        help='path to bax2bam executable')
    parser.add_argument('--pbalign', type=str, default=None,
                        help='path to pbalign executable')
    parser.add_argument('--arrow', type=str, default=None,
                        help='path to arrow executable')
    parser.add_argument('--freebayes', type=str, default=None,
                        help='path to freebayes executable')

    args = parser.parse_args()
    if not args.bax and not args.bam:
        parser.error('either --bax or --bam is required')
    if args.bax and args.bam:
        parser.error('only one of the following can be used: --bax or --bam')
    if args.bam and not os.path.isfile(args.bam + '.pbi'):
        sys.exit('Error: missing ' + args.bam + '.pbi')
    return args


def clean_up():
    files_to_delete = []
    for f in ['pbalign_alignments.bam', 'pbalign_alignments.bam.pbi', 'pbalign_alignments.bam.bai']:
        if os.path.isfile(f):
            files_to_delete.append(f)
    for f in [f for f in os.listdir('.') if os.path.isfile(f) and f.endswith('.fai')]:
        files_to_delete.append(f)
    for f in [f for f in os.listdir('.') if os.path.isfile(f) and f.endswith('.bt2')]:
        files_to_delete.append(f)
    for f in [f for f in os.listdir('.') if os.path.isfile(f) and f.startswith('illumina_align')]:
        files_to_delete.append(f)
    if files_to_delete:
        print()
        print_command(['rm'] + files_to_delete)
        for f in files_to_delete:
            os.remove(f)


def get_tool_paths(args):
    if not args.samtools:
        args.samtools = shutil.which('samtools')
    if not args.samtools:
        sys.exit('Error: could not find samtools')

    if args.bowtie2:
        args.bowtie2_build = args.bowtie2 + '-build'
    else:
        args.bowtie2 = shutil.which('bowtie2')
        args.bowtie2_build = shutil.which('bowtie2-build')
    if not args.bowtie2:
        sys.exit('Error: could not find bowtie2')
    if not args.bowtie2_build:
        sys.exit('Error: could not find bowtie2-build')

    if args.pitchfork:
        def add_to_env_var(key, new_val):
            if key in os.environ:
                os.environ[key] = new_val + ':' + os.environ[key]
            else:
                os.environ[key] = new_val
        add_to_env_var('LD_LIBRARY_PATH', os.path.join(args.pitchfork, 'lib'))
        add_to_env_var('PATH', os.path.join(args.pitchfork, 'bin'))
        add_to_env_var('PYTHONPATH', os.path.join(args.pitchfork, 'lib', 'python2.7',
                                                  'dist-packages'))

    if not args.pbalign:
        args.pbalign = shutil.which('pbalign')
    if not args.pbalign:
        sys.exit('Error: could not find pbalign')

    if not args.arrow:
        args.arrow = shutil.which('arrow')
    if not args.arrow:
        sys.exit('Error: could not find arrow')

    if not args.freebayes:
        args.freebayes = shutil.which('freebayes')
    if not args.freebayes:
        sys.exit('Error: could not find freebayes')


def make_reads_bam(args):
    command = [args.bax2bam] + args.bax
    run_command_no_output(command)
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


def polish_assembly(fasta, round_num, args):
    align_illumina_reads(fasta, args)
    insert_mean, _ = analyse_all_insert_sizes(args)
    align_pacbio_reads(fasta, args)
    raw_variants_gff = '%03d' % round_num + '_raw_variants.gff'
    genomic_consensus(fasta, args, raw_variants_gff)
    raw_variants = load_variants(raw_variants_gff, fasta, args.large)

    for variant in raw_variants:
        variant.assess_against_illumina_alignments(fasta, args, insert_mean)
    filtered_variants = filter_variants(raw_variants, raw_variants_gff, args)
    polished_fasta = '%03d' % round_num + '_polish.fasta'
    apply_variants(fasta, filtered_variants, polished_fasta)
    print_result(raw_variants, filtered_variants, polished_fasta)
    done = len(filtered_variants) == 0
    return done, polished_fasta


def align_pacbio_reads(fasta, args):
    command = [args.pbalign,
               '--nproc', str(args.threads),
               '--minLength', str(args.min_align_length),
               '--algorithmOptions="--minRawSubreadScore 800 --bestn 1"',
               args.bam,
               fasta,
               'pbalign_alignments.bam']
    run_command_no_output(command)
    files = [f for f in os.listdir('.') if os.path.isfile(f)]
    if 'pbalign_alignments.bam' not in files:
        sys.exit('Error: pbalign failed to make pbalign_alignments.bam')
    if 'pbalign_alignments.bam.pbi' not in files:
        sys.exit('Error: pbalign failed to make pbalign_alignments.bam.pbi')
    if 'pbalign_alignments.bam.bai' not in files:
        sys.exit('Error: pbalign failed to make pbalign_alignments.bam.bai')


def align_illumina_reads(fasta, args):
    index = 'bowtie_index'
    bam = 'illumina_alignments.bam'

    run_command_no_output([args.bowtie2_build, fasta, index])

    bowtie2_command = [args.bowtie2,
                       '--end-to-end',
                       '--very-sensitive',
                       '--threads', str(args.threads),
                       '--no-unal',
                       '-I', '0', '-X', '2000',
                       '-x', index,
                       '-1', args.short1, '-2', args.short2]
    samtools_view_command = [args.samtools, 'view', '-hu', '-']
    samtools_sort_command = [args.samtools, 'sort',
                             '-@', str(args.threads),
                             '-o', 'illumina_alignments.bam',
                             '-']
    print_command(bowtie2_command + ['|'] + samtools_view_command + ['|'] + samtools_sort_command)

    bowtie2 = subprocess.Popen(bowtie2_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    samtools_view = subprocess.Popen(samtools_view_command, stdin=bowtie2.stdout,
                                     stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    bowtie2.stdout.close()
    samtools_sort = subprocess.Popen(samtools_sort_command, stdin=samtools_view.stdout,
                                     stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    samtools_view.stdout.close()
    samtools_sort.communicate()

    run_command_no_output([args.samtools, 'index', bam])


def genomic_consensus(fasta, args, raw_variants_filename):
    subprocess.call([args.samtools, 'faidx', fasta])
    command = [args.arrow,
               'pbalign_alignments.bam',
               '-j', str(args.threads),
               '--noEvidenceConsensusCall', 'reference',
               '-r', fasta,
               '-o', raw_variants_filename]
    run_command_no_output(command)
    files = [f for f in os.listdir('.') if os.path.isfile(f)]
    if raw_variants_filename not in files:
        sys.exit('Error: arrow failed to make ' + raw_variants_filename)


def filter_variants(raw_variants, raw_variants_gff, args):
    """
    This function produces a new GFF file without variants in regions where the Illumina mapping
    is quite solid.
    """
    print()
    print(get_out_header())
    filtered_variants = []
    for variant in raw_variants:
        variant_output = variant.get_out_line() + '\t'
        if variant.large:
            passed = variant.illumina_insert_support
        else:  # substitutions and small indels
            passed = (variant.illumina_alt_percent and
                      variant.illumina_alt_percent >= args.illumina_alt and
                      variant.homo_size_before < args.homopolymer)
        if passed:
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
    in_seqs = collections.OrderedDict(load_fasta(in_fasta))
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
    timestamp = '{:%Y-%m-%d %H:%M:%S}'.format(datetime.datetime.now())
    print('\033[1m' + timestamp + '\033[0m' + '  ' + ' '.join(command), flush=True)


def print_round_header(text):
    print()
    print('\033[1m' + '\033[93m' + '\033[4m' + text + '\033[0m', flush=True)


def print_result(raw_variants, filtered_variants, latest_assembly):
    print_result_line('Unfiltered variants:        ' + str(len(raw_variants)))
    print_result_line('Filtered variants:          ' + str(len(filtered_variants)))
    print_result_line('Polished FASTA:             ' + latest_assembly)


def print_result_line(text):
    print('\033[1m' + '\033[32m' + text + '\033[0m', flush=True)


def print_finished(latest_assembly):
    print()
    result = 'No changes in ' + latest_assembly + '. All done!'
    print('\033[1m' + '\033[32m' + result + '\033[0m', flush=True)


def print_warning(warning):
    print('\033[31m' + warning + '\033[0m', flush=True)


# def run_command_print_output(command):
#     print_command(command)
#     try:
#         process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
#         for line in iter(process.stdout.readline, b''):
#             print(line.rstrip().decode(), flush=True)
#         process.wait()
#     except subprocess.CalledProcessError as e:
#         sys.exit(e.output.decode())


def run_command_no_output(command, shell=False):
    print_command(command)
    try:
        if shell:
            os.system(' '.join(command))
        else:
            subprocess.check_output(command, stderr=subprocess.STDOUT, shell=False)
    except subprocess.CalledProcessError as e:
        sys.exit(e.output.decode())


def load_variants(gff_file, fasta, large_var_size):
    reference = dict(load_fasta(fasta))
    variants = []
    with open(gff_file, 'rt') as gff:
        for line in gff:
            line = line.strip()
            if line and not line.startswith('##'):
                variants.append(PacBioVariant(line, reference, large_var_size))
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


class PacBioVariant(object):
    def __init__(self, gff_line, reference, large_var_size):
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

        full_ref_sequence = reference[self.ref_name]

        # Figure out the change if homopolymer length (if applicable) using these rules:
        # Only indels change homopolymer length
        if self.type != 'insertion' and self.type != 'deletion':
            self.homo_size_before = 0
            self.homo_size_after = 0
        # Any indel which contains multiple different bases doesn't change the homopolymer length
        elif has_multiple_bases(self.ref_seq) or has_multiple_bases(self.variant_seq):
            self.homo_size_before = 0
            self.homo_size_after = 0
        # Insertions only change the homopolymer length if they are inserting the same base
        elif self.type == 'insertion' and self.variant_seq[0] != full_ref_sequence[self.start_pos]:
            self.homo_size_before = 0
            self.homo_size_after = 0
        else:
            self.homo_size_before = homopolymer_size(full_ref_sequence, self.start_pos)
            if self.type == 'insertion':
                self.homo_size_after = self.homo_size_before + len(self.variant_seq)
            else:  # deletion
                self.homo_size_after = max(self.homo_size_before - len(self.ref_seq), 0)

        # Categorise indel variants as small or large. These two categories are assessed
        # differently regarding whether or not to apply them.
        if self.type == 'insertion':
            self.large = len(self.variant_seq) >= large_var_size
        elif self.type == 'deletion':
            self.large = len(self.ref_seq) >= large_var_size
        else:
            self.large = False

        self.ro = 0
        self.ao = 0
        self.illumina_alt_percent = None
        self.actual_illumina_insert = None
        self.expected_illumina_insert = None
        self.illumina_insert_support = None

    def assess_against_illumina_alignments(self, reference_fasta, args, insert_mean):

        # The first way to assess a variant against the illumina alignments is to use freebayes to
        # see how many alternate bases are present at the variant location and what fraction of
        # the bases are alternates. A high alternate fraction indicates that something's a bit
        # screwy and we should probably apply the PacBio variant.
        ref_location = self.ref_name + ':' + str(self.start_pos - 5) + '-' + str(self.end_pos + 5)
        freebayes_command = [args.freebayes,
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
                    alt_occurrences = sum(int(x) for x in
                                          line.split(';AO=')[1].split(';')[0].split(','))
                else:
                    alt_occurrences = 0
                total_occurrences = ref_occurrences + alt_occurrences
                if total_occurrences:
                    alt_fraction = alt_occurrences / total_occurrences
                else:
                    alt_fraction = 0.0
                if self.illumina_alt_percent is None or alt_fraction >= self.illumina_alt_percent:
                    self.ao = alt_occurrences
                    self.ro = ref_occurrences
                    self.illumina_alt_percent = alt_fraction * 100.0

        # The second way to assess a variant against the illumina alignments is to see what the
        # insert size looks like around this variant. We compare the mean insert size to two
        # different values: the overall mean insert (what we'd expect if there is no variant) and
        # an adjusted insert (what we'd expect if there is a variant).
        if self.type == 'insertion':
            self.expected_illumina_insert = insert_mean - len(self.variant_seq)
        elif self.type == 'deletion':
            self.expected_illumina_insert = insert_mean + len(self.ref_seq)
        else:
            self.expected_illumina_insert = insert_mean

        read_pairs = get_read_pairs_around_variant(self, args)
        try:
            self.actual_illumina_insert = statistics.mean(x.first_read.insert_size
                                                          for x in read_pairs)
            if self.large and (self.type == 'insertion' or self.type == 'deletion'):
                diff_if_true = abs(self.actual_illumina_insert - self.expected_illumina_insert)
                diff_if_false = abs(self.actual_illumina_insert - insert_mean)
                self.illumina_insert_support = diff_if_true < diff_if_false
        except statistics.StatisticsError:
            pass

    def get_out_line(self):
        if self.homo_size_before > 1:
            variant_type = 'homo ' + str(self.homo_size_before) + ' -> ' + str(
                self.homo_size_after)
        else:
            variant_type = self.type
            if self.large:
                variant_type = 'large ' + variant_type

        ref_seq = self.ref_seq if self.ref_seq else '.'
        variant_seq = self.variant_seq if self.variant_seq else '.'

        alt_percent = 'n/a' if self.illumina_alt_percent is None \
            else '%.1f' % self.illumina_alt_percent
        actual_insert = 'n/a' if self.actual_illumina_insert is None \
            else '%.1f' % self.actual_illumina_insert
        expected_insert = 'n/a' if self.expected_illumina_insert is None \
            else '%.1f' % self.expected_illumina_insert
        if self.illumina_insert_support is None:
            insert_support = 'n/a'
        else:
            insert_support = 'yes' if self.illumina_insert_support else 'no'

        return '\t'.join([self.ref_name,
                          str(self.start_pos + 1),
                          ref_seq,
                          variant_seq,
                          variant_type,
                          str(self.ao),
                          str(self.ro),
                          alt_percent,
                          actual_insert,
                          expected_insert,
                          insert_support])


def get_out_header():
    return '\033[1m' + '\t'.join(['CONTIG',
                                  'POS',
                                  'REF',
                                  'ALT',
                                  'VARIANT TYPE',
                                  'AO',
                                  'RO',
                                  'AO %',
                                  'ACTUAL INSERT',
                                  'EXPECTED INSERT',
                                  'INSERT SUPPORT',
                                  'RESULT']) + '\033[0m'


class AlignedRead(object):
    def __init__(self, sam_line):
        sam_parts = sam_line.split('\t')
        self.name = sam_parts[0]
        self.sam_flags = int(sam_parts[1])
        self.ref_name = sam_parts[2]
        self.cigar_parts = re.findall(r'\d+\w', sam_parts[5])
        self.insert_size = abs(float(sam_parts[8]))
        self.read_seq = sam_parts[9]

        if self.sam_flags & 64:
            self.strand = '-'
        else:
            self.strand = '+'

        self.start_pos = int(sam_parts[3]) - 1
        self.end_pos = self.start_pos + 1
        for cigar_part in self.cigar_parts:
            self.end_pos += get_ref_shift_from_cigar_part(cigar_part)
        self.mid_point = (self.start_pos + self.end_pos) / 2


class AlignedReadPair(object):
    def __init__(self, name, sam_line):
        self.name = name
        self.first_read = None
        self.second_read = None
        self.add_sam_line(sam_line)

    def add_sam_line(self, sam_line):
        sam_line_parts = sam_line.split('\t')
        sam_flags = int(sam_line_parts[1])
        if sam_flags & 64:
            self.first_read = AlignedRead(sam_line)
        elif sam_flags & 128:
            self.second_read = AlignedRead(sam_line)

    def has_both_reads(self):
        return self.first_read is not None and self.second_read is not None

    def one_read_per_strand(self):
        """
        A normally mapped read pair has one read on each strand. This function returns true if
        that's the case for this pair.
        """
        if not self.has_both_reads():
            return False
        if self.first_read.strand == '+' and self.second_read.strand == '-':
            return True
        if self.first_read.strand == '-' and self.second_read.strand == '+':
            return True
        return False

    def order_first_second_by_ref_pos(self):
        """
        This function makes first_read and second_read in order of reference position, instead of
        being based on sam flags.
        """
        if self.has_both_reads() and self.first_read.start_pos > self.second_read.start_pos:
            self.first_read, self.second_read = self.second_read, self.first_read


def get_read_pairs_around_variant(variant, args):
    read_pairs = {}
    margin = 10000  # Just a big distance to either side to ensure we get all relevant read pairs
    pos_range = variant.ref_name + ':' + str(variant.start_pos - margin) + '-' + \
        str(variant.end_pos + margin)
    command = [args.samtools, 'view', 'illumina_alignments.bam', pos_range]

    try:
        process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        for sam_line in iter(process.stdout.readline, b''):
            name = sam_line.decode().split('\t')[0]
            if not name:
                continue
            if name not in read_pairs:
                read_pairs[name] = AlignedReadPair(name, sam_line)
            else:
                read_pairs[name].add_sam_line(sam_line)
        process.wait()
    except subprocess.CalledProcessError as e:
        sys.exit(e.output.decode())

    read_pairs = [x for x in read_pairs.values() if x.one_read_per_strand()]
    for pair in read_pairs:
        pair.order_first_second_by_ref_pos()

    # Find all read pairs where the variant resides in the gap between the reads.
    filtered_read_pairs = [x for x in read_pairs
                           if x.first_read.end_pos <= variant.start_pos and
                           x.second_read.start_pos >= variant.end_pos]

    # If that was too stringent (e.g. the read pairs overlap so we don't have any cases where the
    # variant is in a gap between them), try again allowing the variant to be in the reads.
    if not filtered_read_pairs:
        filtered_read_pairs = [x for x in read_pairs
                               if x.first_read.mid_point <= variant.start_pos and
                               x.second_read.mid_point >= variant.end_pos]
    return filtered_read_pairs


def analyse_all_insert_sizes(args):
    insert_sizes = []
    command = [args.samtools, 'view', 'illumina_alignments.bam']
    try:
        process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        for sam_line in iter(process.stdout.readline, b''):
            try:
                sam_parts = sam_line.decode().split('\t')
                sam_flags = int(sam_parts[1])
                if sam_flags & 2:
                    insert_size = float(sam_parts[8])
                    if 0.0 < insert_size < 10000.0:  # Just a sanity check...
                        insert_sizes.append(insert_size)
            except (ValueError, IndexError):
                pass
        process.wait()
    except subprocess.CalledProcessError as e:
        sys.exit(e.output.decode())
    try:
        return statistics.mean(insert_sizes), statistics.stdev(insert_sizes)
    except statistics.StatisticsError:
        return 0.0, 0.0
