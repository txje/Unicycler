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
import math
from .misc import add_line_breaks_to_sequence, load_fasta, MyHelpFormatter, print_table, \
    get_percentile_sorted, get_pilon_jar_path, colour


def main():
    args = get_arguments()
    get_tool_paths(args)
    clean_up()

    if args.bax:
        print_round_header('Converting bax.h5 reads to BAM format')
        make_reads_bam(args)
        args.bam = 'subreads.bam'

    current = args.fasta
    round_num = 0
    files = [f for f in os.listdir('.') if os.path.isfile(f)]
    while any(f.startswith('%03d' % round_num) for f in files):
        round_num += 1

    # Do a first round of Illumina alignments just for finding the insert size mean and stdev.
    print_round_header('Determining insert size')
    align_illumina_reads(current, args, 0, 2000, make_index=False)
    min_insert, mean_insert, max_insert = analyse_insert_sizes(args)
    if min_insert == 0 or max_insert == 0:
        sys.exit('Error: could not determine Illumina reads insert size')
    clean_up()
    print_insert_sizes(min_insert, mean_insert, max_insert)

    # The main polishing loop:
    #   1) Run Arrow to generate both small and large variants
    #   2) Apply any small variants which are supported by the Illumina alignments
    #   3) If small variants were applied, go back to step 1 and repeat.
    #   4) Run Pilon to look for larger variants.
    #   5) Assess each large variant with ALE and apply the best one (if any).
    #   6) If a large variant was applied, go back to step 1 and repeat.
    shutil.copy(current, '%03d' % round_num + '_starting_sequence.fasta')
    while True:
        large_variants = []
        while True:
            round_num += 1
            current, change_count, large_variants = polish_small_changes(current, round_num, args,
                                                                         min_insert, max_insert)
            if not change_count:
                break
        round_num += 1
        current, change_count = polish_large_changes(current, round_num, args, large_variants,
                                                     min_insert, max_insert)
        if not change_count:
            break

    shutil.copy(current, 'final_polish.fasta')
    print_finished('final_polish.fasta')


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
    parser.add_argument('--samtools', type=str, default='samtools',
                        help='path to samtools executable')
    parser.add_argument('--bowtie2', type=str, default='bowtie2',
                        help='path to bowtie2 executable')
    parser.add_argument('--bax2bam', type=str, default='bax2bam',
                        help='path to bax2bam executable')
    parser.add_argument('--pbalign', type=str, default='pbalign',
                        help='path to pbalign executable')
    parser.add_argument('--arrow', type=str, default='arrow',
                        help='path to arrow executable')
    parser.add_argument('--freebayes', type=str, default='freebayes',
                        help='path to freebayes executable')
    parser.add_argument('--pilon', type=str, default='pilon*.jar',
                        help='path to pilon jar file')
    parser.add_argument('--java', type=str, default='java',
                        help='path to java executable')
    parser.add_argument('--ale', type=str, default='ALE',
                        help='path to ale executable')

    args = parser.parse_args()
    if not args.bax and not args.bam:
        parser.error('either --bax or --bam is required')
    if args.bax and args.bam:
        parser.error('only one of the following can be used: --bax or --bam')
    if args.bam and not os.path.isfile(args.bam + '.pbi'):
        sys.exit('Error: missing ' + args.bam + '.pbi')
    return args


def clean_up(pbalign_alignments=True, illumina_alignments=True, indices=True, large_variants=True,
             ale_scores=True):
    all_files = [f for f in os.listdir('.') if os.path.isfile(f)]
    files_to_delete = []
    if pbalign_alignments:
        files_to_delete += [f for f in all_files if f.startswith('pbalign_align')]
    if illumina_alignments:
        files_to_delete += [f for f in all_files if f.startswith('illumina_align')]
    if indices:
        files_to_delete += [f for f in all_files if f.endswith('.bt2') or f.endswith('.fai')]
    if large_variants:
        files_to_delete += [f for f in all_files if f.startswith('large_variant_')]
    if ale_scores:
        files_to_delete += [f for f in all_files if f.startswith('ale.out')]
    if files_to_delete:
        print_command(['rm'] + files_to_delete)
        for f in files_to_delete:
            os.remove(f)


def get_tool_paths(args):
    args.samtools = shutil.which(args.samtools)
    if not args.samtools:
        sys.exit('Error: could not find samtools')

    args.bowtie2 = shutil.which(args.bowtie2)
    if not args.bowtie2:
        sys.exit('Error: could not find bowtie2')
    args.bowtie2_build = shutil.which(args.bowtie2 + '-build')
    if not args.bowtie2_build:
        sys.exit('Error: could not find bowtie2-build (it should be in the same place as bowtie2)')

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

    args.pbalign = shutil.which(args.pbalign)
    if not args.pbalign:
        sys.exit('Error: could not find pbalign')

    args.arrow = shutil.which(args.arrow)
    if not args.arrow:
        sys.exit('Error: could not find arrow')

    args.freebayes = shutil.which(args.freebayes)
    if not args.freebayes:
        sys.exit('Error: could not find freebayes')

    if args.pilon == 'pilon*.jar':
        args.pilon = get_pilon_jar_path(None)
    else:
        args.pilon = get_pilon_jar_path(args.pilon)
    if not args.pilon:
        sys.exit('Error: could not find pilon jar file')

    args.java = shutil.which(args.java)
    if not args.java:
        sys.exit('Error: could not find java')

    args.ale = shutil.which(args.ale)
    if not args.ale:
        sys.exit('Error: could not find ale')


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


def polish_small_changes(fasta, round_num, args, min_insert, max_insert):
    print_round_header('Round ' + str(round_num) + ': small variants')

    raw_variants_gff = '%03d' % round_num + '_1_raw_variants.gff'
    filtered_variants_gff = '%03d' % round_num + '_2_filtered_variants.gff'
    polished_fasta = '%03d' % round_num + '_3_polish.fasta'

    align_pacbio_reads(fasta, args)
    genomic_consensus(fasta, args, raw_variants_gff)
    raw_variants = load_variants_from_gff(raw_variants_gff, fasta, args.large)

    # For assessing small variants, we want end-to-end alignment (not local).
    align_illumina_reads(fasta, args, min_insert, max_insert, local=False)

    large_variants = [x for x in raw_variants if x.large]
    small_variants = [x for x in raw_variants if not x.large]
    for variant in small_variants:
        variant.assess_against_illumina_alignments(fasta, args)
    clean_up()

    filtered_variants = filter_small_variants(small_variants, raw_variants_gff,
                                              filtered_variants_gff, args)
    apply_variants(fasta, filtered_variants, polished_fasta)

    # if large_variants:
    #     # If small variants were applied, then large variants are discarded because another small
    #     # change round is coming up. If no small variants were applied, then we'll use the large
    #     # variants in the next round.
    #     print()
    #     print('Unused large variants:' if filtered_variants
    #           else 'Large variants for the next round:')
    #     print_simple_large_variant_table(large_variants)

    print_result(filtered_variants, polished_fasta)
    return polished_fasta, len(filtered_variants), large_variants


def polish_large_changes(fasta, round_num, args, large_variants, min_insert, max_insert):
    print_round_header('Round ' + str(round_num) + ': larger variants and Pilon')

    large_variants += get_pilon_variants(fasta, args, min_insert, max_insert, round_num)
    if not large_variants:
        print('\nNo variants found')
        return fasta, 0

    raw_variants_gff = '%03d' % round_num + '_1_raw_variants.gff'
    filtered_variants_gff = '%03d' % round_num + '_2_filtered_variants.gff'
    ale_outputs = '%03d' % round_num + '_3_ALE_output'
    polished_fasta = '%03d' % round_num + '_4_polish.fasta'

    save_large_variants(large_variants, raw_variants_gff)
    open(filtered_variants_gff, 'a').close()
    open(ale_outputs, 'a').close()

    initial_ale_score = run_ale(fasta, args, min_insert, max_insert, ale_outputs)
    best_ale_score = initial_ale_score

    best_modification = None
    applied_variant = []
    for i, variant in enumerate(large_variants):
        modified_assembly = 'large_variant_' + str(i+1) + '.fasta'
        apply_variants(fasta, [variant], modified_assembly)
        variant.ale_score = run_ale(modified_assembly, args, min_insert, max_insert, ale_outputs)
        if variant.ale_score > best_ale_score:
            best_ale_score = variant.ale_score
            best_modification = modified_assembly
            applied_variant = [variant]
            save_large_variants(applied_variant, filtered_variants_gff)

    if best_modification:
        os.rename(best_modification, polished_fasta)
    else:
        shutil.copyfile(fasta, polished_fasta)
    clean_up()

    print_large_variant_table(large_variants, best_ale_score, initial_ale_score)

    print_result(applied_variant, polished_fasta)
    return polished_fasta, len(applied_variant)


def run_ale(fasta, args, min_insert, max_insert, all_ale_outputs):
    """
    ALE is run in --metagenome mode because this polishing script is presumed to be used on
    completed bacterial genomes, where each contig is different replicon (chromosome or plasmid)
    with potentially different depth.
    """
    ale_output = 'ale.out'
    ale_score = float('-inf')
    previous_output_exists = os.path.getsize(all_ale_outputs) > 0

    align_illumina_reads(fasta, args, min_insert, max_insert, local=False, keep_unaligned=True)

    run_command_no_output([args.ale,
                           '--nout',
                           '--metagenome',
                           'illumina_alignments.bam', fasta, ale_output])
    if not os.path.isfile(ale_output):
        sys.exit('Error: ALE did not generate ' + ale_output)

    with open(ale_output, 'rt') as ale_output_file:
        with open(all_ale_outputs, 'at') as all_ale_outputs_file:
            if previous_output_exists:
                all_ale_outputs_file.write('\n\n\n\n\n')
            for line in ale_output_file:
                all_ale_outputs_file.write(line)
                if 'ALE_score:' in line and ale_score == float('-inf'):
                    ale_score = float(line.split('ALE_score:')[1].strip().split()[0])

    clean_up(large_variants=False)
    return ale_score


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


def align_illumina_reads(fasta, args, min_insert, max_insert, make_index=True, local=False,
                         keep_unaligned=False):
    index = 'bowtie_index'
    bam = 'illumina_alignments.bam'

    run_command_no_output([args.bowtie2_build, fasta, index])

    if local:
        bowtie2_command = [args.bowtie2, '--local', '--very-sensitive-local']
    else:
        bowtie2_command = [args.bowtie2, '--end-to-end', '--very-sensitive']
    if not keep_unaligned:
        bowtie2_command += ['--no-unal']
    bowtie2_command += ['--threads', str(args.threads),
                        '-I', str(min_insert), '-X', str(max_insert),
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

    if make_index:
        run_command_no_output([args.samtools, 'index', bam])


def pilon(fasta, args, raw_pilon_changes_filename):
    pilon_command = [args.java, '-jar', args.pilon,
                     '--genome', fasta,
                     '--frags', 'illumina_alignments.bam',
                     '--fix', 'all',
                     '--changes',
                     '--outdir', 'temp_pilon']
    run_command_no_output(pilon_command)

    pilon_changes = os.path.join('temp_pilon', 'pilon.changes')
    if not os.path.isfile(pilon_changes):
        sys.exit('Pilon did not produce pilon.changes')
    shutil.copy(pilon_changes, raw_pilon_changes_filename)
    shutil.rmtree('temp_pilon')


def get_pilon_variants(fasta, args, min_insert, max_insert, round_num):
    # Pilon needs local alignment to help spot misassembly regions and unaligned reads to use
    # when reassembling.
    align_illumina_reads(fasta, args, min_insert, max_insert, local=True, keep_unaligned=True)
    raw_pilon_changes = '%03d' % round_num + '_pilon.changes'
    pilon(fasta, args, raw_pilon_changes)
    clean_up()
    variants = load_variants_from_pilon_changes(raw_pilon_changes, fasta, args.large)
    if not variants:
        os.remove(raw_pilon_changes)
    return variants


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


def filter_small_variants(raw_variants, raw_variants_gff, filtered_variants_gff, args):
    """
    This function produces a new GFF file without variants in regions where the Illumina mapping
    is quite solid.
    """
    filtered_variants = []
    variant_rows = []
    for variant in raw_variants:
        variant_row = variant.get_output_row()
        passed = (variant.illumina_alt_percent and                       # must have an alt%
                  variant.illumina_alt_percent >= args.illumina_alt and  # must have a big alt%
                  variant.homo_size_before < args.homopolymer)           # must not be homopolymer
        if passed:
            filtered_variants.append(variant)
            variant_row.append('PASS')
        else:
            variant_row.append('FAIL')
        variant_rows.append(variant_row)
    print_small_variant_table(variant_rows)

    with open(filtered_variants_gff, 'wt') as new_gff:
        with open(raw_variants_gff, 'rt') as old_gff:
            for line in old_gff:
                if line.startswith('##'):
                    new_gff.write(line)
        for variant in filtered_variants:
            new_gff.write(variant.original_gff_line + '\n')
    return filtered_variants


def save_large_variants(variants, filename):
    with open(filename, 'wt') as output_file:
        for variant in variants:
            output_file.write(variant.get_original_line() + '\n')


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
    command = [os.path.basename(command[0])] + command[1:]  # Remove path from program name
    print('\033[1m' + timestamp + '\033[0m' + '  ' + ' '.join(command), flush=True)


def print_round_header(text):
    print('\n')
    print('\033[1m' + '\033[93m' + '\033[4m' + text + '\033[0m', flush=True)


def print_result(variants, fasta):
    print()
    print(str(len(variants)) + (' variant' if len(variants) == 1 else ' variants') + ' applied')
    print('Result:', fasta, flush=True)


def print_finished(fasta):
    print()
    result = 'All done! Final assembly: ' + fasta
    print('\033[1m' + '\033[32m' + result + '\033[0m', flush=True)


def print_warning(warning):
    print('\033[31m' + warning + '\033[0m', flush=True)


def run_command_no_output(command, shell=False):
    print_command(command)
    try:
        if shell:
            os.system(' '.join(command))
        else:
            subprocess.check_output(command, stderr=subprocess.STDOUT, shell=False)
    except subprocess.CalledProcessError as e:
        sys.exit(e.output.decode())


def load_variants_from_gff(gff_file, fasta, large_var_size):
    reference = dict(load_fasta(fasta))
    variants = []
    with open(gff_file, 'rt') as gff:
        for line in gff:
            line = line.strip()
            if line and not line.startswith('##'):
                variants.append(Variant(reference, large_var_size, gff_line=line))
    return variants


def load_variants_from_pilon_changes(pilon_changes_file, fasta, large_var_size):
    reference = dict(load_fasta(fasta))
    variants = []
    with open(pilon_changes_file, 'rt') as changes:
        for line in changes:
            line = line.strip()
            if line:
                variants.append(Variant(reference, large_var_size, changes_line=line))
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


class Variant(object):
    def __init__(self, reference, large_var_size, gff_line=None, changes_line=None):
        self.original_gff_line = gff_line
        self.original_changes_line = changes_line

        if gff_line:
            # https://github.com/PacificBiosciences/GenomicConsensus/blob/master/doc/VariantsGffSpecification.rst
            line_parts = gff_line.split('\t')
            attributes = {x.split('=')[0]: x.split('=')[1] for x in line_parts[8].split(';')}
            self.source = 'Arrow'
            self.ref_name = line_parts[0]
            self.type = line_parts[2]
            self.ref_seq = attributes['reference'].replace('.', '')
            self.variant_seq = attributes['variantSeq'].replace('.', '')
            self.start_pos = int(line_parts[3]) - 1
            self.end_pos = self.start_pos + len(self.ref_seq)

        elif changes_line:
            # https://github.com/broadinstitute/pilon/wiki/Output-File-Descriptions
            line_parts = changes_line.split(' ')
            self.source = 'Pilon'
            self.ref_name = line_parts[0].split(':')[0]
            self.ref_seq = line_parts[2].replace('.', '')
            self.variant_seq = line_parts[3].replace('.', '')
            self.start_pos = int(line_parts[0].split(':')[1].split('-')[0]) - 1
            self.end_pos = self.start_pos + len(self.ref_seq)
            if len(self.ref_seq) > len(self.variant_seq):
                self.type = 'deletion'
            elif len(self.ref_seq) < len(self.variant_seq):
                self.type = 'insertion'
            else:
                self.type = 'substitution'

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
        # differently regarding whether or not to apply them. All Pilon changes are 'large',
        # because it is run only asking for 'local' fixes, not 'bases' fixes.
        if self.source == 'Pilon':
            self.large = True
        if self.type == 'insertion':
            self.large = len(self.variant_seq) >= large_var_size
        elif self.type == 'deletion':
            self.large = len(self.ref_seq) >= large_var_size
        else:
            self.large = False

        self.ro = 0
        self.ao = 0
        self.illumina_alt_percent = None
        self.ale_score = float('-inf')

    def assess_against_illumina_alignments(self, reference_fasta, args):
        """
        To assess a variant against the illumina alignments, we use freebayes to see how many
        alternate bases are present at the variant location and what fraction of the bases are
        alternates. A high alternate fraction indicates that something's a bit screwy and we
        should probably apply the PacBio polishing suggestion.
        """
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

    def get_output_row(self):
        if self.homo_size_before > 1:
            variant_type = 'homo ' + str(self.homo_size_before) + ' -> ' + str(
                self.homo_size_after)
        else:
            variant_type = self.type
            if self.large:
                variant_type = 'large ' + variant_type

        ref_seq = self.ref_seq if self.ref_seq else '.'
        variant_seq = self.variant_seq if self.variant_seq else '.'
        # noinspection PyStringFormat
        alt_percent = 'n/a' if self.illumina_alt_percent is None \
            else '%.1f' % self.illumina_alt_percent

        return [self.ref_name, str(self.start_pos + 1), ref_seq, variant_seq, variant_type,
                str(self.ao), str(self.ro), alt_percent]

    def get_original_line(self):
        if self.original_gff_line:
            return self.original_gff_line
        else:
            return self.original_changes_line


def print_small_variant_table(rows):
    print()
    header = ['Contig', 'Pos', 'Ref', 'Alt', 'Type', 'AO', 'RO', 'AO%', 'Result']
    print_table([header] + rows, alignments='LRLLLRRRR')


def print_simple_large_variant_table(variants):
    table = [['Contig', 'Pos', 'Ref', 'Alt']]
    for v in variants:
        ref_seq = v.ref_seq if v.ref_seq else '.'
        variant_seq = v.variant_seq if v.variant_seq else '.'
        table.append([v.ref_name, str(v.start_pos), ref_seq, variant_seq])
    print_table(table, alignments='LRLL')


def print_large_variant_table(variants, best_ale_score, initial_ale_score):
    print()
    table = [['Source', 'Contig', 'Pos', 'Ref', 'Alt', 'ALE score']]
    text_colour = 'green' if initial_ale_score == best_ale_score else 'red'
    table.append(['No variant', '', '', '', '', colour('%.6f' % initial_ale_score, text_colour)])
    for v in variants:
        ref_seq = v.ref_seq if v.ref_seq else '.'
        variant_seq = v.variant_seq if v.variant_seq else '.'
        text_colour = 'green' if v.ale_score == best_ale_score else 'red'
        ale_score_str = colour('%.6f' % v.ale_score, text_colour)
        table.append([v.source, v.ref_name, str(v.start_pos), ref_seq, variant_seq, ale_score_str])
    print_table(table, alignments='LLRLLR')


def analyse_insert_sizes(args):
    insert_sizes = []
    command = [args.samtools, 'view', 'illumina_alignments.bam']
    try:
        process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        for sam_line in iter(process.stdout.readline, b''):
            sam_line = sam_line.decode()
            try:
                sam_parts = sam_line.split('\t')
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
    insert_sizes = sorted(insert_sizes)

    min_insert = math.floor(get_percentile_sorted(insert_sizes, 2.5))
    mean_insert = statistics.mean(insert_sizes)
    max_insert = math.ceil(get_percentile_sorted(insert_sizes, 97.5))

    return min_insert, mean_insert, max_insert


def print_insert_sizes(min_insert, mean_insert, max_insert):
    print()
    print(' 2.5th percentile:', min_insert)
    print(' mean insert size:', '\033[1m' + '%.1f' % mean_insert + '\033[0m')
    print('97.5th percentile:', max_insert, flush=True)
