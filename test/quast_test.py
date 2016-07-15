#!/usr/bin/env python3
"""
Tool for testing QUAST accuracy.
It fragments an input genome into random pieces and then gives those to QUAST to assess. Since it
did not introduce any errors, QUAST should not report any misassemblies or errors.

Usage:
assembler_comparison.py --reference path/to/reference.fasta

Author: Ryan Wick
email: rrwick@gmail.com
"""

import random
import os
import argparse
import subprocess
import shutil


def main():
    args = get_args()
    references = load_fasta(args.reference)
    quast_results = create_quast_results_table()
    ref_length, _ = get_fasta_length_and_seq_count(args.reference)

    while True:
        frag_filename = os.path.abspath('fragments.fasta')
        frag_file = open(frag_filename, 'wt')
        total_fragment_count = 0
        frag_size = random.randint(1000, ref_length // 8)
        print('\nFragmenting to approximately ' + str(frag_size) + ' bp')

        for ref in references:
            ref_seq = ref[1]

            # Randomly rotate the sequence.
            random_start = random.randint(0, len(ref_seq) - 1)
            rotated = ref_seq[random_start:] + ref_seq[:random_start]

            # Randomly fragment the sequence.
            fragment_count = int(round(len(rotated) / frag_size))
            fragment_positions = [0] + [random.randint(1, len(rotated) - 2) for _ in
                                        range(fragment_count)] + [len(rotated)]
            fragment_positions = sorted(list(set(fragment_positions)))
            for i, start_pos in enumerate(fragment_positions[:-1]):
                end_pos = fragment_positions[i + 1]
                frag_file.write('>' + str(i+1) + '\n')
                frag_file.write(rotated[start_pos:end_pos] + '\n')
                total_fragment_count += 1

        frag_file.close()

        run_quast(frag_filename, args, quast_results, total_fragment_count, frag_size)
        for item in os.listdir('.'):
            if item.startswith('fragments.fasta'):
                os.remove(item)


def get_args():
    """
    Specifies the command line arguments required by the script.
    """
    parser = argparse.ArgumentParser(description='QUAST tester')
    parser.add_argument('--reference', type=str, required=True,
                        help='The reference genome to fragment and assess in QUAST')
    parser.add_argument('--threads', type=int, required=False, default=8,
                        help='Number of CPU threads')

    args = parser.parse_args()
    args.reference = os.path.abspath(args.reference)
    return args


def run_quast(assembly, args, all_quast_results, fragment_count, frag_size):
    reference_name = get_reference_name_from_filename(args.reference)

    ref_length, ref_count = get_fasta_length_and_seq_count(args.reference)
    quast_line = [reference_name, str(ref_length), str(ref_count),
                  'fragmented reference', str(fragment_count), '', '', '', '']

    run_name = assembly + ', ' + str(frag_size) + 'bp fragments'
    run_dir_name = run_name.replace(', ', '_').replace(' ', '_').replace('.', '_')

    print('\nRunning QUAST for', reference_name, flush=True)
    quast_dir = os.path.join('quast_results', run_dir_name)
    this_quast_results = os.path.join(quast_dir, 'transposed_report.tsv')

    quast_command = ['quast.py',
                     assembly,
                     '-R', args.reference,
                     '-o', quast_dir,
                     '-l', '"' + run_name.replace(',', '') + '"',
                     '--threads', str(args.threads),
                     '--no-snps']
    print(' '.join(quast_command))
    try:
        subprocess.check_output(quast_command, stderr=subprocess.STDOUT)
    except subprocess.CalledProcessError as e:
        print('QUAST encountered an error:\n' + e.output.decode(), flush=True)
        quast_line += [''] * 46
    else:
        with open(this_quast_results, 'rt') as results:
            results.readline()  # header line
            quast_line += results.readline().split('\t')[1:]

    with open(all_quast_results, 'at') as all_results:
        all_results.write('\t'.join(quast_line))

    shutil.rmtree(quast_dir)


def get_reference_name_from_filename(reference):
    return reference.split('/')[-1].split('.')[0]


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


def create_quast_results_table():
    quast_results_filename = 'quast_results.tsv'
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
                            "Largest alignment\t"
                            "NA50\t"
                            "NGA50\t"
                            "NA75\t"
                            "NGA75\t"
                            "LA50\t"
                            "LGA50\t"
                            "LA75\t"
                            "LGA75\n")
        quast_results.close()
    return os.path.abspath(quast_results_filename)


if __name__ == '__main__':
    main()
