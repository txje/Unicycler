"""
Functions relating to Pilon.

Author: Ryan Wick
email: rrwick@gmail.com
"""

import os
import subprocess
from .misc import float_to_str


class CannotPolish(Exception):
    pass


def polish_with_pilon(graph, bowtie2_path, bowtie2_build_path, pilon_path, samtools_path,
                      min_polish_size, polish_dir, verbosity, short_1, short_2, threads):
    """
    Runs Pilon on the graph to hopefully fix up small mistakes.
    """
    segments_to_polish = [x for x in graph.segments.values() if x.get_length() >= min_polish_size]
    if not segments_to_polish:
        raise CannotPolish

    polish_fasta_filename = os.path.join(polish_dir, 'polish.fasta')
    polish_fasta = open(polish_fasta_filename, 'w')
    for segment in segments_to_polish:
        polish_fasta.write('>' + str(segment.number) + '\n')
        polish_fasta.write(segment.forward_sequence)
        polish_fasta.write('\n')
    polish_fasta.close()

    # Prepare the FASTA for Bowtie2 alignment.
    command = [bowtie2_build_path, polish_fasta_filename, polish_fasta_filename]
    process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    _, err = process.communicate()
    index_made = any(x.endswith('.bt2') for x in os.listdir(polish_dir))
    if not index_made:
        if verbosity > 1:
            if err:
                print('\nbowtie2-build encountered an error:\n' + err.decode())
            else:
                print('\nbowtie2-build failed to build an index')
        raise CannotPolish

    # Perform the alignment by piping the Bowtie output into Samtools.
    if verbosity > 1:
        print('Aligning short reads to assembly using Bowtie2...')
    bam = os.path.join(polish_dir, 'alignments')
    bowtie2_command = [bowtie2_path, '-x', polish_fasta_filename, '-1', short_1, '-2', short_2,
                       '--end-to-end', '--very-sensitive', '--threads', str(threads)]
    samtools_view_command = [samtools_path, 'view', '-@', str(threads), '-F', '12', '-b', '-']
    samtools_sort_command = [samtools_path, 'sort', '-@', str(threads), '-', bam]
    bowtie2 = subprocess.Popen(bowtie2_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    samtools_view = subprocess.Popen(samtools_view_command, stdin=bowtie2.stdout,
                                     stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    samtools_sort = subprocess.Popen(samtools_sort_command, stdin=samtools_view.stdout,
                                     stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    _, err = samtools_sort.communicate()
    if err:
        if verbosity > 1:
            print('\nAn error occurred during alignment:\n' + err.decode())
        raise CannotPolish
    samtools_index_command = [samtools_path, 'index', bam + '.bam']
    samtools_index = subprocess.Popen(samtools_index_command, stdin=bowtie2.stdout,
                                      stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    _, err = samtools_index.communicate()
    if err:
        if verbosity > 1:
            print('\nsamtools encountered an error:\n' + err.decode())
        raise CannotPolish
