"""
Functions relating to Pilon.

Author: Ryan Wick
email: rrwick@gmail.com
"""

import os
import subprocess
import statistics
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
    raw_sam_filename = os.path.join(polish_dir, 'alignments_raw.sam')
    bowtie2_command = [bowtie2_path, '-x', polish_fasta_filename, '-1', short_1, '-2', short_2,
                       '--end-to-end', '--very-sensitive', '--threads', str(threads),
                       '--no-discordant', '--no-mixed', '--no-unal', '-I', '0', '-X', '1000']
    samtools_view_command = [samtools_path, 'view', '-@', str(threads), '-f', '2', '-h', '-',
                             '-o', raw_sam_filename]
    bowtie2 = subprocess.Popen(bowtie2_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    samtools_view = subprocess.Popen(samtools_view_command, stdin=bowtie2.stdout,
                                     stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    _, err = samtools_view.communicate()

    # # CHECK TO MAKE SURE THAT alignments.sam EXISTS AND LOOKS CORRECT.
    # if err:
    #     if verbosity > 1:
    #         print('\nAn error occurred during alignment:\n' + err.decode())
    #     raise CannotPolish

    # Loop through alignments.sam once to collect the insert sizes.
    insert_sizes = []
    raw_sam = open(raw_sam_filename, 'r')
    for sam_line in raw_sam:
        if sam_line.startswith('@'):
            continue
        sam_parts = sam_line.split('\t')
        if len(sam_parts) > 8:
            insert_size = int(sam_parts[8])
            if insert_size > 0:
                insert_sizes.append(insert_size)
    raw_sam.close()
    insert_mean = statistics.mean(insert_sizes)
    insert_std_dev = statistics.stdev(insert_sizes)
    min_insert = max(insert_mean - (2 * insert_std_dev), 0.0)
    max_insert = insert_mean + (2 * insert_std_dev)

    # Produce a new sam file with only the pairs with an appropriate insert size.
    filtered_sam_filename = os.path.join(polish_dir, 'alignments_filtered.sam')
    filtered_sam = open(filtered_sam_filename, 'w')
    raw_sam = open(raw_sam_filename, 'r')
    for sam_line in raw_sam:
        if sam_line.startswith('@'):
            filtered_sam.write(sam_line)
        sam_parts = sam_line.split('\t')
        if len(sam_parts) > 8:
            insert_size = abs(int(sam_parts[8]))
            if min_insert <= insert_size <= max_insert:
                filtered_sam.write(sam_line)
    raw_sam.close()
    filtered_sam.close()

    # Sort and index the alignments.
    bam = os.path.join(polish_dir, 'alignments.bam')
    samtools_sort_command = [samtools_path, 'sort', '-@', str(threads), filtered_sam_filename,
                             '-o', bam]
    samtools_sort = subprocess.Popen(samtools_sort_command, stdout=subprocess.PIPE,
                                     stderr=subprocess.PIPE)
    _, err = samtools_sort.communicate()
    if err:
        if verbosity > 1:
            print('\nsamtools encountered an error:\n' + err.decode())
        raise CannotPolish
    samtools_index_command = [samtools_path, 'index', bam]
    samtools_index = subprocess.Popen(samtools_index_command, stdout=subprocess.PIPE,
                                      stderr=subprocess.PIPE)
    _, err = samtools_index.communicate()
    if err:
        if verbosity > 1:
            print('\nsamtools encountered an error:\n' + err.decode())
        raise CannotPolish

    # TO DO: RUN PILON!
    # TO DO: RUN PILON!
    # TO DO: RUN PILON!
    # TO DO: RUN PILON!
    # TO DO: RUN PILON!
    # TO DO: RUN PILON!
    # TO DO: RUN PILON!
    # TO DO: RUN PILON!
    # TO DO: RUN PILON!
    # TO DO: RUN PILON!
    # TO DO: RUN PILON!
    # TO DO: RUN PILON!
    # TO DO: RUN PILON!
    # TO DO: RUN PILON!




