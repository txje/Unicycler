#!/usr/bin/env python
'''
My totally awesome hybrid assembler.
'''

from __future__ import print_function
from __future__ import division
import argparse
import subprocess
import os
import sys
import shutil
import gzip
import math
from assembly_graph import AssemblyGraph, load_sam_alignments
from semi_global_long_read_aligner import quit_with_error, add_aligning_arguments, \
                                          fix_up_arguments


spades_path = '/Users/Ryan/Applications/SPAdes-3.7.1-Darwin/bin/spades.py'
graphmap_path = '/Users/Ryan/Applications/graphmap/bin/Mac/graphmap'
starting_kmer_fraction = 0.2 # Relative to median read length
max_kmer_fraction = 0.9 # Relative to median read length
kmer_count = 10
read_depth_filter_cutoff = 0.25


def main():
    '''
    Script execution starts here.
    '''
    args = get_arguments()
    make_output_directory(args.out)
    log_filepath = make_log_file(args.out)
    log_file = open(log_filepath, 'w')
    check_file_exists(args.short1)
    check_file_exists(args.short2)
    check_file_exists(args.long)

    # Produce a SPAdes assembly graph with a k-mer that balances contig length and connectivity.
    assembly_graph, kmer = get_best_spades_graph(args.short1, args.short2, args.out, log_file,
                                                 args.read_depth_filter)
    spades_assembly_path = os.path.join(args.out, '01_spades_assembly')
    assembly_graph.save_to_gfa(spades_assembly_path + '.gfa')
    assembly_graph.save_to_fasta(spades_assembly_path + '.fasta')

    # Get a list of single-copy nodes.
    assembly_graph.determine_copy_depth(1000)
    alignment_path = os.path.join(args.out, '02_long_read_alignment')
    assembly_graph.save_to_gfa(os.path.join(alignment_path, 'unbridged_graph.gfa'))
    # TO DO: GET SINGLE COPY NODES

    # Produce an initial set of bridges using SPAdes contig paths.
    # TO DO
    # TO DO
    # TO DO
    # TO DO
    # TO DO

    # If all long reads are available now, then we do the entire process in one pass.
    if args.long:
        # TO DO: ALIGN LONG READS TO GRAPH
        # TO DO: PRODUCE BRIDGES USING LONG READ ALIGNMENTS
        # TO DO: APPLY THE BRIDGES TO THE GRAPH
        # TO DO: SAVE THE RESULTS

    # If we are getting long reads incrementally, then we do the process iteratively.
    elif args.long_dir:
        finished = False
        while not finished:
            # TO DO: WAIT FOR NEW READS TO BECOME AVAILABLE
            # TO DO: ALIGN LONG READS TO GRAPH
            # TO DO: PRODUCE BRIDGES USING LONG READ ALIGNMENTS
            # TO DO: APPLY THE BRIDGES TO THE GRAPH
            # TO DO: SAVE THE RESULTS
            # TO DO: ASSESS WHETHER THE PROCESS IS COMPLETE

def get_arguments():
    '''
    Parse the command line arguments.
    '''
    parser = argparse.ArgumentParser(description='Hybrid Assembler')

    parser.add_argument('-1', required=True, default=argparse.SUPPRESS,
                        help='FASTQ file of short reads (first reads in each pair).')
    parser.add_argument('-2', required=True, default=argparse.SUPPRESS,
                        help='FASTQ file of short reads (second reads in each pair).')
    parser.add_argument('--long', required=False, default=argparse.SUPPRESS,
                        help='FASTQ file of long reads, if all reads are available at start.')
    parser.add_argument('--long_dir', required=False, default=argparse.SUPPRESS,
                        help='Directory where FASTQ files will be deposited.')
    parser.add_argument('--out', required=True, default=argparse.SUPPRESS,
                        help='Output directory')
    parser.add_argument('--read_depth_filter', type=float, required=False, default=0.25,
                        help='Minimum allowed read depth, expressed as a fraction of the median'
                             'read depth. Graph segments with less depth will be removed.')

    # Add the arguments for the aligner, but suppress the help text.
    add_aligning_arguments(parser, True)

    args = parser.parse_args()
    fix_up_arguments(args)

    try:
        args.long
    except AttributeError:
        args.long = None
    try:
        args.long_dir
    except AttributeError:
        args.long_dir = None

    if not args.long and not args.long_dir:
        quit_with_error('either --long or --long_dir is required')
    if args.long and args.long_dir:
        quit_with_error('--long and --long_dir cannot both be used')

    return args


def make_output_directory(outdir):
    '''
    Creates the output directory, if it doesn't already exist.  If it does exist, that's okay as
    long as it's empty, in which case this function does nothing.  If it exists and isn't empty,
    this function quits with an error.
    '''
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    elif os.listdir(outdir):
        print('Error: the output directory (' + outdir + ') is not empty.', file=sys.stderr)
        quit()

def check_file_exists(filename):
    '''
    If the given file exists, this function does nothing.  If it doesn't exist, it displays an
    error and quits.
    '''
    if not os.path.isfile(filename):
        print('Error: could not find ' + filename, file=sys.stderr)
        quit()

def get_best_spades_graph(short1, short2, outdir, log_file, read_depth_filter):
    '''
    This function tries a SPAdes assembly at different kmers and returns the best.
    'The best' is defined as the smallest dead-end count after low-depth filtering.  If multiple
    graphs have the same dead-end count (e.g. zero!) then the highest kmer is used.
    '''
    kmer_range = get_kmer_range(short1, short2, log_file)
    spades_dir = os.path.join(outdir, '01_spades_assembly')
    read_correction_dir = os.path.join(spades_dir, 'read_correction')
    read_files = spades_read_correction(short1, short2, read_correction_dir, log_file)
    shutil.rmtree(read_correction_dir)
    assem_dir = os.path.join(spades_dir, 'assembly')
    print('Conducting SPAdes assemblies...\n')
    print('kmer\tsegments\tdead ends\tconnected_components\tscore')
    best_score = 0.0
    for i, kmer in enumerate(kmer_range):
        graph_file = spades_assembly(read_files, assem_dir, log_file, kmer_range[:i+1])
        assembly_graph = AssemblyGraph(graph_file, kmer)
        clean_graph(assembly_graph, read_depth_filter)
        assembly_graph.save_to_fastg(os.path.join(spades_dir,
                                                  'k' + str(kmer) + '_assembly_graph_cleaned.gfa'))
        dead_ends, connected_components, segment_count = get_graph_info(assembly_graph)
        score = 1.0 / (segment_count * ((dead_ends + 1) ** 2))
        print(str(kmer) + '\t' + str(segment_count) + '\t' + str(dead_ends) +
              '\t' + str(connected_components) + '\t' + str(score))
        if score >= best_score:
            best_kmer = kmer
            best_score = score
            best_assembly_graph = assembly_graph
    shutil.rmtree(assem_dir)
    print('\nBest kmer: ' + str(best_kmer))
    return best_assembly_graph, best_kmer

def clean_graph(assembly_graph, read_depth_filter):
    '''
    For the given graph, this function normalises the read depths, filters based on depth and
    simplifies some graph structures.
    '''
    assembly_graph.repair_four_way_junctions()
    # TO DO: FIX UP SELF REV COMP NODE STRUCTURES
    assembly_graph.filter_by_read_depth(read_depth_filter)
    assembly_graph.filter_homopolymer_loops()
    assembly_graph.normalise_read_depths()

def get_graph_info(assembly_graph):
    '''
    This function returns some graph information that will be useful for determining which kmer is
    best.
    '''
    dead_ends = assembly_graph.total_dead_end_count()
    connected_components = len(assembly_graph.get_connected_components())
    segment_count = len(assembly_graph.segments)
    return dead_ends, connected_components, segment_count

def spades_read_correction(short1, short2, outdir, log_file):
    '''
    This runs SPAdes with the --only-error-correction option.
    '''
    short1_no_extension = strip_read_extensions(short1)
    short2_no_extension = strip_read_extensions(short2)
    print('SPAdes read error correction... ', end='')
    sys.stdout.flush()
    log_file.write('SPAdes read error correction\n')
    log_file.write('----------------------------\n')
    command = [spades_path, '-1', short1, '-2', short2, '--only-error-correction', '-o', outdir]
    process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out, err = process.communicate()
    if len(err) > 0:
        print('\nSPAdes encountered an error:', file=sys.stderr)
        print(err, file=sys.stderr)
        quit()
    log_file.write(out)
    log_file.write('\n\n\n\n\n\n\n\n\n\n\n')
    corrected_dir = os.path.join(outdir, 'corrected')
    parent_dir = os.path.dirname(outdir)
    corrected_1 = os.path.join(parent_dir, 'corrected_reads_1.fastq.gz')
    corrected_2 = os.path.join(parent_dir, 'corrected_reads_2.fastq.gz')
    corrected_u = os.path.join(parent_dir, 'corrected_reads_u.fastq.gz')
    files = os.listdir(corrected_dir)
    for spades_file in files:
        file_path = os.path.join(corrected_dir, spades_file)
        if short1_no_extension in spades_file:
            shutil.move(file_path, corrected_1)
        elif short2_no_extension in spades_file:
            shutil.move(file_path, corrected_2)
        elif '_unpaired' in spades_file:
            shutil.move(file_path, corrected_u)
    print('done')
    sys.stdout.flush()
    return (corrected_1, corrected_2, corrected_u)

def spades_assembly(read_files, outdir, log_file, kmers):
    '''
    This runs a SPAdes assembly, possibly continuing from a previous assembly.
    '''
    short1 = read_files[0]
    short2 = read_files[1]
    unpaired = read_files[2]
    kmer_string = ','.join([str(x) for x in kmers])
    this_kmer = 'k' + str(kmers[-1])
    message = 'SPAdes assembly: k=' + kmer_string
    log_file.write(message + '\n')
    dashes = '-' * len(message)
    log_file.write(dashes + '\n')
    command = [spades_path]
    command += ['--disable-rr', '-o', outdir, '-k', kmer_string]
    if len(kmers) > 1:
        last_kmer = 'k' + str(kmers[-2])
        command += ['--restart-from', last_kmer]
    else:
        command += ['--only-assembler', '-1', short1, '-2', short2]
        if unpaired:
            command += ['-s', unpaired]
    process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out, err = process.communicate()
    if len(err) > 0:
        print('\nSPAdes encountered an error:', file=sys.stderr)
        print(err, file=sys.stderr)
        quit()
    log_file.write(out)
    log_file.write('\n\n\n\n\n\n\n\n\n\n\n')
    graph_file = os.path.join(outdir, 'assembly_graph.fastg')
    paths_file = os.path.join(outdir, 'contigs.paths')

    parent_dir = os.path.dirname(outdir)
    moved_graph_file = os.path.join(parent_dir, this_kmer + '_assembly_graph.fastg')
    shutil.move(graph_file, moved_graph_file)
    moved_paths_file = os.path.join(parent_dir, this_kmer + '_contigs.paths')
    shutil.move(paths_file, moved_paths_file)


    sys.stdout.flush()
    return moved_graph_file

def strip_read_extensions(read_file_name):
    '''
    This function removes certain file extensions from a file name.
    '''
    no_extensions = read_file_name
    lower = no_extensions.lower()
    if lower.endswith('.gz'):
        no_extensions = no_extensions[:-3]
    lower = no_extensions.lower()
    if lower.endswith('.fasta') or lower.endswith('.fastq') or lower.endswith('.fastg'):
        no_extensions = no_extensions[:-6]
    return no_extensions

def make_log_file(outdir):
    '''
    Creates an empty log file named log.txt in the given directory.
    '''
    log_file_path = os.path.join(outdir, 'log.txt')
    open(log_file_path, 'w')
    return log_file_path

def get_kmer_range(reads_1_filename, reads_2_filename, log_file):
    '''
    Uses the read lengths to determine the k-mer range to be used in the SPAdes assembly.
    '''
    read_lengths = get_read_lengths(reads_1_filename) + get_read_lengths(reads_2_filename)
    read_lengths = sorted(read_lengths)
    median_read_length = read_lengths[len(read_lengths) // 2]
    starting_kmer = round_to_nearest_odd(starting_kmer_fraction * median_read_length)
    max_kmer = round_to_nearest_odd(max_kmer_fraction * median_read_length)
    if starting_kmer < 11:
        starting_kmer = 11
    if max_kmer > 127:
        max_kmer = 127
    interval = 2
    while True:
        kmer_range = range(starting_kmer, max_kmer, interval) + [max_kmer]
        if len(kmer_range) <= kmer_count:
            break
        interval += 2
    log_file.write('K-mer sizes for SPAdes assembly\n')
    log_file.write('-------------------------------\n')
    log_file.write('Median read length: ' + str(median_read_length) + '\n')
    log_file.write('Starting k-mer:     ' + str(starting_kmer) + '\n')
    log_file.write('Maximum k-mer:      ' + str(max_kmer) + '\n')
    log_file.write('k-mer range:        ' + ', '.join([str(x) for x in kmer_range]) + '\n')
    log_file.write('\n\n\n\n\n\n\n\n\n\n\n')
    return kmer_range


def round_to_nearest_odd(num):
    '''
    Rounds a float to an odd integer.
    '''
    round_up = int(math.ceil(num))
    if round_up % 2 == 0:
        round_up += 1
    round_down = int(math.floor(num))
    if round_down % 2 == 0:
        round_down -= 1
    up_diff = round_up - num
    down_diff = num - round_down
    if up_diff > down_diff:
        return round_down
    else:
        return round_up


def get_read_lengths(reads_filename):
    '''
    Returns a list of the read lengths for the given read file.
    '''
    # TO DO: check if .fastq.gz or just .fastq
    reads = gzip.open(reads_filename, 'rb')
    read_lengths = []
    i = 0
    for line in reads:
        if i == 1:
            read_lengths.append(len(line.strip()))
        i += 1
        if i == 4:
            i = 0
    return read_lengths


if __name__ == '__main__':
    main()
