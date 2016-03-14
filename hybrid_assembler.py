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
from assembly_graph import AssemblyGraph


spades_path = '/Users/Ryan/Applications/SPAdes-3.7.1-Darwin/bin/spades.py'
starting_kmer = 11
max_kmer = 99
kmer_interval = 2
read_depth_filter_cutoff = 0.25


def main():
    args = get_arguments()
    make_output_directory(args.out)
    log_filepath = make_log_file(args.out)
    log_file = open(log_filepath, 'w')
    check_file_exists(args.short1)
    check_file_exists(args.short2)
    # check_file_exists(args.long)
    assembly_graph = get_best_spades_graph(args.short1, args.short2, args.out, log_file)



def get_arguments():
    '''
    Parse the command line arguments.
    '''
    parser = argparse.ArgumentParser(description='Hybrid Assembler')
    parser.add_argument('short1', help='Illumina reads 1')
    parser.add_argument('short2', help='Illumina reads 2')
    parser.add_argument('long', help='Long reads')
    parser.add_argument('out', help='Output directory')
    return parser.parse_args()


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


def get_best_spades_graph(short1, short2, outdir, log_file):
    '''
    This function tries a SPAdes assembly at different kmers and returns the best.
    'The best' is defined as the smallest dead-end count after low-depth filtering.  If multiple
    graphs have the same dead-end count (e.g. zero!) then the highest kmer is used.
    '''
    spades_dir = os.path.join(outdir, '1_spades_assembly')
    read_correction_dir = os.path.join(spades_dir, 'read_correction')
    read_files = spades_read_correction(short1, short2, read_correction_dir, log_file)
    shutil.rmtree(read_correction_dir)
    assem_dir = os.path.join(spades_dir, 'assembly')
    print('Conducting SPAdes assemblies...\n')
    print('kmer\tsegments\tdead ends\tconnected_components')
    kmer_range = range(starting_kmer, max_kmer, kmer_interval) + [max_kmer]
    for i, kmer in enumerate(kmer_range):
        graph_file = spades_assembly(read_files, assem_dir, log_file, kmer_range[:i+1])
        assembly_graph = AssemblyGraph(graph_file, kmer)
        filter_graph(assembly_graph, spades_dir, kmer)
        dead_ends, connected_components, segment_count = get_graph_info(assembly_graph)
        
        print(str(kmer) + '\t' + str(segment_count) + '\t' + str(dead_ends) + '\t' + str(connected_components))
    shutil.rmtree(assem_dir)

def filter_graph(assembly_graph, graph_dir, kmer):
    '''
    For the given graph, this function normalises the read depths, filters based on depth and saves
    it, with '_filtered' in the filename.
    '''
    assembly_graph.normalise_read_depths()
    assembly_graph.filter_by_read_depth(read_depth_filter_cutoff)
    filtered_path = os.path.join(graph_dir, 'assembly_graph_k' + str(kmer) + '_filtered.fastg')
    assembly_graph.save_to_fastg(filtered_path)

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
    parent_dir = os.path.dirname(outdir)
    moved_graph_file = os.path.join(parent_dir, 'assembly_graph_' + this_kmer + '.fastg')
    shutil.move(graph_file, moved_graph_file)
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


if __name__ == '__main__':
    main()
