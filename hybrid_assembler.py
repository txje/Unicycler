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
import copy
from multiprocessing import cpu_count

SCIRPT_DIR = os.path.dirname(os.path.realpath(__file__))
sys.path.append(os.path.join(SCIRPT_DIR, 'lib'))
from assembly_graph import AssemblyGraph
from bridge import SpadesContigBridge, LongReadBridge, LoopUnrollingBridge, \
                   create_spades_contig_bridges, find_contig_bridges, \
                   create_long_read_bridges, create_loop_unrolling_bridges
from misc import int_to_str, float_to_str, quit_with_error, check_file_exists, check_graphmap, \
                 get_percentile

sys.dont_write_bytecode = True
from semi_global_aligner import add_aligning_arguments, fix_up_arguments, \
                                semi_global_align_long_reads, Read, Reference, \
                                load_references, load_long_reads, AlignmentScoringScheme, \
                                load_sam_alignments

# TO DO: make these parameters, not globals.
spades_path = '/Users/Ryan/Applications/SPAdes-3.8.0-Darwin/bin/spades.py'
graphmap_path = '/Users/Ryan/Applications/graphmap/bin/Mac/graphmap'
starting_kmer_fraction = 0.2 # Relative to median read length
max_kmer_fraction = 0.9 # Relative to median read length
kmer_count = 10
read_depth_filter_cutoff = 0.25

def main():
    '''
    Script execution starts here.
    '''
    full_command = ' '.join(sys.argv)
    args = get_arguments()

    verbosity = args.verbosity
    if verbosity > 0:
        print()

    check_file_exists(args.short1)
    check_file_exists(args.short2)
    check_file_exists(args.long)
    if not args.no_graphmap:
        check_graphmap(args.graphmap_path)

    make_output_directory(args.out, verbosity)
    unbridged_graph = os.path.join(args.out, '002_unbridged_graph.gfa')
    spades_bridged_graph_unmerged = os.path.join(args.out, '003_spades_bridges_unmerged.gfa')
    spades_bridged_graph_merged = os.path.join(args.out, '004_spades_bridges_merged.gfa')

    # Produce a SPAdes assembly graph with a k-mer that balances contig length and connectivity.
    assembly_graph = get_best_spades_graph(args.short1, args.short2, args.out,
                                           args.read_depth_filter, verbosity,
                                           args.threads, args.keep_temp)

    # Determine copy number and get single-copy segments.
    single_copy_segments = get_single_copy_segments(assembly_graph, verbosity, 0)
    assembly_graph.save_to_gfa(unbridged_graph, verbosity, save_copy_depth_info=True)

    # Make an initial set of bridges using the SPAdes contig paths.
    bridges = create_spades_contig_bridges(assembly_graph, single_copy_segments, verbosity)
    bridges += create_loop_unrolling_bridges(assembly_graph, single_copy_segments, verbosity)
    bridged_graph = copy.deepcopy(assembly_graph)
    bridged_graph.apply_bridges(bridges, verbosity, args.min_bridge_qual)
    bridged_graph.save_to_gfa(spades_bridged_graph_unmerged, verbosity, save_seg_type_info=True)
    bridged_graph.merge_all_possible()
    bridged_graph.save_to_gfa(spades_bridged_graph_merged, verbosity)
    if verbosity > 0:
        print()

    # Prepare the directory for long read alignment.
    alignment_dir = os.path.join(args.out, '005_read_alignment')
    if not os.path.exists(alignment_dir):
        os.makedirs(alignment_dir)
    graph_fasta = os.path.join(alignment_dir, 'graph_segments.fasta')
    alignments_sam = os.path.join(alignment_dir, 'long_read_alignments.sam')
    temp_alignment_dir = os.path.join(alignment_dir, 'temp')
    assembly_graph.save_to_fasta(graph_fasta)

    # If all long reads are available now, then we do the entire process in one pass.
    if args.long:
        references = load_references(graph_fasta, 0)
        reference_dict = {x.name: x for x in references}
        read_dict, read_names = load_long_reads(args.long, 1)
        scoring_scheme = AlignmentScoringScheme(args.scores)

        # Load existing alignments if available.
        if os.path.isfile(alignments_sam) and sam_references_match(alignments_sam, assembly_graph):
            print('SAM file already exists. Will use these alignments instead of conducting a new '
                  'alignment:')
            print('  ' + alignments_sam)
            alignments = load_sam_alignments(alignments_sam, read_dict, reference_dict,
                                             scoring_scheme)
            for alignment in alignments:
                read_dict[alignment.read.name].alignments.append(alignment)

        # Conduct alignment if existing alignments are not available.
        else:
            alignments_sam_in_progress = alignments_sam + '.incomplete'
            min_alignment_length = assembly_graph.overlap * 2 # TO DO: make this a parameter?
            allowed_overlap = int(round(assembly_graph.overlap * 1.1)) # TO DO: adjust?

            semi_global_align_long_reads(references, graph_fasta, read_dict, read_names,
                                         args.long, temp_alignment_dir, args.graphmap_path,
                                         args.threads, scoring_scheme, args.low_score,
                                         not args.no_graphmap, False, args.kmer,
                                         min_alignment_length, alignments_sam_in_progress,
                                         full_command, allowed_overlap, verbosity)
            shutil.move(alignments_sam_in_progress, alignments_sam)

        # Use the long reads which aligned entirely within contigs (which are most likely correct)
        # to determine a minimum score
        contained_reads = [x for x in read_dict.itervalues() if x.has_one_contained_alignment()]
        contained_scores = []
        for read in contained_reads:
            contained_scores += [x.scaled_score for x in read.alignments]
        min_scaled_score_percentile = 5.0 # TO DO: make this a parameter?
        min_scaled_score = get_percentile(contained_scores, min_scaled_score_percentile)

        if verbosity > 1:
            print()
            print('Setting the minimum scaled score to the ' +
                  float_to_str(min_scaled_score_percentile, 1) +
                  '% percentile of full read alignments:', float_to_str(min_scaled_score, 2))
            print()

        # Make the long read bridges and apply them - this is the good part!
        bridges += create_long_read_bridges(assembly_graph, read_dict, read_names,
                                            single_copy_segments, verbosity, bridges,
                                            min_scaled_score)




        quit() # TEMP




        bridged_graph = copy.deepcopy(assembly_graph)
        bridged_graph.apply_bridges(bridges, verbosity, args.min_bridge_qual)
        bridged_graph.save_to_gfa(os.path.join(args.out, '006_long_read_bridges_unmerged.gfa'),
                                  verbosity, save_seg_type_info=True)
        bridged_graph.merge_all_possible()
        bridged_graph.save_to_gfa(os.path.join(args.out, '007_long_read_bridges_merged.gfa'),
                                  verbosity)

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
            finished = True # TEMP

def get_arguments():
    '''
    Parse the command line arguments.
    '''
    parser = argparse.ArgumentParser(description='Hybrid Assembler')

    parser.add_argument('--short1', required=True, default=argparse.SUPPRESS,
                        help='FASTQ file of short reads (first reads in each pair).')
    parser.add_argument('--short2', required=True, default=argparse.SUPPRESS,
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
    parser.add_argument('--threads', type=int, required=False, default=argparse.SUPPRESS,
                        help='Number of CPU threads used to align (default: the number of '
                             'available CPUs)')
    parser.add_argument('--keep_temp', action='store_true', default=argparse.SUPPRESS,
                        help='Keep all temporary files in output directory (default: delete most '
                             'temporary files to save space)')
    parser.add_argument('--min_bridge_qual', type=float, default=0.0,
                        help='Bridges with a quality below this value will not be applied)')
    parser.add_argument('--verbosity', type=int, required=False, default=1,
                        help='Level of stdout information (0 to 2)')

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
    try:
        args.threads
    except AttributeError:
        args.threads = cpu_count()
        if args.verbosity > 2:
            print('\nThread count set to', args.threads)
    try:
        args.keep_temp
    except AttributeError:
        args.keep_temp = False

    if not args.long and not args.long_dir:
        quit_with_error('either --long or --long_dir is required')
    if args.long and args.long_dir:
        quit_with_error('--long and --long_dir cannot both be used')

    # Change the output directory to a more useful full path.
    args.out = os.path.abspath(args.out)

    return args

def make_output_directory(outdir, verbosity):
    '''
    Creates the output directory, if it doesn't already exist.  If it does exist, that's okay as
    long as it's empty, in which case this function does nothing.  If it exists and isn't empty,
    this function quits with an error.
    '''
    if verbosity > 1:
        print('Making output directory')
        print('-----------------------')
        print(outdir)
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    elif os.listdir(outdir) and verbosity > 1:
        print('The directory already exists and files may be reused and/or overwritten.')
    if verbosity > 1:
        print()

def get_best_spades_graph(short1, short2, outdir, read_depth_filter, verbosity, threads,
                          keep_temp):
    '''
    This function tries a SPAdes assembly at different kmers and returns the best.
    'The best' is defined as the smallest dead-end count after low-depth filtering.  If multiple
    graphs have the same dead-end count (e.g. zero!) then the highest kmer is used.
    '''
    spades_dir = os.path.join(outdir, '001_spades_assembly')
    reads = spades_read_correction(short1, short2, spades_dir, verbosity, threads)

    kmer_range = get_kmer_range(short1, short2, spades_dir, verbosity)

    assem_dir = os.path.join(spades_dir, 'assembly')

    if verbosity > 0:
        print('Conducting SPAdes assemblies')
        print('----------------------------')

    # Check to see if the SPAdes assembies already exist. If so, we'll use them instead of doing
    # the assembly again.
    files_exist = True
    file_list = []
    for kmer in kmer_range:
        graph_file = os.path.join(spades_dir, 'k' + str(kmer) + '_assembly_graph.gfa')
        file_list.append(graph_file)
        if not os.path.isfile(graph_file):
            files_exist = False
    if files_exist and verbosity > 0:
        print('Assemblies already exist. Will use these graphs instead of running SPAdes '
              'assemblies:')
        print(''.join(['  ' + x + '\n' for x in file_list]))

    # Conduct a SPAdes assembly for each k-mer (or load existing ones) and score them to choose
    # the best.
    if verbosity == 1:
        print('  k-mer   segments   dead ends         score')
    best_score = 0.0
    best_assembly_graph = None
    for i, kmer in enumerate(kmer_range):
        clean_graph_filename = os.path.join(spades_dir, 'k' + str(kmer) +
                                            '_assembly_graph.gfa')
        if files_exist:
            assembly_graph = AssemblyGraph(clean_graph_filename, kmer)
        else:
            graph_file, paths_file = spades_assembly(reads, assem_dir, kmer_range[:i+1],
                                                     verbosity, threads)
            assembly_graph = AssemblyGraph(graph_file, kmer, paths_file=paths_file)
            assembly_graph.clean(read_depth_filter)
            assembly_graph.save_to_gfa(os.path.join(spades_dir, clean_graph_filename), 0)
            if not keep_temp:
                os.remove(graph_file)
                os.remove(paths_file)
        dead_ends, connected_components, segment_count = get_graph_info(assembly_graph)
        if segment_count == 0:
            score = 0.0
        else:
            score = 1.0 / (segment_count * ((dead_ends + 1) ** 2))
        if verbosity == 1:
            print(int_to_str(kmer).rjust(7) + int_to_str(segment_count).rjust(11) +
                  int_to_str(dead_ends).rjust(12) + '{:.2e}'.format(score).rjust(14))
        if verbosity > 1:
            print('Summary for k' + int_to_str(kmer) + ':')
            print('segments:            ', int_to_str(segment_count))
            print('total length:        ', int_to_str(assembly_graph.get_total_length()), 'bp')
            print('dead ends:           ', int_to_str(dead_ends))
            print('connected components:', int_to_str(connected_components))
            print('score:               ', '{:.2e}'.format(score))
            print()
        if score >= best_score:
            best_kmer = kmer
            best_score = score
            best_assembly_graph = assembly_graph

    # Clean up SPAdes files.
    if not keep_temp and os.path.isdir(assem_dir):
        shutil.rmtree(assem_dir)

    if verbosity == 1:
        print()

    if best_score == 0.0:
        quit_with_error('none of the SPAdes assemblies produced assembled sequence')

    if verbosity > 0:
        print('Best kmer: ' + str(best_kmer))
        print()

    return best_assembly_graph

def get_graph_info(assembly_graph):
    '''
    This function returns some graph information that will be useful for determining which kmer is
    best.
    '''
    dead_ends = assembly_graph.total_dead_end_count()
    connected_components = len(assembly_graph.get_connected_components())
    segment_count = len(assembly_graph.segments)
    return dead_ends, connected_components, segment_count

def spades_read_correction(short1, short2, spades_dir, verbosity, threads):
    '''
    This runs SPAdes with the --only-error-correction option.
    '''
    if verbosity > 0:
        print('SPAdes read error correction')
        print('----------------------------')
        sys.stdout.flush()

    # If the corrected reads already exist, then we just use them and proceed.
    corrected_1 = os.path.join(spades_dir, 'corrected_1.fastq.gz')
    corrected_2 = os.path.join(spades_dir, 'corrected_2.fastq.gz')
    corrected_u = os.path.join(spades_dir, 'corrected_u.fastq.gz')
    if os.path.isfile(corrected_1) and os.path.isfile(corrected_2):
        if verbosity > 0:
            print('Corrected reads already exist. Will use these reads instead of running SPAdes '
                  'error correction:')
            print('  ' + corrected_1)
            print('  ' + corrected_2)
            if os.path.isfile(corrected_u):
                print('  ' + corrected_u)
            print()
        return (corrected_1, corrected_2, corrected_u)

    # If the corrected reads don't exist, then we run SPAdes in error correction only mode.
    read_correction_dir = os.path.join(spades_dir, 'read_correction')
    command = [spades_path, '-1', short1, '-2', short2, '-o', read_correction_dir,
               '--threads', str(threads), '--only-error-correction']
    process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    while process.poll() is None:
        spades_output = process.stdout.readline().rstrip()
        print_line = False
        if verbosity > 1:
            print_line = True
        elif verbosity > 0:
            if 'Command line:' in spades_output:
                spades_output = ' '.join(spades_output.split())
                print_line = True
        if spades_output and print_line:
            print(spades_output)
            sys.stdout.flush()

    spades_error = process.stderr.readline().strip()
    if spades_error:
        quit_with_error('SPAdes encountered an error: ' + spades_error)

    # Read error correction should be done now, so copy the correct read files to a more permanent
    # location.
    short1_no_extension = strip_read_extensions(short1)
    short2_no_extension = strip_read_extensions(short2)
    corrected_dir = os.path.join(read_correction_dir, 'corrected')
    files = os.listdir(corrected_dir)
    for spades_file in files:
        file_path = os.path.join(corrected_dir, spades_file)
        if short1_no_extension in spades_file:
            shutil.move(file_path, corrected_1)
        elif short2_no_extension in spades_file:
            shutil.move(file_path, corrected_2)
        elif '_unpaired' in spades_file:
            shutil.move(file_path, corrected_u)
    shutil.rmtree(read_correction_dir)

    if not os.path.isfile(corrected_1) or not os.path.isfile(corrected_2):
        quit_with_error('SPAdes read error correction failed')

    if verbosity > 0:
        print()
        print('Corrected reads:')
        print('  ' + corrected_1)
        print('  ' + corrected_2)
        if os.path.isfile(corrected_u):
            print('  ' + corrected_u)
        print()

    return (corrected_1, corrected_2, corrected_u)

def spades_assembly(read_files, outdir, kmers, verbosity, threads):
    '''
    This runs a SPAdes assembly, possibly continuing from a previous assembly.
    '''
    short1 = read_files[0]
    short2 = read_files[1]
    unpaired = read_files[2]
    kmer_string = ','.join([str(x) for x in kmers])
    this_kmer = 'k' + str(kmers[-1])
    command = [spades_path, '-o', outdir, '-k', kmer_string, '--threads', str(threads)]
    if len(kmers) > 1:
        last_kmer = 'k' + str(kmers[-2])
        command += ['--restart-from', last_kmer]
    else:
        command += ['--only-assembler', '-1', short1, '-2', short2]
        if unpaired:
            command += ['-s', unpaired]
    process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    while process.poll() is None:
        spades_output = process.stdout.readline().rstrip()
        if spades_output and verbosity > 1:
            print(spades_output)
            sys.stdout.flush()

    spades_error = process.stderr.readline().strip()
    if spades_error:
        quit_with_error('SPAdes encountered an error: ' + spades_error)

    graph_file = os.path.join(outdir, 'assembly_graph.fastg')
    paths_file = os.path.join(outdir, 'contigs.paths')

    parent_dir = os.path.dirname(outdir)
    moved_graph_file = os.path.join(parent_dir, this_kmer + '_assembly_graph.fastg')
    shutil.move(graph_file, moved_graph_file)
    moved_paths_file = os.path.join(parent_dir, this_kmer + '_contigs.paths')
    shutil.move(paths_file, moved_paths_file)

    return moved_graph_file, moved_paths_file

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

def get_kmer_range(reads_1_filename, reads_2_filename, spades_dir, verbosity):
    '''
    Uses the read lengths to determine the k-mer range to be used in the SPAdes assembly.
    '''
    if verbosity > 0:
        print('Choosing k-mer range for assembly')
        print('---------------------------------')
        sys.stdout.flush()

    # If the k-mer range file already exists, we use its values and proceed.
    kmer_range_filename = os.path.join(spades_dir, 'kmer_range')
    if os.path.isfile(kmer_range_filename):
        with open(kmer_range_filename, 'r') as kmer_range_file:
            kmer_range = kmer_range_file.readline().strip().split(', ')
        if len(kmer_range) == kmer_count:
            try:
                kmer_range = [int(x) for x in kmer_range]
                if verbosity > 0:
                    print('K-mer range already exists:')
                    print('  ' + kmer_range_filename)
                    print('Will use this existing range:')
                    print('  ' + ', '.join([str(x) for x in kmer_range]))
                    print()
                return kmer_range
            except ValueError:
                pass

    # If the code got here, then the k-mer range doesn't already exist and we'll create one by
    # examining the read lengths.
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
    kmer_range_str = ', '.join([str(x) for x in kmer_range])

    if verbosity > 0:
        print('Median read length: ' + str(median_read_length))
        print('Starting k-mer:     ' + str(starting_kmer))
        print('Maximum k-mer:      ' + str(max_kmer))
        print('k-mer range:        ' + kmer_range_str)
        print()

    kmer_range_file = open(kmer_range_filename, 'w')
    kmer_range_file.write(kmer_range_str)
    kmer_range_file.close()
    return kmer_range

def get_single_copy_segments(graph, verbosity, min_single_copy_length):
    '''
    Returns a list of the graph segments determined to be single-copy.
    '''
    if verbosity > 0:
        print('Finding single-copy segments')
        print('----------------------------')
        sys.stdout.flush()

    graph.determine_copy_depth(verbosity)
    single_copy_segments = [x for x in graph.get_single_copy_segments() \
                            if x.get_length() >= min_single_copy_length]

    if verbosity > 1:
        print()
    if verbosity > 0:
        total_single_copy_length = sum([x.get_length() for x in single_copy_segments])
        print(int_to_str(len(single_copy_segments)),
              'single copy segments (' + int_to_str(total_single_copy_length) + ' bp) out of',
              int_to_str(len(graph.segments)),
              'total segments (' + int_to_str(graph.get_total_length()) + ' bp)')
        sys.stdout.flush()

    return single_copy_segments

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

def sam_references_match(sam_filename, assembly_graph):
    '''
    Returns True if the references in the SAM header exactly match the graph segment numbers.
    '''
    sam_file = open(sam_filename, 'r')
    ref_numbers_in_sam = []
    for line in sam_file:
        if not line.startswith('@'):
            break
        if not line.startswith('@SQ'):
            continue
        line_parts = line.strip().split()
        if len(line_parts) < 2:
            continue
        ref_name_parts = line_parts[1].split(':')
        if len(ref_name_parts) < 2:
            continue
        try:
            ref_numbers_in_sam.append(int(ref_name_parts[1]))
        except ValueError:
            pass

    ref_numbers_in_sam = sorted(ref_numbers_in_sam)
    seg_numbers_in_graph = sorted(assembly_graph.segments.keys())
    return ref_numbers_in_sam == seg_numbers_in_graph

if __name__ == '__main__':
    main()
