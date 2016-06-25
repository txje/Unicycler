#!/usr/bin/env python3
'''
Unicycler - bacterial genome assembler for hybrid read sets

Author: Ryan Wick
email: rrwick@gmail.com
'''

import argparse
import os
import sys
import shutil
import copy
from multiprocessing import cpu_count
from .assembly_graph import AssemblyGraph
from .bridge import create_spades_contig_bridges, \
                    create_long_read_bridges, create_loop_unrolling_bridges
from .misc import int_to_str, float_to_str, quit_with_error, get_percentile, \
                  print_section_header, check_files_and_programs
from .spades_func import get_best_spades_graph
from .antpath import add_aligning_arguments, fix_up_arguments, semi_global_align_long_reads, \
                     load_references, load_long_reads, AlignmentScoringScheme, load_sam_alignments

def main():
    '''
    Script execution starts here.
    '''
    full_command = ' '.join(sys.argv)
    args = get_arguments()
    verbosity = args.verbosity
    check_files_and_programs([args.short1, args.short2, args.long], spades_path=args.spades_path,
                             graphmap_path=(None if args.no_graphmap else args.graphmap_path))
    make_output_directory(args.out, verbosity)

    unbridged_graph = os.path.join(args.out, '001_unbridged_graph.gfa')
    spades_bridged_graph_unmerged = os.path.join(args.out, '002_spades_bridges_applied.gfa')
    spades_bridged_graph_merged = os.path.join(args.out, '003_spades_bridges_merged.gfa')

    # Produce a SPAdes assembly graph with a k-mer that balances contig length and connectivity.
    if os.path.isfile(unbridged_graph):
        if verbosity > 0:
            print()
            print('Unbridged graph already exists. Will use this graph instead of running SPAdes:')
            print('  ' + unbridged_graph)
        assembly_graph = AssemblyGraph(unbridged_graph, None)
    else:
        assembly_graph = get_best_spades_graph(args.short1, args.short2, args.out,
                                               args.read_depth_filter, verbosity, args.spades_path,
                                               args.threads, args.keep_temp, args.kmer_count,
                                               args.min_kmer_frac, args.max_kmer_frac)

    # Determine copy number and get single-copy segments.
    single_copy_segments = get_single_copy_segments(assembly_graph, verbosity, 0)
    assembly_graph.save_to_gfa(unbridged_graph, verbosity, save_copy_depth_info=True)

    # Make an initial set of bridges using the SPAdes contig paths.
    print_section_header('Bridging graph with SPAdes contig paths', verbosity)
    bridges = create_spades_contig_bridges(assembly_graph, single_copy_segments, verbosity)
    bridges += create_loop_unrolling_bridges(assembly_graph, single_copy_segments, verbosity)
    bridged_graph = copy.deepcopy(assembly_graph)
    bridged_graph.apply_bridges(bridges, verbosity, args.min_bridge_qual, single_copy_segments)
    bridged_graph.save_to_gfa(spades_bridged_graph_unmerged, verbosity, save_seg_type_info=True,
                              single_copy_segments=single_copy_segments)
    bridged_graph.remove_small_components(args.min_component_size, verbosity)
    bridged_graph.merge_all_possible()
    bridged_graph.save_to_gfa(spades_bridged_graph_merged, verbosity)

    # Prepare the directory for long read alignment.
    alignment_dir = os.path.join(args.out, 'read_alignment_temp')
    if not os.path.exists(alignment_dir):
        os.makedirs(alignment_dir)
    graph_fasta = os.path.join(alignment_dir, 'graph_segments.fasta')
    alignments_sam = os.path.join(alignment_dir, 'long_read_alignments.sam')
    temp_alignment_dir = os.path.join(alignment_dir, 'temp')
    assembly_graph.save_to_fasta(graph_fasta)

    scoring_scheme = AlignmentScoringScheme(args.scores)
    min_alignment_length = assembly_graph.overlap * 2 # TO DO: make this a parameter?

    # If all long reads are available now, then we do the entire process in one pass.
    if args.long:
        references = load_references(graph_fasta, 0)
        reference_dict = {x.name: x for x in references}
        read_dict, read_names = load_long_reads(args.long, verbosity)

        # Load existing alignments if available.
        if os.path.isfile(alignments_sam) and sam_references_match(alignments_sam, assembly_graph):
            if verbosity > 0:
                print('\nSAM file already exists. Will use these alignments instead of conducting '
                      'a new alignment:')
                print('  ' + alignments_sam)
            alignments = load_sam_alignments(alignments_sam, read_dict, reference_dict,
                                             scoring_scheme)
            for alignment in alignments:
                read_dict[alignment.read.name].alignments.append(alignment)

        # Conduct alignment if existing alignments are not available.
        else:
            alignments_sam_in_progress = alignments_sam + '.incomplete'
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
        contained_reads = [x for x in read_dict.values() if x.has_one_contained_alignment()]
        contained_scores = []
        for read in contained_reads:
            contained_scores += [x.scaled_score for x in read.alignments]
        min_scaled_score_percentile = 5.0 # TO DO: make this a parameter?
        min_scaled_score = get_percentile(contained_scores, min_scaled_score_percentile)

        if verbosity > 1:
            print('\nSetting the minimum scaled score to the ' +
                  float_to_str(min_scaled_score_percentile, 1) +
                  '% percentile of full read alignments:', float_to_str(min_scaled_score, 2))

        # Do the long read bridging - this is the good part!
        read_bridged_graph_unmerged = os.path.join(args.out, '006_long_read_bridges_applied.gfa')
        read_bridged_graph_merged = os.path.join(args.out, '007_long_read_bridges_merged.gfa')
        print_section_header('Building long read bridges', verbosity)
        bridges = create_long_read_bridges(assembly_graph, read_dict, read_names,
                                           single_copy_segments, verbosity, bridges,
                                           min_scaled_score, args.threads, scoring_scheme,
                                           min_alignment_length)
        bridged_graph = copy.deepcopy(assembly_graph)
        print_section_header('Bridging graph with long reads', verbosity)
        bridged_graph.apply_bridges(bridges, verbosity, args.min_bridge_qual, single_copy_segments)
        bridged_graph.save_to_gfa(read_bridged_graph_unmerged, verbosity, save_seg_type_info=True,
                                  single_copy_segments=single_copy_segments)

        bridged_graph.remove_small_components(args.min_component_size, verbosity)
        bridged_graph.merge_all_possible()
        bridged_graph.remove_small_dead_ends(args.min_dead_end_size, verbosity)
        bridged_graph.merge_all_possible()

        bridged_graph.save_to_gfa(read_bridged_graph_merged, verbosity)
        if verbosity > 0:
            print()
            print(bridged_graph.get_summary('Final assembly graph'))

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
    parser.add_argument('--read_depth_filter', type=float, required=False, default=0.5,
                        help='Minimum allowed read depth, expressed as a fraction of the median'
                             'read depth. Graph segments with less depth will be removed.')
    parser.add_argument('--threads', type=int, required=False, default=argparse.SUPPRESS,
                        help='Number of CPU threads used to align (default: the number of '
                             'available CPUs)')
    parser.add_argument('--keep_temp', type=int, default=1,
                        help='0 = keep only main checkpoints, 1 = keep some temporary files, '
                             'including alignment SAM, 2 = keep all temporary files')
    parser.add_argument('--min_bridge_qual', type=float, default=5.0,
                        help='Bridges with a quality below this value will not be applied)')
    parser.add_argument('--min_component_size', type=int, default=1000,
                        help='Unbridged graph components smaller than this size will be removed '
                             'from the final graph')
    parser.add_argument('--min_dead_end_size', type=int, default=1000,
                        help='Graph dead ends smaller than this size will be removed from the '
                             'final graph')
    parser.add_argument('--spades_path', type=str, default='spades.py',
                        help='Path to the SPAdes executable')
    parser.add_argument('--min_kmer_frac', type=float, default=0.2,
                        help='Lowest k-mer size for SPAdes assembly, expressed as a fraction of '
                             'the read length')
    parser.add_argument('--max_kmer_frac', type=float, default=0.9,
                        help='Highest k-mer size for SPAdes assembly, expressed as a fraction of '
                             'the read length')
    parser.add_argument('--kmer_count', type=int, default=10,
                        help='Number of k-mer steps to use in SPAdes assembly')
    parser.add_argument('--verbosity', type=int, required=False, default=1,
                        help='Level of stdout information (0 to 2)')

    # Add the arguments for the aligner, but suppress the help text.
    add_aligning_arguments(parser, True)

    args = parser.parse_args()
    fix_up_arguments(args)

    if args.keep_temp < 0 or args.keep_temp > 2:
        quit_with_error('--keep_temp must be between 0 and 2 (inclusive)')

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
    Creates the output directory, if it doesn't already exist.
    '''
    if verbosity > 1:
        print_section_header('Making output directory', verbosity)
    if not os.path.exists(outdir):
        os.makedirs(outdir)
        print(outdir)
    elif os.listdir(outdir) and verbosity > 1:
        print('The directory already exists and files may be reused and/or overwritten:')
        print('  ' + outdir)

def get_single_copy_segments(graph, verbosity, min_single_copy_length):
    '''
    Returns a list of the graph segments determined to be single-copy.
    '''
    print_section_header('Finding single-copy segments', verbosity)
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
              'total segments (' + int_to_str(graph.get_total_length()) + ' bp)', flush=True)

    return single_copy_segments

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
