#!/usr/bin/env python3
"""
Unicycler - bacterial genome assembler for hybrid read sets

Author: Ryan Wick
email: rrwick@gmail.com
"""

import argparse
import os
import sys
import shutil
import copy
import random
from multiprocessing import cpu_count
from .assembly_graph import AssemblyGraph
from .bridge import create_spades_contig_bridges, \
    create_long_read_bridges, create_loop_unrolling_bridges
from .misc import int_to_str, float_to_str, quit_with_error, get_percentile, \
    print_section_header, check_files_and_programs
from .spades_func import get_best_spades_graph
from .blast_func import find_start_gene, CannotFindStart
from .antpath import add_aligning_arguments, fix_up_arguments, semi_global_align_long_reads, \
    load_references, load_long_reads, AlignmentScoringScheme, load_sam_alignments, \
    print_alignment_summary_table
from .pilon_func import polish_with_pilon, CannotPolish


def main():
    """
    Script execution starts here.
    """
    # Fix the random seed so the program produces the same output every time it's run.
    random.seed(0)

    full_command = ' '.join(sys.argv)
    args = get_arguments()
    verbosity = args.verbosity
    files = [args.short1, args.short2]
    if not args.no_long:
        files.append(args.long)
    check_files_and_programs(files,
                             spades_path=args.spades_path,
                             graphmap_path=(None if (args.no_graphmap or args.no_long)
                                            else args.graphmap_path),
                             makeblastdb_path=(None if args.no_rotate else args.makeblastdb_path),
                             tblastn_path=(None if args.no_rotate else args.tblastn_path),
                             gene_db_path=(None if args.no_rotate else args.start_genes),
                             pilon_path=(None if args.no_pilon else args.pilon_path),
                             samtools_path=(None if args.no_pilon else args.samtools_path),
                             bowtie2_path=(None if args.no_pilon else args.bowtie2_path),
                             bowtie2_build_path=(None if args.no_pilon
                                                 else args.bowtie2_build_path))
    make_output_directory(args.out, verbosity)

    file_num = 1
    unbridged_graph_filename = os.path.join(args.out, str(file_num).zfill(3) +
                                            '_unbridged_graph.gfa')

    # Produce a SPAdes assembly graph with a k-mer that balances contig length and connectivity.
    if os.path.isfile(unbridged_graph_filename):
        if verbosity > 0:
            print()
            print('Unbridged graph already exists. Will use this graph instead of running SPAdes:')
            print('  ' + unbridged_graph_filename)
        unbridged_graph = AssemblyGraph(unbridged_graph_filename, None)
    else:
        unbridged_graph = get_best_spades_graph(args.short1, args.short2, args.out,
                                                args.read_depth_filter, verbosity, args.spades_path,
                                                args.threads, args.keep_temp, args.kmer_count,
                                                args.min_kmer_frac, args.max_kmer_frac,
                                                args.no_spades_correct)

    # Determine copy number and get single-copy segments.
    single_copy_segments = get_single_copy_segments(unbridged_graph, verbosity, 0)
    unbridged_graph.save_to_gfa(unbridged_graph_filename, verbosity, save_copy_depth_info=True)

    # Make an initial set of bridges using the SPAdes contig paths.
    print_section_header('Bridging graph with SPAdes contig paths', verbosity,
                         last_newline=(verbosity > 1))
    bridges = create_spades_contig_bridges(unbridged_graph, single_copy_segments)
    bridges += create_loop_unrolling_bridges(unbridged_graph)
    graph = copy.deepcopy(unbridged_graph)
    seg_nums_used_in_bridges = graph.apply_bridges(bridges, verbosity, args.min_bridge_qual)
    file_num += 1
    spades_bridged_graph_unmerged = os.path.join(args.out, str(file_num).zfill(3) +
                                                 '_spades_bridges_applied.gfa')
    graph.save_to_gfa(spades_bridged_graph_unmerged, verbosity, save_seg_type_info=True,
                      single_copy_segments=single_copy_segments)

    # Clean up unnecessary segments after bridging.
    graph.clean_up_after_bridging(single_copy_segments, seg_nums_used_in_bridges,
                                  args.min_component_size, args.min_dead_end_size, verbosity)
    file_num += 1
    spades_bridged_graph_cleaned = os.path.join(args.out, str(file_num).zfill(3) + '_cleaned.gfa')
    graph.save_to_gfa(spades_bridged_graph_cleaned, verbosity, save_seg_type_info=True,
                      single_copy_segments=single_copy_segments)

    # Merge the segments to simplify the graph.
    graph.merge_all_possible()
    file_num += 1
    spades_bridged_graph_merged = os.path.join(args.out, str(file_num).zfill(3) + '_merged.gfa')
    graph.save_to_gfa(spades_bridged_graph_merged, verbosity)

    # Prepare for long read alignment.
    alignment_dir = os.path.join(args.out, 'read_alignment_temp')
    graph_fasta = os.path.join(alignment_dir, 'all_segments.fasta')
    single_copy_segments_fasta = os.path.join(alignment_dir, 'single_copy_segments.fasta')
    alignments_sam = os.path.join(alignment_dir, 'long_read_alignments.sam')
    temp_alignment_dir = os.path.join(alignment_dir, 'temp')
    scoring_scheme = AlignmentScoringScheme(args.scores)
    min_alignment_length = unbridged_graph.overlap * 2  # TO DO: make this a parameter?
    if not args.no_long:
        if not os.path.exists(alignment_dir):
            os.makedirs(alignment_dir)
        unbridged_graph.save_to_fasta(graph_fasta)
        unbridged_graph.save_specific_segments_to_fasta(single_copy_segments_fasta,
                                                        single_copy_segments)

    # If all long reads are available now, then we do the entire process in one pass.
    if args.long:
        references = load_references(graph_fasta, 0)
        reference_dict = {x.name: x for x in references}
        read_dict, read_names = load_long_reads(args.long, verbosity)

        # Load existing alignments if available.
        if os.path.isfile(alignments_sam) and sam_references_match(alignments_sam, unbridged_graph):
            if verbosity > 0:
                print('\nSAM file already exists. Will use these alignments instead of conducting '
                      'a new alignment:')
                print('  ' + alignments_sam)
            alignments = load_sam_alignments(alignments_sam, read_dict, reference_dict,
                                             scoring_scheme, args.threads, verbosity)
            for alignment in alignments:
                read_dict[alignment.read.name].alignments.append(alignment)
            print_alignment_summary_table(read_dict, verbosity)

        # Conduct the alignment if an existing SAM is not available.
        else:
            alignments_1_sam = os.path.join(alignment_dir, 'long_read_alignments_pass_1.sam')
            alignments_1_in_progress = alignments_1_sam + '.incomplete'
            alignments_2_sam = os.path.join(alignment_dir, 'long_read_alignments_pass_2.sam')
            alignments_2_in_progress = alignments_2_sam + '.incomplete'

            allowed_overlap = int(round(unbridged_graph.overlap * 1.1))  # TO DO: adjust?
            low_score_threshold = [args.low_score]
            semi_global_align_long_reads(references, graph_fasta, read_dict, read_names,
                                         args.long, temp_alignment_dir, args.graphmap_path,
                                         args.threads, scoring_scheme, low_score_threshold,
                                         not args.no_graphmap, False, args.kmer,
                                         min_alignment_length, alignments_1_in_progress,
                                         full_command, allowed_overlap, False, verbosity,
                                         stdout_header='Aligning reads (first pass)')
            shutil.move(alignments_1_in_progress, alignments_1_sam)

            # Run partially aligned reads again using more sensitive settings.
            overlapping_read_names = [x.name for x in read_dict.values()
                                      if x.get_fraction_aligned() < 1.0]
            if overlapping_read_names:
                semi_global_align_long_reads(references, single_copy_segments_fasta, read_dict,
                                             overlapping_read_names, args.long, temp_alignment_dir,
                                             args.graphmap_path, args.threads, scoring_scheme,
                                             low_score_threshold, False, False, args.kmer,
                                             min_alignment_length, alignments_2_in_progress,
                                             full_command, allowed_overlap, True, verbosity,
                                             stdout_header='Aligning reads (second pass)',
                                             display_low_score=False)
                shutil.move(alignments_2_in_progress, alignments_2_sam)

                # Now we have to put together a final SAM file. If a read is in the second pass,
                # then we use the alignments from that SAM. Otherwise we take the alignments from
                # the first SAM.
                overlapping_read_names = set(overlapping_read_names)
                with open(alignments_sam, 'wt') as alignments_file:
                    with open(alignments_1_sam, 'rt') as alignments_1:
                        for line in alignments_1:
                            if line.startswith('@'):
                                alignments_file.write(line)
                            else:
                                read_name = line.split('\t', 1)[0]
                                if read_name not in overlapping_read_names:
                                    alignments_file.write(line)
                    with open(alignments_2_sam, 'rt') as alignments_2:
                        for line in alignments_2:
                            if not line.startswith('@'):
                                alignments_file.write(line)

            # If there are no overlapping alignments (that's unfortunate) we just rename the
            # first pass SAM to the final SAM.
            else:
                shutil.move(alignments_1_sam, alignments_sam)

            if args.keep_temp < 1:
                shutil.rmtree(alignment_dir)
            if args.keep_temp < 2 and os.path.isfile(alignments_1_sam):
                os.remove(alignments_1_sam)
            if args.keep_temp < 2 and os.path.isfile(alignments_2_sam):
                os.remove(alignments_2_sam)
            if args.keep_temp < 2 and os.path.isfile(graph_fasta):
                os.remove(graph_fasta)
            if args.keep_temp < 2 and os.path.isfile(single_copy_segments_fasta):
                os.remove(single_copy_segments_fasta)

        # Use the long reads which aligned entirely within contigs (which are most likely correct)
        # to determine a minimum score.
        contained_reads = [x for x in read_dict.values() if x.has_one_contained_alignment()]
        contained_scores = []
        for read in contained_reads:
            contained_scores += [x.scaled_score for x in read.alignments]
        min_scaled_score_percentile = 5.0  # TO DO: make this a parameter?
        min_scaled_score = get_percentile(contained_scores, min_scaled_score_percentile)

        if verbosity > 1:
            print('\nSetting the minimum scaled score to the ' +
                  float_to_str(min_scaled_score_percentile, 1) +
                  'th percentile of full read alignments:', float_to_str(min_scaled_score, 2))

        # Do the long read bridging - this is the good part!
        print_section_header('Building long read bridges', verbosity, last_newline=(verbosity == 1))
        bridges = create_long_read_bridges(unbridged_graph, read_dict, read_names,
                                           single_copy_segments, verbosity, bridges,
                                           min_scaled_score, args.threads, scoring_scheme,
                                           min_alignment_length)
        graph = copy.deepcopy(unbridged_graph)
        print_section_header('Bridging graph with long reads', verbosity,
                             last_newline=(verbosity > 1))
        seg_nums_used_in_bridges = graph.apply_bridges(bridges, verbosity, args.min_bridge_qual)
        file_num += 1
        read_bridged_graph_unmerged = os.path.join(args.out, str(file_num).zfill(3) +
                                                   '_long_read_bridges_applied.gfa')
        graph.save_to_gfa(read_bridged_graph_unmerged, verbosity, save_seg_type_info=True,
                          single_copy_segments=single_copy_segments)

        # Clean up unnecessary segments after bridging.
        graph.clean_up_after_bridging(single_copy_segments, seg_nums_used_in_bridges,
                                      args.min_component_size, args.min_dead_end_size, verbosity)
        file_num += 1
        read_bridged_graph_cleaned = os.path.join(args.out,
                                                  str(file_num).zfill(3) + '_cleaned.gfa')
        graph.save_to_gfa(read_bridged_graph_cleaned, verbosity, save_seg_type_info=True,
                          single_copy_segments=single_copy_segments)

        # Merge the segments to simplify the graph.
        graph.merge_all_possible()
        file_num += 1
        read_bridged_graph_merged = os.path.join(args.out, str(file_num).zfill(3) + '_merged.gfa')
        graph.save_to_gfa(read_bridged_graph_merged, verbosity)

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
            finished = True  # TEMP

    # Perform a final clean on the graph, including overlap removal.
    graph.final_clean(verbosity)
    if verbosity > 0:
        print_section_header('Bridged assembly graph', verbosity)
        print(graph.get_summary(verbosity, show_components=True), end='')
    file_num += 1
    cleaned_graph = os.path.join(args.out, str(file_num).zfill(3) + '_cleaned.gfa')
    graph.save_to_gfa(cleaned_graph, verbosity)

    # Rotate completed replicons in the graph to a standard starting gene.
    completed_replicons = graph.completed_circular_replicons()
    if not args.no_rotate and len(completed_replicons) > 0:
        print_section_header('Rotating completed replicons', verbosity, last_newline=False)
        blast_dir = os.path.join(args.out, 'blast_temp')
        if not os.path.exists(blast_dir):
            os.makedirs(blast_dir)
        completed_replicons = sorted(completed_replicons, reverse=True,
                                     key=lambda x: graph.segments[x].get_length())
        rotation_count = 0
        for i, completed_replicon in enumerate(completed_replicons):
            segment = graph.segments[completed_replicon]
            sequence = segment.forward_sequence
            if graph.overlap > 0:
                sequence = sequence[:-graph.overlap]
            depth = segment.depth
            if verbosity > 0:
                print('\nReplicon ' + str(i+1) + ' (' + int_to_str(len(sequence)), 'bp, ' +
                      float_to_str(depth, 2) + 'x):')
            try:
                start_gene, start_pos, flip = find_start_gene(sequence, args.start_genes,
                                                              args.start_gene_id,
                                                              args.start_gene_cov, blast_dir,
                                                              args.makeblastdb_path,
                                                              args.tblastn_path, args.threads,
                                                              verbosity)
            except CannotFindStart:
                if verbosity > 0:
                    print('  Unable to find a starting gene - not rotating sequence')
            else:
                print('  Starting gene:', start_gene)
                strand_str = '(reverse strand)' if flip else '(forward strand)'
                print('  Starting position:', int_to_str(start_pos), strand_str)
                segment.rotate_sequence(start_pos, flip, graph.overlap)
                rotation_count += 1
        if rotation_count:
            file_num += 1
            rotated_graph = os.path.join(args.out, str(file_num).zfill(3) + '_rotated.gfa')
            graph.save_to_gfa(rotated_graph, verbosity)
        if args.keep_temp < 2 and os.path.exists(blast_dir):
            shutil.rmtree(blast_dir)

    # Polish the final assembly!
    if not args.no_pilon:
        print_section_header('Polishing assembly with Pilon', verbosity)
        polish_dir = os.path.join(args.out, 'polish_temp')
        if not os.path.exists(polish_dir):
            os.makedirs(polish_dir)
        try:
            polish_with_pilon(graph, args.bowtie2_path, args.bowtie2_build_path, args.pilon_path,
                              args.samtools_path, args.min_polish_size, polish_dir, verbosity,
                              args.short1, args.short2, args.threads)
        except CannotPolish as e:
            if verbosity > 0:
                print('Unable to polish assembly using Pilon: ', end='')
                print(e.message)
        else:
            file_num += 1
            polished_graph = os.path.join(args.out, str(file_num).zfill(3) + '_polished.gfa')
            graph.save_to_gfa(polished_graph, verbosity)
        if args.keep_temp < 2 and os.path.exists(polish_dir):
            shutil.rmtree(polish_dir)

    # Save the final state as both a GFA and FASTA file.
    if verbosity > 0:
        print_section_header('Complete', verbosity, last_newline=False)
    graph.save_to_gfa(os.path.join(args.out, 'assembly.gfa'), verbosity)
    graph.save_to_fasta(os.path.join(args.out, 'assembly.fasta'), verbosity)


def get_arguments():
    """
    Parse the command line arguments.
    """
    this_script_dir = os.path.dirname(os.path.realpath(__file__))

    parser = argparse.ArgumentParser(description='Hybrid Assembler')

    parser.add_argument('--short1', required=True, default=argparse.SUPPRESS,
                        help='FASTQ file of short reads (first reads in each pair).')
    parser.add_argument('--short2', required=True, default=argparse.SUPPRESS,
                        help='FASTQ file of short reads (second reads in each pair).')
    parser.add_argument('--long', required=False, default=argparse.SUPPRESS,
                        help='FASTQ or FASTA file of long reads, if all reads are available at '
                             'start.')
    parser.add_argument('--long_dir', required=False, default=argparse.SUPPRESS,
                        help='Directory where FASTQ or FASTA read files will be deposited.')
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
    parser.add_argument('--no_rotate', action='store_true',
                        help='Do not rotate completed replicons to start at a standard gene')
    parser.add_argument('--start_genes', type=str,
                        default=os.path.join(this_script_dir, 'gene_data', 'start_genes.fasta'),
                        help='FASTA file of genes for start point of rotated replicons')
    parser.add_argument('--start_gene_id', type=float, default=90.0,
                        help='The minimum required BLAST percent identity for a start gene search')
    parser.add_argument('--start_gene_cov', type=float, default=95.0,
                        help='The minimum required BLAST percent coverage for a start gene search')
    parser.add_argument('--makeblastdb_path', type=str, default='makeblastdb',
                        help='Path to the makeblastdb executable')
    parser.add_argument('--tblastn_path', type=str, default='tblastn',
                        help='Path to the tblastn executable')
    parser.add_argument('--no_pilon', action='store_true',
                        help='Do not use Pilon to polish the final assembly')
    parser.add_argument('--bowtie2_path', type=str, default='bowtie2',
                        help='Path to the bowtie2 executable')
    parser.add_argument('--bowtie2_build_path', type=str, default='bowtie2-build',
                        help='Path to the bowtie2_build executable')
    parser.add_argument('--samtools_path', type=str, default='samtools',
                        help='Path to the samtools executable')
    parser.add_argument('--pilon_path', type=str, default='pilon.jar',
                        help='Path to the executable Pilon Java archive file')
    parser.add_argument('--min_polish_size', type=int, default=10000,
                        help='Sequences shorter than this value will not be polished using Pilon')
    parser.add_argument('--no_long', action='store_true',
                        help='Do not use any long reads - assemble with short reads only')
    parser.add_argument('--no_spades_correct', action='store_true',
                        help='Skip SPAdes error correction step')

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

    # Make sure that only one of the three long read arguments was used.
    long_arg_count = 1 if args.long else 0
    long_arg_count += 1 if args.long_dir else 0
    long_arg_count += 1 if args.no_long else 0
    if long_arg_count == 0:
        quit_with_error('One of the following options is required: '
                        '--long, --long_dir or --no_long')
    if long_arg_count > 1:
        quit_with_error('Only one of the following options can be used: '
                        '--long, --long_dir or --no_long')

    # Change some arguments to full paths.
    args.out = os.path.abspath(args.out)
    args.short1 = os.path.abspath(args.short1)
    args.short2 = os.path.abspath(args.short2)
    if args.long:
        args.long = os.path.abspath(args.long)
    if args.long_dir:
        args.long_dir = os.path.abspath(args.long_dir)

    return args


def make_output_directory(out_dir, verbosity):
    """
    Creates the output directory, if it doesn't already exist.
    """
    if verbosity > 1:
        print_section_header('Making output directory', verbosity)
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
        print(out_dir)
    elif os.listdir(out_dir) and verbosity > 1:
        print('The directory already exists and files may be reused and/or overwritten:')
        print('  ' + out_dir)


def get_single_copy_segments(graph, verbosity, min_single_copy_length):
    """
    Returns a list of the graph segments determined to be single-copy.
    """
    print_section_header('Finding single-copy segments', verbosity)
    graph.determine_copy_depth(verbosity)
    single_copy_segments = [x for x in graph.get_single_copy_segments()
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
    """
    Returns True if the references in the SAM header exactly match the graph segment numbers.
    """
    sam_file = open(sam_filename, 'rt')
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
