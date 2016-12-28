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
import multiprocessing
from .assembly_graph import AssemblyGraph
from .bridge import create_spades_contig_bridges, create_long_read_bridges, \
    create_loop_unrolling_bridges
from .misc import int_to_str, float_to_str, quit_with_error, get_percentile, print_v, bold, \
    print_section_header, check_files_and_programs, MyHelpFormatter, print_table, get_ascii_art
from .spades_func import get_best_spades_graph
from .blast_func import find_start_gene, CannotFindStart
from .unicycler_align import add_aligning_arguments, fix_up_arguments, AlignmentScoringScheme, \
    semi_global_align_long_reads, load_references, load_long_reads, load_sam_alignments, \
    print_alignment_summary_table
from .pilon_func import polish_with_pilon, CannotPolish
from . import settings
from .version import __version__


def main():
    """
    Script execution starts here.
    """
    # Fix the random seed so the program produces the same output every time it's run.
    random.seed(0)

    full_command = ' '.join(sys.argv)
    args = get_arguments()
    verbosity = args.verbosity

    files = []
    if args.short1 and args.short2:
      files.extend([args.short1, args.short2])
    if args.unpaired:
        files.append(args.unpaired)
    if args.long:
        files.append(args.long)
    check_files_and_programs(files, args,
                             spades_path=args.spades_path,
                             makeblastdb_path=(None if args.no_rotate else args.makeblastdb_path),
                             tblastn_path=(None if args.no_rotate else args.tblastn_path),
                             gene_db_path=(None if args.no_rotate else args.start_genes),
                             pilon_path=(None if args.no_pilon else args.pilon_path),
                             java_path=(None if args.no_pilon else args.java_path),
                             samtools_path=(None if args.no_pilon else args.samtools_path),
                             bowtie2_path=(None if args.no_pilon else args.bowtie2_path),
                             bowtie2_build_path=(None if args.no_pilon
                                                 else args.bowtie2_build_path))

    print_intro_message(args, verbosity, full_command)
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
        unbridged_graph = get_best_spades_graph(args.short1, args.short2, args.unpaired,
                                                args.out, settings.READ_DEPTH_FILTER, verbosity,
                                                args.spades_path, args.threads, args.keep_temp,
                                                args.kmer_count, args.min_kmer_frac,
                                                args.max_kmer_frac, args.no_correct,
                                                args.expected_linear_seqs)

    # Determine copy number and get single-copy segments.
    single_copy_segments = get_single_copy_segments(unbridged_graph, verbosity, 0)
    unbridged_graph.save_to_gfa(unbridged_graph_filename, verbosity, save_copy_depth_info=True,
                                leading_newline=False)

    # Make an initial set of bridges using the SPAdes contig paths. This step is skipped when
    # using conservative bridging mode (in that case we don't trust SPAdes contig paths at all).
    if args.mode == 0:
        bridges = []
        graph = copy.deepcopy(unbridged_graph)
    else:
        print_section_header('Bridging graph with SPAdes contigs', verbosity)
        bridges = create_spades_contig_bridges(unbridged_graph, single_copy_segments)
        bridges += create_loop_unrolling_bridges(unbridged_graph)
        graph = copy.deepcopy(unbridged_graph)
        if not bridges:
            print_v('none found', verbosity, 1)
        else:
            seg_nums_used_in_bridges = graph.apply_bridges(bridges, verbosity, args.min_bridge_qual,
                                                           unbridged_graph)
            file_num += 1
            spades_bridged_graph_unmerged = os.path.join(args.out, str(file_num).zfill(3) +
                                                         '_spades_bridges_applied.gfa')
            graph.save_to_gfa(spades_bridged_graph_unmerged, verbosity, save_seg_type_info=True,
                              save_copy_depth_info=True)

            # Clean up unnecessary segments after bridging.
            graph.clean_up_after_bridging_1(single_copy_segments, seg_nums_used_in_bridges,
                                            verbosity)
            graph.clean_up_after_bridging_2(seg_nums_used_in_bridges, args.min_component_size,
                                            args.min_dead_end_size, verbosity, unbridged_graph,
                                            single_copy_segments)
            if args.keep_temp > 1:
                file_num += 1
                graph.save_to_gfa(os.path.join(args.out, str(file_num).zfill(3) +
                                               '_cleaned.gfa'), verbosity,
                                  save_seg_type_info=True, save_copy_depth_info=True)
            graph.merge_all_possible(single_copy_segments, args.mode)
            if args.keep_temp > 1:
                file_num += 1
                graph.save_to_gfa(os.path.join(args.out, str(file_num).zfill(3) + '_merged.gfa'),
                                  verbosity)

    # Prepare for long read alignment.
    alignment_dir = os.path.join(args.out, 'read_alignment_temp')
    graph_fasta = os.path.join(alignment_dir, 'all_segments.fasta')
    single_copy_segments_fasta = os.path.join(alignment_dir, 'single_copy_segments.fasta')
    single_copy_segment_names = set(str(x.number) for x in single_copy_segments)
    alignments_sam = os.path.join(alignment_dir, 'long_read_alignments.sam')
    scoring_scheme = AlignmentScoringScheme(args.scores)
    min_alignment_length = unbridged_graph.overlap * \
        settings.MIN_ALIGNMENT_LENGTH_RELATIVE_TO_GRAPH_OVERLAP
    if args.long:
        if not os.path.exists(alignment_dir):
            os.makedirs(alignment_dir)
        unbridged_graph.save_to_fasta(graph_fasta)
        unbridged_graph.save_specific_segments_to_fasta(single_copy_segments_fasta,
                                                        single_copy_segments)

    # If all long reads are available now, then we do the entire process in one pass.
    if args.long:
        references = load_references(graph_fasta, 0)
        reference_dict = {x.name: x for x in references}
        read_dict, read_names, long_read_filename = load_long_reads(args.long, verbosity)

        # Load existing alignments if available.
        if os.path.isfile(alignments_sam) and sam_references_match(alignments_sam, unbridged_graph):
            if verbosity > 0:
                print('\nSAM file already exists. Will use these alignments instead of conducting '
                      'a new alignment:')
                print('  ' + alignments_sam)
            alignments = load_sam_alignments(alignments_sam, read_dict, reference_dict,
                                             scoring_scheme, verbosity)
            for alignment in alignments:
                read_dict[alignment.read.name].alignments.append(alignment)
            print_alignment_summary_table(read_dict, verbosity, False)

        # Conduct the alignment if an existing SAM is not available.
        else:
            alignments_1_sam = os.path.join(alignment_dir, 'long_read_alignments_pass_1.sam')
            alignments_1_in_progress = alignments_1_sam + '.incomplete'
            alignments_2_sam = os.path.join(alignment_dir, 'long_read_alignments_pass_2.sam')
            alignments_2_in_progress = alignments_2_sam + '.incomplete'

            allowed_overlap = int(round(unbridged_graph.overlap *
                                        settings.ALLOWED_ALIGNMENT_OVERLAP))
            low_score_threshold = [args.low_score]
            semi_global_align_long_reads(references, graph_fasta, read_dict, read_names,
                                         long_read_filename,
                                         args.threads, scoring_scheme, low_score_threshold,
                                         False, min_alignment_length, alignments_1_in_progress,
                                         full_command, allowed_overlap, 0, args.contamination,
                                         verbosity, stdout_header='Aligning reads (first pass)',
                                         single_copy_segment_names=single_copy_segment_names)
            shutil.move(alignments_1_in_progress, alignments_1_sam)

            # Some reads are aligned again with a more sensitive mode: those with multiple
            # alignments or unaligned parts.
            retry_read_names = [x.name for x in read_dict.values() if
                                (x.get_fraction_aligned() < settings.MIN_READ_FRACTION_ALIGNED or
                                 len(x.alignments) > 1) and x.get_length() >= min_alignment_length]
            if retry_read_names:
                semi_global_align_long_reads(references, single_copy_segments_fasta, read_dict,
                                             retry_read_names, long_read_filename,
                                             args.threads, scoring_scheme, low_score_threshold,
                                             False, min_alignment_length, alignments_2_in_progress,
                                             full_command, allowed_overlap, 3,
                                             args.contamination, verbosity,
                                             stdout_header='Aligning reads (second pass)',
                                             display_low_score=False,
                                             single_copy_segment_names=single_copy_segment_names)
                shutil.move(alignments_2_in_progress, alignments_2_sam)

                # Now we have to put together a final SAM file. If a read is in the second pass,
                # then we use the alignments from that SAM. Otherwise we take the alignments from
                # the first SAM.
                retry_read_names = set(retry_read_names)
                with open(alignments_sam, 'wt') as alignments_file:
                    with open(alignments_1_sam, 'rt') as alignments_1:
                        for line in alignments_1:
                            if line.startswith('@'):
                                alignments_file.write(line)
                            else:
                                read_name = line.split('\t', 1)[0]
                                if read_name not in retry_read_names:
                                    alignments_file.write(line)
                    with open(alignments_2_sam, 'rt') as alignments_2:
                        for line in alignments_2:
                            if not line.startswith('@'):
                                alignments_file.write(line)

            # If there are no low fraction reads, we can just rename the first pass SAM to the
            # final SAM.
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

        # Discard any reads that mostly align to known contamination.
        if args.contamination:
            filtered_read_names = []
            filtered_read_dict = {}
            contaminant_read_count = 0
            for read_name in read_names:
                if read_dict[read_name].mostly_aligns_to_contamination():
                    contaminant_read_count += 1
                else:
                    filtered_read_names.append(read_name)
                    filtered_read_dict[read_name] = read_dict[read_name]
            read_names = filtered_read_names
            read_dict = filtered_read_dict
            if verbosity > 1:
                print('\nDiscarded', contaminant_read_count, 'reads as contamination')

        # Use the long reads which aligned entirely within contigs (which are most likely correct)
        # to determine a minimum score.
        contained_reads = [x for x in read_dict.values() if x.has_one_contained_alignment()]
        contained_scores = []
        for read in contained_reads:
            contained_scores += [x.scaled_score for x in read.alignments]
        min_scaled_score = get_percentile(contained_scores, settings.MIN_SCALED_SCORE_PERCENTILE)

        if verbosity > 1:
            print('Setting the minimum scaled score to the ' +
                  float_to_str(settings.MIN_SCALED_SCORE_PERCENTILE, 1) +
                  'th percentile of full read alignments:', float_to_str(min_scaled_score, 2))

        # Do the long read bridging - this is the good part!
        print_section_header('Building long read bridges', verbosity)
        expected_linear_seqs = args.expected_linear_seqs > 0
        bridges = create_long_read_bridges(unbridged_graph, read_dict, read_names,
                                           single_copy_segments, verbosity, bridges,
                                           min_scaled_score, args.threads, scoring_scheme,
                                           min_alignment_length, expected_linear_seqs,
                                           args.min_bridge_qual)
        graph = copy.deepcopy(unbridged_graph)
        print_section_header('Bridging graph with long reads', verbosity)
        seg_nums_used_in_bridges = graph.apply_bridges(bridges, verbosity, args.min_bridge_qual,
                                                       unbridged_graph)
        file_num += 1
        read_bridged_graph_unmerged = os.path.join(args.out, str(file_num).zfill(3) +
                                                   '_long_read_bridges_applied.gfa')
        graph.save_to_gfa(read_bridged_graph_unmerged, verbosity, save_seg_type_info=True,
                          save_copy_depth_info=True)

        # Clean up unnecessary segments after bridging.
        graph.clean_up_after_bridging_1(single_copy_segments, seg_nums_used_in_bridges, verbosity)
        graph.clean_up_after_bridging_2(seg_nums_used_in_bridges, args.min_component_size,
                                        args.min_dead_end_size, verbosity, unbridged_graph,
                                        single_copy_segments)
        if args.keep_temp > 1:
            file_num += 1
            graph.save_to_gfa(os.path.join(args.out, str(file_num).zfill(3) +
                                           '_cleaned.gfa'), verbosity,
                              save_seg_type_info=True, save_copy_depth_info=True)
        graph.merge_all_possible(single_copy_segments, args.mode)
        if args.keep_temp > 1:
            file_num += 1
            graph.save_to_gfa(os.path.join(args.out, str(file_num).zfill(3) + '_merged.gfa'),
                              verbosity)

    # # If we are getting long reads incrementally, then we do the process iteratively.
    # elif args.long_dir:
    #     finished = False
    #     while not finished:
    #         # TO DO: WAIT FOR NEW READS TO BECOME AVAILABLE
    #         # TO DO: ALIGN LONG READS TO GRAPH
    #         # TO DO: PRODUCE BRIDGES USING LONG READ ALIGNMENTS
    #         # TO DO: APPLY THE BRIDGES TO THE GRAPH
    #         # TO DO: SAVE THE RESULTS
    #         # TO DO: ASSESS WHETHER THE PROCESS IS COMPLETE
    #         finished = True  # TEMP

    # Perform a final clean on the graph, including overlap removal.
    graph.final_clean(verbosity)
    print_section_header('Bridged assembly graph', verbosity)
    if verbosity > 0:
        graph.print_component_table()
    file_num += 1
    cleaned_graph = os.path.join(args.out, str(file_num).zfill(3) + '_final_clean.gfa')
    graph.save_to_gfa(cleaned_graph, verbosity)

    # Rotate completed replicons in the graph to a standard starting gene.
    completed_replicons = graph.completed_circular_replicons()
    if not args.no_rotate and len(completed_replicons) > 0:
        print_section_header('Rotating completed replicons', verbosity)

        rotation_result_table = [['Replicon', 'Length', 'Depth', 'Starting gene', 'Position',
                                  'Strand', 'Identity', 'Coverage']]
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
            rotation_result_row = [str(i+1), int_to_str(len(sequence)),
                                   float_to_str(depth, 2) + 'x']
            try:
                blast_hit = find_start_gene(sequence, args.start_genes, args.start_gene_id,
                                            args.start_gene_cov, blast_dir, args.makeblastdb_path,
                                            args.tblastn_path, args.threads, verbosity)
            except CannotFindStart:
                rotation_result_row += ['none found', '', '', '', '']
            else:
                rotation_result_row += [blast_hit.qseqid, int_to_str(blast_hit.start_pos),
                                        'reverse' if blast_hit.flip else 'forward',
                                        '%.1f' % blast_hit.pident + '%',
                                        '%.1f' % blast_hit.query_cov + '%']
                segment.rotate_sequence(blast_hit.start_pos, blast_hit.flip, graph.overlap)
                rotation_count += 1
            rotation_result_table.append(rotation_result_row)

        if verbosity > 0:
            print_table(rotation_result_table, alignments='LRRLRLRR', indent=0,
                        sub_colour={'none found': 'red'})
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
                              args.java_path, args.samtools_path, args.min_polish_size, polish_dir,
                              verbosity, args.short1, args.short2, args.unpaired, args.threads)
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
        print_section_header('Complete', verbosity)
    graph.save_to_gfa(os.path.join(args.out, 'assembly.gfa'), verbosity, leading_newline=False)
    graph.save_to_fasta(os.path.join(args.out, 'assembly.fasta'), verbosity, leading_newline=False,
                        min_length=args.min_fasta_length)


def get_arguments():
    """
    Parse the command line arguments.
    """
    description = bold('Unicycler: a hybrid assembly pipeline for bacterial genomes')
    this_script_dir = os.path.dirname(os.path.realpath(__file__))

    if '--helpall' in sys.argv or '--allhelp' in sys.argv or '--all_help' in sys.argv:
        sys.argv.append('--help_all')
    show_all_args = '--help_all' in sys.argv

    # Show the ASCII art if the terminal is wide enough for it.
    terminal_width = shutil.get_terminal_size().columns
    if terminal_width >= 70:
        full_description = 'R|' + get_ascii_art() + '\n\n' + description
    else:
        full_description = description
    parser = argparse.ArgumentParser(description=full_description, formatter_class=MyHelpFormatter,
                                     add_help=False)

    # Help options
    help_group = parser.add_argument_group('Help')
    help_group.add_argument('-h', '--help', action='help',
                            help='Show this help message and exit')
    help_group.add_argument('--help_all', action='help',
                            help='Show a help message with all program options')
    help_group.add_argument('--version', action='version', version='Unicycler v' + __version__,
                            help="Show Unicycler's version number")

    # Short read input options
    input_group = parser.add_argument_group('Input')
    input_group.add_argument('-1', '--short1', required=False, default=argparse.SUPPRESS,
                             help='FASTQ file of short reads (first reads in each pair)')
    input_group.add_argument('-2', '--short2', required=False, default=argparse.SUPPRESS,
                             help='FASTQ file of short reads (second reads in each pair)')
    input_group.add_argument('-s', '--unpaired', required=False, default=argparse.SUPPRESS,
                             help='FASTQ file of unpaired short reads')

    # Long read input options
    input_group.add_argument('-l', '--long', required=False, default=argparse.SUPPRESS,
                             help='FASTQ or FASTA file of long reads, if all reads are available '
                                  'at start.')

    # Output options
    output_group = parser.add_argument_group('Output')
    output_group.add_argument('-o', '--out', required=True, default=argparse.SUPPRESS,
                              help='Output directory')
    output_group.add_argument('--verbosity', type=int, required=False, default=1,
                              help='R|Level of stdout information (0 to 3, default: 1)\n  '
                                   '0 = no stdout, 1 = basic progress indicators, '
                                   '2 = extra info, 3 = debugging info')
    output_group.add_argument('--min_fasta_length', type=int, required=False, default=1,
                              help='Exclude contigs from the FASTA file which are shorter than '
                                   'this length (default: 1)')
    output_group.add_argument('--keep_temp', type=int, default=0,
                              help='R|Level of file retention (0 to 2, default: 0)\n  '
                                   '0 = only keep files at main checkpoints, '
                                   '1 = keep some temp files including SAM, '
                                   '2 = keep all temp files')

    other_group = parser.add_argument_group('Other')
    other_group.add_argument('-t', '--threads', type=int, required=False, default=argparse.SUPPRESS,
                             help='Number of threads used (default: number of CPUs, up to ' +
                                  str(settings.MAX_AUTO_THREAD_COUNT) + ')')
    other_group.add_argument('--mode', choices=['conservative', 'normal', 'bold'], default='normal',
                             help='B|Bridging mode (default: normal)\n'
                                  '  conservative = smaller contigs, lowest misassembly rate\n'
                                  '  normal = moderate contig size and misassembly rate\n'
                                  '  bold = longest contigs, higher misassembly rate')
    other_group.add_argument('--min_bridge_qual', type=float, default=argparse.SUPPRESS,
                             help='R|Do not apply bridges with a quality below this value\n'
                                  '  conservative mode default: ' +
                                  str(settings.CONSERVATIVE_MIN_BRIDGE_QUAL) + '\n'
                                  '  normal mode default: ' +
                                  str(settings.NORMAL_MIN_BRIDGE_QUAL) + '\n'
                                  '  bold mode default: ' +
                                  str(settings.BOLD_MIN_BRIDGE_QUAL)
                                  if show_all_args else argparse.SUPPRESS)
    other_group.add_argument('--expected_linear_seqs', type=int, required=False, default=0,
                             help='The expected number of linear (i.e. non-circular) sequences in '
                                  'the underlying sequence')

    # SPAdes assembly options
    spades_group = parser.add_argument_group('SPAdes assembly',
                                             'These options control the short read SPAdes '
                                             'assembly at the beginning of the Unicycler pipeline.'
                                             if show_all_args else argparse.SUPPRESS)
    spades_group.add_argument('--spades_path', type=str, default='spades.py',
                              help='Path to the SPAdes executable'
                                   if show_all_args else argparse.SUPPRESS)
    spades_group.add_argument('--no_correct', action='store_true',
                              help='Skip SPAdes error correction step'
                                   if show_all_args else argparse.SUPPRESS)
    spades_group.add_argument('--min_kmer_frac', type=float, default=0.2,
                              help='Lowest k-mer size for SPAdes assembly, expressed as a '
                                   'fraction of the read length'
                                   if show_all_args else argparse.SUPPRESS)
    spades_group.add_argument('--max_kmer_frac', type=float, default=0.95,
                              help='Highest k-mer size for SPAdes assembly, expressed as a '
                                   'fraction of the read length'
                                   if show_all_args else argparse.SUPPRESS)
    spades_group.add_argument('--kmer_count', type=int, default=10,
                              help='Number of k-mer steps to use in SPAdes assembly'
                                   if show_all_args else argparse.SUPPRESS)

    # Rotation options
    rotation_group = parser.add_argument_group('Assembly rotation',
                                               'These options control the rotation of completed '
                                               'circular sequence near the end of the Unicycler '
                                               'pipeline.'
                                               if show_all_args else argparse.SUPPRESS)
    rotation_group.add_argument('--no_rotate', action='store_true',
                                help='Do not rotate completed replicons to start at a standard gene'
                                     if show_all_args else argparse.SUPPRESS)
    rotation_group.add_argument('--start_genes', type=str,
                                default=os.path.join(this_script_dir, 'gene_data',
                                                     'start_genes.fasta'),
                                help='FASTA file of genes for start point of rotated replicons'
                                     if show_all_args else argparse.SUPPRESS)
    rotation_group.add_argument('--start_gene_id', type=float, default=90.0,
                                help='The minimum required BLAST percent identity for a start gene '
                                     'search'
                                     if show_all_args else argparse.SUPPRESS)
    rotation_group.add_argument('--start_gene_cov', type=float, default=95.0,
                                help='The minimum required BLAST percent coverage for a start gene '
                                     'search'
                                     if show_all_args else argparse.SUPPRESS)
    rotation_group.add_argument('--makeblastdb_path', type=str, default='makeblastdb',
                                help='Path to the makeblastdb executable'
                                     if show_all_args else argparse.SUPPRESS)
    rotation_group.add_argument('--tblastn_path', type=str, default='tblastn',
                                help='Path to the tblastn executable'
                                     if show_all_args else argparse.SUPPRESS)

    # Polishing options
    polish_group = parser.add_argument_group('Pilon polishing',
                                             'These options control the final assembly polish '
                                             'using Pilon at the end of the Unicycler pipeline.'
                                             if show_all_args else argparse.SUPPRESS)
    polish_group.add_argument('--no_pilon', action='store_true',
                              help='Do not use Pilon to polish the final assembly'
                                   if show_all_args else argparse.SUPPRESS)
    polish_group.add_argument('--bowtie2_path', type=str, default='bowtie2',
                              help='Path to the bowtie2 executable'
                                   if show_all_args else argparse.SUPPRESS)
    polish_group.add_argument('--bowtie2_build_path', type=str, default='bowtie2-build',
                              help='Path to the bowtie2_build executable'
                                   if show_all_args else argparse.SUPPRESS)
    polish_group.add_argument('--samtools_path', type=str, default='samtools',
                              help='Path to the samtools executable'
                                   if show_all_args else argparse.SUPPRESS)
    polish_group.add_argument('--pilon_path', type=str, default='pilon',
                              help='Path to a Pilon executable or the Pilon Java archive file'
                                   if show_all_args else argparse.SUPPRESS)
    polish_group.add_argument('--java_path', type=str, default='java',
                              help='Path to the java executable'
                                   if show_all_args else argparse.SUPPRESS)
    polish_group.add_argument('--min_polish_size', type=int, default=10000,
                              help='Sequences shorter than this value will not be polished using '
                                   'Pilon'
                                   if show_all_args else argparse.SUPPRESS)

    # Graph cleaning options
    cleaning_group = parser.add_argument_group('Graph cleaning',
                                               'These options control the removal of small '
                                               'leftover sequences after bridging is complete.'
                                               if show_all_args else argparse.SUPPRESS)
    cleaning_group.add_argument('--min_component_size', type=int, default=1000,
                                help='Unbridged graph components smaller than this size will be '
                                     'removed from the final graph'
                                     if show_all_args else argparse.SUPPRESS)
    cleaning_group.add_argument('--min_dead_end_size', type=int, default=1000,
                                help='Graph dead ends smaller than this size will be removed from '
                                     'the final graph'
                                     if show_all_args else argparse.SUPPRESS)

    # Add the arguments for the aligner, but suppress the help text.
    align_group = parser.add_argument_group('Long read alignment',
                                            'These options control the alignment of long reads to '
                                            'the assembly graph.'
                                            if show_all_args else argparse.SUPPRESS)
    add_aligning_arguments(align_group, show_all_args)

    args = parser.parse_args()
    fix_up_arguments(args)

    if args.keep_temp < 0 or args.keep_temp > 2:
        quit_with_error('--keep_temp must be between 0 and 2 (inclusive)')

    try:
        args.short1
    except AttributeError:
        args.short1 = None
    try:
        args.short2
    except AttributeError:
        args.short2 = None
    try:
        args.unpaired
    except AttributeError:
        args.unpaired = None
    try:
        args.long
    except AttributeError:
        args.long = None
    try:
        args.threads
    except AttributeError:
        args.threads = min(multiprocessing.cpu_count(), settings.MAX_AUTO_THREAD_COUNT)
        if args.verbosity > 2:
            print('\nThread count set to', args.threads)
    try:
        args.keep_temp
    except AttributeError:
        args.keep_temp = False

    # Set up bridging mode related stuff.
    user_set_bridge_qual = hasattr(args, 'min_bridge_qual')
    bridging_mode_val = 1
    if args.mode == 'conservative':
        bridging_mode_val = 0
        if not user_set_bridge_qual:
            args.min_bridge_qual = settings.CONSERVATIVE_MIN_BRIDGE_QUAL
    elif args.mode == 'normal':
        bridging_mode_val = 1
        if not user_set_bridge_qual:
            args.min_bridge_qual = settings.NORMAL_MIN_BRIDGE_QUAL
    elif args.mode == 'bold':
        bridging_mode_val = 2
        if not user_set_bridge_qual:
            args.min_bridge_qual = settings.BOLD_MIN_BRIDGE_QUAL
    args.mode = bridging_mode_val

    # Change some arguments to full paths.
    args.out = os.path.abspath(args.out)
    if args.short1:
        args.short1 = os.path.abspath(args.short1)
    if args.short2:
        args.short2 = os.path.abspath(args.short2)
    if args.unpaired:
        args.unpaired = os.path.abspath(args.unpaired)
    if args.short1 is None and args.short2 is None and args.unpaired is None:
        quit_with_error('Some short reads, either --short1, --short2 or --unpaired must be included')
    if args.long:
        args.long = os.path.abspath(args.long)
    # if args.long_dir:
    #     args.long_dir = os.path.abspath(args.long_dir)

    return args


def make_output_directory(out_dir, verbosity):
    """
    Creates the output directory, if it doesn't already exist.
    """
    print_v('', verbosity, 1)
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
        print_v('Making output directory:', verbosity, 1)
        print_v('  ' + out_dir, verbosity, 1)
    elif os.listdir(out_dir) and verbosity > 1:
        print_v('The output directory already exists and files may be reused and/or overwritten:',
                verbosity, 1)
        print_v('  ' + out_dir, verbosity, 1)
    else:  # directory exists but is empty
        print_v('The output directory already exists:', verbosity, 1)
        print_v('  ' + out_dir, verbosity, 1)


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


def print_intro_message(args, verbosity, full_command):
    """
    Prints a message at the start of the program's execution.
    """
    if verbosity == 0:
        return

    print_section_header('Starting Unicycler', verbosity)
    print('Command: ' + bold(full_command))
    if verbosity > 1:
        print()
        if args.mode == 0:
            print('Bridging mode: conservative')
            if args.min_bridge_qual == settings.CONSERVATIVE_MIN_BRIDGE_QUAL:
                print('  using default conservative bridge quality cutoff: ', end='')
            else:
                print('  using user-specified bridge quality cutoff: ', end='')
        elif args.mode == 1:
            print('Bridging mode: normal')
            if args.min_bridge_qual == settings.NORMAL_MIN_BRIDGE_QUAL:
                print('  using default normal bridge quality cutoff: ', end='')
            else:
                print('  using user-specified bridge quality cutoff: ', end='')
        else:  # args.mode == 2
            print('Bridging mode: bold')
            if args.min_bridge_qual == settings.BOLD_MIN_BRIDGE_QUAL:
                print('  using default bold bridge quality cutoff: ', end='')
            else:
                print('  using user-specified bridge quality cutoff: ', end='')
        print(float_to_str(args.min_bridge_qual, 2))
