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
from .unicycler_align import add_aligning_arguments, fix_up_arguments, \
    semi_global_align_long_reads, \
    load_references, load_long_reads, AlignmentScoringScheme, load_sam_alignments, \
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
    print_intro_message(args, verbosity, full_command)

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
                                                settings.READ_DEPTH_FILTER, verbosity,
                                                args.spades_path, args.threads, args.keep_temp,
                                                args.kmer_count, args.min_kmer_frac,
                                                args.max_kmer_frac, args.no_spades_correct)

    # Determine copy number and get single-copy segments.
    single_copy_segments = get_single_copy_segments(unbridged_graph, verbosity, 0)
    unbridged_graph.save_to_gfa(unbridged_graph_filename, verbosity, save_copy_depth_info=True)

    # Make an initial set of bridges using the SPAdes contig paths. This step is skipped when
    # using conservative bridging mode (in that case we don't trust SPAdes contig paths at all).
    if args.mode == 0:
        bridges = []
        graph = copy.deepcopy(unbridged_graph)
    else:
        print_section_header('Bridging graph with SPAdes contig paths', verbosity,
                             last_newline=(verbosity > 1))
        bridges = create_spades_contig_bridges(unbridged_graph, single_copy_segments)
        bridges += create_loop_unrolling_bridges(unbridged_graph)
        graph = copy.deepcopy(unbridged_graph)
        seg_nums_used_in_bridges = graph.apply_bridges(bridges, verbosity, args.min_bridge_qual,
                                                       unbridged_graph)
        file_num += 1
        spades_bridged_graph_unmerged = os.path.join(args.out, str(file_num).zfill(3) +
                                                     '_spades_bridges_applied.gfa')
        graph.save_to_gfa(spades_bridged_graph_unmerged, verbosity, save_seg_type_info=True,
                          save_copy_depth_info=True)

        # Clean up unnecessary segments after bridging.
        graph.clean_up_after_bridging_1(single_copy_segments, seg_nums_used_in_bridges, verbosity)
        if args.keep_temp > 1:
            file_num += 1
            graph.save_to_gfa(os.path.join(args.out, str(file_num).zfill(3) +
                                           '_cleaned_first_pass.gfa'), verbosity,
                              save_seg_type_info=True, save_copy_depth_info=True)
        graph.clean_up_after_bridging_2(seg_nums_used_in_bridges, args.min_component_size,
                                        args.min_dead_end_size, verbosity, unbridged_graph,
                                        single_copy_segments)
        if args.keep_temp > 1:
            file_num += 1
            graph.save_to_gfa(os.path.join(args.out, str(file_num).zfill(3) +
                                           '_cleaned_second_pass.gfa'), verbosity,
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
    alignments_sam = os.path.join(alignment_dir, 'long_read_alignments.sam')
    temp_alignment_dir = os.path.join(alignment_dir, 'temp')
    scoring_scheme = AlignmentScoringScheme(args.scores)
    min_alignment_length = unbridged_graph.overlap * \
        settings.MIN_ALIGNMENT_LENGTH_RELATIVE_TO_GRAPH_OVERLAP
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
        read_dict, read_names, long_read_filename = load_long_reads(args.long, verbosity)

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

            allowed_overlap = int(round(unbridged_graph.overlap *
                                        settings.ALLOWED_ALIGNMENT_OVERLAP))
            low_score_threshold = [args.low_score]
            semi_global_align_long_reads(references, graph_fasta, read_dict, read_names,
                                         long_read_filename, temp_alignment_dir, args.graphmap_path,
                                         args.threads, scoring_scheme, low_score_threshold,
                                         not args.no_graphmap, False, args.kmer,
                                         min_alignment_length, alignments_1_in_progress,
                                         full_command, allowed_overlap, False, verbosity,
                                         stdout_header='Aligning reads (first pass)')
            shutil.move(alignments_1_in_progress, alignments_1_sam)

            # Reads with a lot of unaligned parts are tried again, this time on extra sensitive
            # mode.
            retry_read_names = [x.name for x in read_dict.values()
                                if x.get_fraction_aligned() < settings.MIN_READ_FRACTION_ALIGNED]
            if retry_read_names:
                semi_global_align_long_reads(references, single_copy_segments_fasta, read_dict,
                                             retry_read_names, long_read_filename,
                                             temp_alignment_dir, args.graphmap_path,
                                             args.threads, scoring_scheme, low_score_threshold,
                                             False, False, args.kmer, min_alignment_length,
                                             alignments_2_in_progress, full_command,
                                             allowed_overlap, True, verbosity,
                                             stdout_header='Aligning reads (second pass)',
                                             display_low_score=False)
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

        # Use the long reads which aligned entirely within contigs (which are most likely correct)
        # to determine a minimum score.
        contained_reads = [x for x in read_dict.values() if x.has_one_contained_alignment()]
        contained_scores = []
        for read in contained_reads:
            contained_scores += [x.scaled_score for x in read.alignments]
        min_scaled_score = get_percentile(contained_scores, settings.MIN_SCALED_SCORE_PERCENTILE)

        if verbosity > 1:
            print('\nSetting the minimum scaled score to the ' +
                  float_to_str(settings.MIN_SCALED_SCORE_PERCENTILE, 1) +
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
        seg_nums_used_in_bridges = graph.apply_bridges(bridges, verbosity, args.min_bridge_qual,
                                                       unbridged_graph)
        file_num += 1
        read_bridged_graph_unmerged = os.path.join(args.out, str(file_num).zfill(3) +
                                                   '_long_read_bridges_applied.gfa')
        graph.save_to_gfa(read_bridged_graph_unmerged, verbosity, save_seg_type_info=True,
                          save_copy_depth_info=True)

        # Clean up unnecessary segments after bridging.
        graph.clean_up_after_bridging_1(single_copy_segments, seg_nums_used_in_bridges, verbosity)
        if args.keep_temp > 1:
            file_num += 1
            graph.save_to_gfa(os.path.join(args.out, str(file_num).zfill(3) +
                                           '_cleaned_first_pass.gfa'), verbosity,
                              save_seg_type_info=True, save_copy_depth_info=True)
        graph.clean_up_after_bridging_2(seg_nums_used_in_bridges, args.min_component_size,
                                        args.min_dead_end_size, verbosity, unbridged_graph,
                                        single_copy_segments)
        if args.keep_temp > 1:
            file_num += 1
            graph.save_to_gfa(os.path.join(args.out, str(file_num).zfill(3) +
                                           '_cleaned_second_pass.gfa'), verbosity,
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
    if verbosity > 0:
        print_section_header('Bridged assembly graph', verbosity)
        print(graph.get_summary(verbosity, show_components=True), end='')
    file_num += 1
    cleaned_graph = os.path.join(args.out, str(file_num).zfill(3) + '_final_clean.gfa')
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


class MyHelpFormatter(argparse.RawDescriptionHelpFormatter):
    """
    This is a custom formatter class for argparse. It allows for some custom formatting,
    in particular for the help texts with multiple options (like bridging mode and verbosity level).
    http://stackoverflow.com/questions/3853722
    """
    def __init__(self, prog):
        terminal_width = shutil.get_terminal_size().columns
        max_help_position = min(max(24, terminal_width // 3), 40)
        super().__init__(prog, max_help_position=max_help_position)

    def _split_lines(self, text, width):
        if text.startswith('B|') or text.startswith('R|'):
            text_lines = text[2:].splitlines()
            wrapped_text_lines = []
            for line in text_lines:
                if len(line) <= width:
                    wrapped_text_lines.append(line)
                else:
                    wrap_column = 2

                    # The bridging mode help text should wrap each line around to the column of
                    # the equals sign.
                    if text.startswith('B|'):
                        line_parts = line.split()
                        wrap_column += line.find('=')
                        join = ''
                        current_line = '  ' + line_parts[0]

                    # The other multi-option help texts should wrap an entire option at a time.
                    else:  # text.startswith('R|')
                        line_parts = line.split(', ')
                        join = ','
                        current_line = line_parts[0]
                    for part in line_parts[1:]:
                        if len(current_line) + len(join) + 1 + len(part) <= width:
                            current_line += join + ' ' + part
                        else:
                            wrapped_text_lines.append(current_line + join)
                            current_line = ' ' * wrap_column + part
                    wrapped_text_lines.append(current_line)
            return wrapped_text_lines
        else:
            return argparse.HelpFormatter._split_lines(self, text, width)

    def _fill_text(self, text, width, indent):
        if text.startswith('R|'):
            return argparse.RawDescriptionHelpFormatter._fill_text(self, text[2:], width, indent)
        else:
            return argparse.HelpFormatter._fill_text(self, text, width, indent)


def get_arguments():
    """
    Parse the command line arguments.
    """
    description = 'Unicycler: a hybrid assembly pipeline for bacterial genomes'
    this_script_dir = os.path.dirname(os.path.realpath(__file__))

    if '--helpall' in sys.argv or '--allhelp' in sys.argv or '--all_help' in sys.argv:
        sys.argv.append('--help_all')
    show_all_args = '--help_all' in sys.argv

    # Show the ASCII art if the terminal is wide enough for it.
    terminal_width = shutil.get_terminal_size().columns
    os.environ['COLUMNS'] = str(terminal_width)
    ascii_art = ("R|       __\n"
                 "       \ \___\n"
                 "        \ ___\\\n"
                 "        //\n"
                 "   ____//      _    _         _                     _\n"
                 " //_  //\\\\    | |  | |       |_|                   | |\n"
                 "//  \//  \\\\   | |  | | _ __   _   ___  _   _   ___ | |  ___  _ __\n"
                 "||  (O)  ||   | |  | || '_ \ | | / __|| | | | / __|| | / _ \| '__|\n"
                 "\\\\    \_ //   | |__| || | | || || (__ | |_| || (__ | ||  __/| |\n"
                 " \\\\_____//     \____/ |_| |_||_| \___| \__, | \___||_| \___||_|\n"
                 "                                        __/ |\n"
                 "                                       |___/\n\n")

    if terminal_width >= 70:
        full_description = ascii_art + description
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
    input_group.add_argument('-1', '--short1', required=True, default=argparse.SUPPRESS,
                             help='FASTQ file of short reads (first reads in each pair).')
    input_group.add_argument('-2', '--short2', required=True, default=argparse.SUPPRESS,
                             help='FASTQ file of short reads (second reads in each pair).')

    # Long read input options
    input_group.add_argument('-l', '--long', required=False, default=argparse.SUPPRESS,
                             help='FASTQ or FASTA file of long reads, if all reads are available '
                                  'at start.')
    # input_group.add_argument('-d', '--long_dir', required=False, default=argparse.SUPPRESS,
    #                          help='Directory where FASTQ or FASTA read files will be deposited.')
    input_group.add_argument('--no_long', action='store_true',
                             help='Do not use any long reads (assemble with short reads only)')

    # Output options
    output_group = parser.add_argument_group('Output')
    output_group.add_argument('-o', '--out', required=True, default=argparse.SUPPRESS,
                              help='Output directory')
    output_group.add_argument('--verbosity', type=int, required=False, default=1,
                              help='R|Level of stdout information (0 to 3, default: 1)\n  '
                                   '0 = no stdout, 1 = basic progress indicators, '
                                   '2 = extra info, 3 = debugging info')
    output_group.add_argument('--keep_temp', type=int, default=1,
                              help='R|Level of file retention (0 to 2, default: 1)\n  '
                                   '0 = only keep files at main checkpoints, '
                                   '1 = keep some temp files including SAM, '
                                   '2 = keep all temp files')

    other_group = parser.add_argument_group('Other')
    other_group.add_argument('-t', '--threads', type=int, required=False, default=argparse.SUPPRESS,
                             help='Number of threads used to align (default: number of CPUs)')
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

    # SPAdes assembly options
    spades_group = parser.add_argument_group('SPAdes assembly',
                                             'These options control the short read SPAdes '
                                             'assembly at the beginning of the Unicycler pipeline.'
                                             if show_all_args else argparse.SUPPRESS)
    spades_group.add_argument('--spades_path', type=str, default='spades.py',
                              help='Path to the SPAdes executable'
                                   if show_all_args else argparse.SUPPRESS)
    spades_group.add_argument('--no_spades_correct', action='store_true',
                              help='Skip SPAdes error correction step'
                                   if show_all_args else argparse.SUPPRESS)
    spades_group.add_argument('--min_kmer_frac', type=float, default=0.2,
                              help='Lowest k-mer size for SPAdes assembly, expressed as a '
                                   'fraction of the read length'
                                   if show_all_args else argparse.SUPPRESS)
    spades_group.add_argument('--max_kmer_frac', type=float, default=0.9,
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
    polish_group.add_argument('--pilon_path', type=str, default='pilon.jar',
                              help='Path to the executable Pilon Java archive file'
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
                                            'the assembly graph using Graphmap and/or '
                                            'Unicycler-align'
                                            if show_all_args else argparse.SUPPRESS)
    add_aligning_arguments(align_group, show_all_args)

    args = parser.parse_args()
    fix_up_arguments(args)

    if args.keep_temp < 0 or args.keep_temp > 2:
        quit_with_error('--keep_temp must be between 0 and 2 (inclusive)')

    try:
        args.long
    except AttributeError:
        args.long = None
    # try:
    #     args.long_dir
    # except AttributeError:
    #     args.long_dir = None

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
    # long_arg_count += 1 if args.long_dir else 0
    long_arg_count += 1 if args.no_long else 0
    if long_arg_count == 0:
        quit_with_error('One of the following options is required: '
                        '--long or --no_long')
    if long_arg_count > 1:
        quit_with_error('Only one of the following options can be used: '
                        '--long or --no_long')

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
    args.short1 = os.path.abspath(args.short1)
    args.short2 = os.path.abspath(args.short2)
    if args.long:
        args.long = os.path.abspath(args.long)
    # if args.long_dir:
    #     args.long_dir = os.path.abspath(args.long_dir)

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
    else:  # directory exists but is empty
        print('The directory already exists:')
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


def print_intro_message(args, verbosity, full_command):
    """
    Prints a message at the start of the program's execution.
    """
    if verbosity == 0:
        return

    print_section_header('Starting Unicycler', verbosity)
    print('Command:', full_command)
    if verbosity > 1:
        print()
        if args.mode == 0:
            print('Bridging mode: conservative')
            if args.min_bridge_qual == settings.CONSERVATIVE_MIN_BRIDGE_QUAL:
                print('  Using default conservative bridge quality cutoff: ', end='')
            else:
                print('  Using user-specified bridge quality cutoff: ', end='')
        elif args.mode == 1:
            print('Bridging mode: normal')
            if args.min_bridge_qual == settings.CONSERVATIVE_MIN_BRIDGE_QUAL:
                print('  Using default normal bridge quality cutoff: ', end='')
            else:
                print('  Using user-specified bridge quality cutoff: ', end='')
        else:  # args.mode == 2
            print('Bridging mode: bold')
            if args.min_bridge_qual == settings.CONSERVATIVE_MIN_BRIDGE_QUAL:
                print('  Using default bold bridge quality cutoff: ', end='')
            else:
                print('  Using user-specified bridge quality cutoff: ', end='')
        print(float_to_str(args.min_bridge_qual, 2))
