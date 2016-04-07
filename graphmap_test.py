from __future__ import print_function
from __future__ import division
import os
from assembly_graph import AssemblyGraph
from assembly_graph import Segment
from long_read_alignment import Alignment
from long_read_alignment import run_graphmap_alignment_all_segments
from long_read_alignment import run_graphmap_alignment_one_segment_at_a_time
from long_read_alignment import run_blasr_alignment_all_segments
from long_read_alignment import load_alignments
from long_read_alignment import load_long_reads
from long_read_alignment import filter_by_end_gaps


graph_filename = '/Users/Ryan/Dropbox/Uni_research/Projects/Hybrid_assembler/synthetic_test_case/assem_out/assembly_graph_k69_filtered.fastg'
pacbio_reads_fastq = '/Users/Ryan/Dropbox/Uni_research/Projects/Hybrid_assembler/synthetic_test_case/synthetic_pacbio_reads.fastq'
pacbio_reads_fasta = '/Users/Ryan/Dropbox/Uni_research/Projects/Hybrid_assembler/synthetic_test_case/synthetic_pacbio_reads.fasta'
graphmap_path = '/Users/Ryan/Applications/graphmap-semi-global/bin/Mac/graphmap'
blasr_path = '/Users/Ryan/Applications/blasr_install/blasr/blasr'


graph = AssemblyGraph(graph_filename, 69)
graph.repair_four_way_junctions()
graph.filter_by_read_depth(0.25)
graph.filter_homopolymer_loops()
graph.normalise_read_depths()

temp_dir = '/Users/Ryan/Desktop/alignment_test'
if not os.path.exists(temp_dir):
    os.makedirs(temp_dir)

sam_file = os.path.join(temp_dir, 'alignments.sam')
# run_graphmap_alignment_all_segments(graph, pacbio_reads_fastq, sam_file, graphmap_path, temp_dir)
# run_graphmap_alignment_one_segment_at_a_time(graph, pacbio_reads_fastq, sam_file, graphmap_path, temp_dir)
# run_blasr_alignment_all_segments(graph, pacbio_reads_fasta, sam_file, blasr_path, temp_dir)

segment_lengths = {x.get_fastg_header(True): x.get_length() for x in graph.segments.values()}
all_alignments = load_alignments(sam_file, segment_lengths)


print('Alignments:', len(all_alignments))
long_reads = load_long_reads(pacbio_reads_fastq)
print('Reads:', len(long_reads))
good_alignments, bad_alignments = filter_by_end_gaps(all_alignments, 10)
print('Good alignments:', len(good_alignments))
print('Bad alignments:', len(bad_alignments))


for alignment in bad_alignments:
    print(alignment)
    print('  read start gap:   ', alignment.read_start_pos)
    print('  read end gap:     ', alignment.read_end_gap)
    print('  segment start gap:', alignment.segment_start_pos)
    print('  segment end gap:  ', alignment.segment_end_gap)
    print()


for alignment in good_alignments:
    long_reads[alignment.read_name].add_alignment(alignment)

# for read in long_reads.values():
#     if read.name == 'S1_7_11909':
#         print()
#         print(read.name + ' ('  + str(len(read.sequence)) + ' bp)')
#         print('----------------------------')
#         for alignment in read.alignments:
#             print(alignment)
#             print('  read start gap:   ', alignment.read_start_pos)
#             print('  read end gap:     ', alignment.read_end_gap)
#             print('  segment start gap:', alignment.segment_start_pos)
#             print('  segment end gap:  ', alignment.segment_end_gap)
#             print()
#         print()


# contained_alignments = segment_with_alignments.get_contained_alignments()

# for alignment in bad_alignments:
#     print('Read name:', alignment.read_name)
#     print('  Read gaps:   ', alignment.read_start_pos, alignment.read_end_gap)
#     print('  Segment gaps:', alignment.segment_start_pos, alignment.segment_end_gap)

# print()
# print(len(all_alignments))
# print(len(contained_alignments))

