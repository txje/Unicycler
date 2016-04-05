from __future__ import print_function
from __future__ import division
import os
from assembly_graph import AssemblyGraph
from assembly_graph import Segment
from long_read import SegmentWithAlignments
from long_read import Alignment


graph_filename = '/Users/Ryan/Dropbox/Uni_research/Projects/Hybrid_assembler/synthetic_test_case/assem_out/assembly_graph_k69_filtered.fastg'
pacbio_reads = '/Users/Ryan/Dropbox/Uni_research/Projects/Hybrid_assembler/synthetic_test_case/synthetic_pacbio_reads.fastq'
graphmap_path = '/Users/Ryan/Applications/graphmap/bin/Mac/graphmap'


graph = AssemblyGraph(graph_filename, 69)
graph.repair_four_way_junctions()
graph.filter_by_read_depth(0.25)
graph.filter_homopolymer_loops()
graph.normalise_read_depths()


longest_segment = sorted(graph.segments.values(), key=lambda x: x.get_length())[-1]

temp_dir = '/Users/Ryan/Desktop/graphmap_test'
if not os.path.exists(temp_dir):
    os.makedirs(temp_dir)

segment_with_alignments = SegmentWithAlignments(longest_segment, pacbio_reads, graphmap_path, temp_dir, 10)


all_alignments = segment_with_alignments.alignments
contained_alignments = segment_with_alignments.get_contained_alignments()

for alignment in all_alignments:
    print('Read name:', alignment.read_name)
    print('  Read length:', len(alignment.read_sequence),
          '  Read start fraction:   ', alignment.read_start_fraction,
          '  Read end fraction:     ', alignment.read_end_fraction)
    print('  Segment start fraction:', alignment.segment_start_fraction,
          '  Segment end fraction:  ', alignment.segment_end_fraction)

print()
print(len(all_alignments))
print(len(contained_alignments))

