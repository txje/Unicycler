from __future__ import print_function
from __future__ import division
import os
import time
from assembly_graph import AssemblyGraph
from assembly_graph import Segment
from long_read_alignment import Alignment
from long_read_alignment import run_graphmap_alignment_one_segment_at_a_time
from long_read_alignment import load_alignments
from long_read_alignment import load_long_reads



def main():

    graph_filename = ('/Users/Ryan/Dropbox/Uni_research/Projects/Hybrid_assembler/'
                      'synthetic_test_case/assem_out/assembly_graph_k69_filtered.fastg')
    pacbio_reads_fastq = ('/Users/Ryan/Dropbox/Uni_research/Projects/Hybrid_assembler/'
                          'synthetic_test_case/synthetic_pacbio_reads.fastq')
    graphmap_path = ('/Users/Ryan/Applications/graphmap-extanchorend/bin/Mac/graphmap')


    graph = AssemblyGraph(graph_filename, 69)
    graph.repair_four_way_junctions()
    graph.filter_by_read_depth(0.25)
    graph.filter_homopolymer_loops()
    graph.normalise_read_depths()

    temp_dir = '/Users/Ryan/Desktop/alignment_test'
    if not os.path.exists(temp_dir):
        os.makedirs(temp_dir)

    sam_file = os.path.join(temp_dir, 'alignments.sam')
    # run_graphmap_alignment_one_segment_at_a_time(graph, pacbio_reads_fastq, sam_file, graphmap_path, temp_dir)
    references = {x.get_fastg_header(True): x.forward_sequence for x in graph.segments.itervalues()}
    all_alignments = load_alignments(sam_file, references)




if __name__ == '__main__':
    main()
