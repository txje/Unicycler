from __future__ import print_function
from __future__ import division
import os
from assembly_graph import AssemblyGraph
from assembly_graph import Segment
from long_read_alignment import Alignment
from long_read_alignment import semi_global_align_long_reads
from long_read_alignment import LongRead



def main():
    '''
    Script execution starts here.
    '''

    graph_filename = ('/Users/Ryan/Dropbox/Uni_research/Projects/Hybrid_assembler/'
                      'synthetic_test_case/assem_out/assembly_graph_k69_filtered.fastg')
    pacbio_reads_fastq = ('/Users/Ryan/Dropbox/Uni_research/Projects/Hybrid_assembler/'
                          'synthetic_test_case/synthetic_pacbio_reads_subsampled.fastq')
    graphmap_path = ('/Users/Ryan/Applications/graphmap-extanchorend/bin/Mac/graphmap')

    graph = AssemblyGraph(graph_filename, 69)
    graph.repair_four_way_junctions()
    graph.filter_by_read_depth(0.25)
    graph.filter_homopolymer_loops()
    graph.normalise_read_depths()

    temp_dir = '/Users/Ryan/Desktop/alignment_test'
    if not os.path.exists(temp_dir):
        os.makedirs(temp_dir)

    sam_raw = os.path.join(temp_dir, 'alignments.sam')
    sam_filtered = os.path.join(temp_dir, 'filtered_alignments.sam')
    graph_fasta = os.path.join(temp_dir, 'graph.fasta')
    graph.save_to_fasta(graph_fasta)
    long_reads = semi_global_align_long_reads(graph_fasta, pacbio_reads_fastq, sam_raw,
                                              sam_filtered, temp_dir, graphmap_path, False, 8)
    os.remove(graph_fasta)

    for read in long_reads.itervalues():
        for alignment in read.alignments:
            if alignment.percent_identity < 80.0:
                print(alignment)




if __name__ == '__main__':
    main()
