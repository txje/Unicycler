from copy_depth_graph import CopyDepthGraph
graph = CopyDepthGraph('/Users/Ryan/Desktop/aligner_test/Klebs_INF008/out_dir/assembly_graph.fastg', 113, 0.2)
graph.repair_four_way_junctions()
graph.filter_by_read_depth(0.25)
graph.filter_homopolymer_loops()
graph.normalise_read_depths()


# print(graph.get_n_segment_length(99))
# quit()


graph.determine_copy_depth()
graph.save_to_gfa('/Users/Ryan/Desktop/test.gfa')
