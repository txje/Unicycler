import unittest
import os
import unicycler.assembly_graph


class TestAssemblyGraphFunctionsFastg(unittest.TestCase):
    """
    Tests various AssemblyGraph functions on a graph loaded from a SPAdes FASTG file.
    """

    def setUp(self):
        test_fastg = os.path.join(os.path.dirname(__file__), 'test.fastg')
        self.graph = unicycler.assembly_graph.AssemblyGraph(test_fastg, 25)

    def test_attributes(self):
        self.assertEqual(self.graph.overlap, 25)

    def test_segments(self):
        self.assertEqual(len(self.graph.segments), 336)

    def test_forward_links(self):
        link_count = 0
        for link_list in self.graph.forward_links.values():
            link_count += len(link_list)
        self.assertEqual(link_count, 904)

    def test_reverse_links(self):
        link_count = 0
        for link_list in self.graph.reverse_links.values():
            link_count += len(link_list)
        self.assertEqual(link_count, 904)

    def test_links_match(self):
        """
        Tests that each segment's forward links match up with its reverse complement's reverse
        links.
        """
        for segment, forward_link_list in self.graph.forward_links.items():
            reverse_link_list = self.graph.reverse_links[-segment]
            flipped_reverse_link_list = [-x for x in reverse_link_list]
            self.assertEqual(sorted(forward_link_list), sorted(flipped_reverse_link_list))

    def test_load_spades_paths(self):
        self.assertEqual(len(self.graph.paths), 0)
        paths_file = os.path.join(os.path.dirname(__file__), 'test.fastg.paths')
        self.graph.load_spades_paths(paths_file)
        self.assertEqual(len(self.graph.paths), 53)

    def test_get_median_read_depth(self):
        diff = abs(self.graph.get_median_read_depth() - 40.2)
        self.assertTrue(diff < 0.1)

    def test_normalise_read_depths(self):
        self.graph.normalise_read_depths()
        self.assertAlmostEqual(self.graph.get_median_read_depth(), 1.0)

    def test_get_total_length(self):
        self.assertEqual(self.graph.get_total_length(), 187896)

    def test_get_total_length_no_overlaps(self):
        self.assertEqual(self.graph.get_total_length_no_overlaps(), 179496)

    def test_total_dead_end_count(self):
        self.assertEqual(self.graph.total_dead_end_count(), 0)

    def test_dead_end_count(self):
        for seg_num in self.graph.segments:
            self.assertEqual(self.graph.dead_end_count(seg_num), 0)


# class TestAssemblyGraphFunctionsGfa(unittest.TestCase):
#     """
#     Tests various AssemblyGraph functions on a graph loaded from a GFA file.
#     """
#
#     def setUp(self):
#         pass
#
#     def test_segments(self):
#         pass
#
#     def test_forward_links(self):
#         pass
#
#     def test_reverse_links(self):
#         pass
