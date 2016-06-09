from __future__ import print_function
from __future__ import division

class Bridge(object):
    '''
    This class describes a bridge between two single-copy segments in an assembly graph.
    '''
    def __init__(self, spades_contig_path=None, graph=None):
        self.start_segment = None
        self.end_segment = None
        self.graph_path = []
        self.bridge_sequence = ''

        if spades_contig_path and graph:
            self.graph_path = spades_contig_path
            self.start_segment = self.graph_path.pop(0)
            self.end_segment = self.graph_path.pop()
            self.bridge_sequence = graph.get_path_sequence(self.graph_path)





















