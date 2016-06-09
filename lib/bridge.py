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
        self.bridge_type = ''

        if spades_contig_path and graph:
            self.graph_path = spades_contig_path
            self.start_segment = self.graph_path.pop(0)
            self.end_segment = self.graph_path.pop()
            self.bridge_sequence = graph.get_path_sequence(self.graph_path)
            self.bridge_type = 'spades_contig_bridge'


        # self.bridge_type = 'long_read_bridge_through_graph'
        # self.bridge_type = 'long_read_bridge_not_through_graph'


    def __repr__(self):
        return str(self.start_segment) + ' -> ' + ', '.join([str(x) for x in self.graph_path]) + \
               ' -> ' + str(self.end_segment)


    def get_quality(self):
        '''
        This function gives a quality score to bridges. This is used to sort bridges so we can apply
        the best ones first.
        '''
        # TO DO: Make this function more robust! When two bridges are of the same type, they should
        #        give different quality scores.

        if self.bridge_type == 'spades_contig_bridge':
            return 1.0

        if self.bridge_type == 'long_read_bridge_not_through_graph':
            return 2.0

        if self.bridge_type == 'long_read_bridge_through_graph':
            return 3.0

        return 0.0





















