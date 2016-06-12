from __future__ import print_function
from __future__ import division
import sys

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



def create_spades_contig_bridges(graph, single_copy_segments, verbosity):
    '''
    Builds graph bridges using the SPAdes contig paths.
    '''
    if verbosity > 0:
        print()
        print('Bridging graph with SPAdes contig paths')
        print('---------------------------------------')
        sys.stdout.flush()

    bridge_path_set = set()
    single_copy_numbers = [x.number for x in single_copy_segments]
    for segment in single_copy_segments:
        for path in graph.paths.itervalues():
            flipped_path = [-x for x in reversed(path)]
            contig_bridges = find_contig_bridges(segment.number, path, single_copy_numbers)
            contig_bridges += find_contig_bridges(segment.number, flipped_path, single_copy_numbers)
            for contig_bridge in contig_bridges:
                flipped_contig_bridge = [-x for x in reversed(contig_bridge)]
                contig_bridge_str = ','.join([str(x) for x in contig_bridge])
                flipped_contig_bridge_str = ','.join([str(x) for x in flipped_contig_bridge])
                if contig_bridge_str not in bridge_path_set and \
                   flipped_contig_bridge_str not in bridge_path_set:
                    if contig_bridge[0] < 0 and contig_bridge[-1] < 0:
                        bridge_path_set.add(flipped_contig_bridge_str)
                    else:
                        bridge_path_set.add(contig_bridge_str)

    bridge_path_list = sorted(list([[int(y) for y in x.split(',')] for x in bridge_path_set]))

    # If multiple bridge paths start with or end with the same segment, that implies a conflict
    # between SPADes' paths and our single-copy determination. Throw these bridges out.
    bridge_paths_by_start = {}
    bridge_paths_by_end = {}
    for path in bridge_path_list:
        start = path[0]
        end = path[-1]
        if start not in bridge_paths_by_start:
            bridge_paths_by_start[start] = []
        if end not in bridge_paths_by_end:
            bridge_paths_by_end[end] = []
        if -end not in bridge_paths_by_start:
            bridge_paths_by_start[-end] = []
        if -start not in bridge_paths_by_end:
            bridge_paths_by_end[-start] = []
        bridge_paths_by_start[start].append(path)
        bridge_paths_by_end[end].append(path)
        bridge_paths_by_start[-end].append(path)
        bridge_paths_by_end[-start].append(path)
    conflicting_paths = []
    for grouped_paths in bridge_paths_by_start.itervalues():
        if len(grouped_paths) > 1:
            conflicting_paths += grouped_paths
    for grouped_paths in bridge_paths_by_end.itervalues():
        if len(grouped_paths) > 1:
            conflicting_paths += grouped_paths
    conflicting_paths_no_dups = []
    for path in conflicting_paths:
        if path not in conflicting_paths_no_dups:
            conflicting_paths_no_dups.append(path)
    conflicting_paths = conflicting_paths_no_dups
    if verbosity > 1:
        print('Bridge paths in conflict with single-copy segments: ', end='')
        if conflicting_paths:
            print(', '.join([str(x) for x in conflicting_paths]))
        else:
            print('none')
        print()

    final_bridge_paths = [x for x in bridge_path_list if x not in conflicting_paths]
    if verbosity > 1:
        print('Final SPAdes contig bridge paths: ', end='')
        if final_bridge_paths:
            print(', '.join([str(x) for x in final_bridge_paths]))
        else:
            print('none')
        print()

    return [Bridge(spades_contig_path=x, graph=graph) for x in final_bridge_paths]

def find_contig_bridges(segment_num, path, single_copy_numbers):
    '''
    This function returns a list of lists: every part of the path which starts on the segment_num
    and ends on any of the single_copy_numbers.
    '''
    bridge_paths = []
    indices = [i for i, x in enumerate(path) if abs(x) == segment_num]
    for index in indices:
        bridge_path = [path[index]]
        for i in range(index+1, len(path)):
            bridge_path.append(path[i])
            if path[i] in single_copy_numbers or -path[i] in single_copy_numbers:
                break
        else:
            bridge_path = []
        if bridge_path:
            bridge_paths.append(bridge_path)
    return bridge_paths




