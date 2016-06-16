'''
This module describes bridges - links between two single-copy segments in an assembly graph.
Bridges can come from multiple sources, so there are a few different classes which all share
the important members/methods (for duck-typing purposes).
'''

from __future__ import print_function
from __future__ import division
import sys
from misc import float_to_str

class SpadesContigBridge(object):
    '''
    This class describes a bridge created from the contigs.paths file made by SPAdes.
    '''
    def __init__(self, graph, spades_contig_path):

        # The two single-copy segments which are being bridged.
        self.start_segment = None
        self.end_segment = None

        # The path through the unbridged graph.
        self.graph_path = []

        # The bridge sequence, gotten from the graph path.
        self.bridge_sequence = ''

        # A score used to determine the order of bridge application.
        self.quality = 1.0


        self.graph_path = spades_contig_path
        self.start_segment = self.graph_path.pop(0)
        self.end_segment = self.graph_path.pop()
        self.bridge_sequence = graph.get_path_sequence(self.graph_path)

        # TO DO: QUALITY CALCULATION
        # TO DO: QUALITY CALCULATION
        # TO DO: QUALITY CALCULATION
        # TO DO: QUALITY CALCULATION
        # TO DO: QUALITY CALCULATION
        # TO DO: QUALITY CALCULATION

    def __repr__(self):
        return 'SPAdes contig bridge: ' + str(self.start_segment) + ' -> ' + \
               ', '.join([str(x) for x in self.graph_path]) + ' -> ' + str(self.end_segment) + \
               ' (quality = ' + float_to_str(self.quality, 2) + ')'


class LongReadBridge(object):
    '''
    This class describes a bridge created from long read alignments.
    '''
    def __init__(self, spades_contig_path=None, graph=None):

        # The two single-copy segments which are being bridged.
        self.start_segment = None
        self.end_segment = None

        # The individual reads contributing to the bridge. The sequences/qualities are not for the
        # entire read, just the part in the bridge. The lengths are stored separately because
        # negative bridge lengths are possible (when the start and end segment alignments overlap).
        # In these cases we don't have any read sequences but will still need to know the lengths
        # (which will be negative).
        self.read_names = []
        self.read_sequences = []
        self.read_qualities = []
        self.read_lengths = []

        # The bridge can also contain incomplete read sequences which don't bridge the entire span
        # between the start and end segments. These are still useful as they can contribute to the
        # consensus sequence.
        self.start_partial_read_names = []
        self.start_partial_sequences = []
        self.start_partial_qualities = []
        self.end_partial_read_names = []
        self.end_partial_sequences = []
        self.end_partial_qualities = []

        # The consensus of all read sequences. If there is only one read, then this is the same as
        # that read's sequence. If there are multiple reads, this is hopefully a more accurate
        # sequence. If this is a case where the start and end segments overlap, then there will not
        # be a consensus sequence and consensus_length will be the mean of read_lengths (a negative
        # value).
        self.consensus_sequence = ''
        self.consensus_length = 0.0

        # The path through the unbridged graph, if one was found.
        self.graph_path = []

        # The bridge sequence, gotten from the graph path if a good path was found. Otherwise it's
        # from the consensus sequence.
        self.bridge_sequence = ''

        # A score used to determine the order of bridge application.
        self.quality = 1.0


        # TO DO
        # TO DO
        # TO DO
        # TO DO
        # TO DO
        # TO DO
        # TO DO
        # TO DO
        # TO DO
        # TO DO
        # TO DO
        # TO DO
        # TO DO
        # TO DO
        # TO DO
        # TO DO
        # TO DO
        # TO DO
        # TO DO
        # TO DO
        # TO DO
        # TO DO
        # TO DO
        # TO DO
        # TO DO
        # TO DO
        # TO DO
        # TO DO
        # TO DO
        # TO DO
        # TO DO
        # TO DO
        # TO DO
        # TO DO
        # TO DO
        # TO DO
        # TO DO
        # TO DO
        # TO DO



    def __repr__(self):
        return 'long read bridge: ' + str(self.start_segment) + ' -> ' + \
               ', '.join([str(x) for x in self.graph_path]) + ' -> ' + str(self.end_segment) + \
               ' (quality = ' + float_to_str(self.quality, 2) + ')'




class LoopUnrollingBridge(object):
    '''
    This class describes a bridge created from unrolling an assembly graph loop.
    '''
    def __init__(self, graph, start, end, middle, repeat):
        '''
        This constructor assumes the the start, end, middle and repeat segments form a simple loop
        in the graph supported by either a SPAdes contig or a long read alignment. It will use
        segment depths to determine the loop count and score the bridge's quality.
        '''
        # The two single-copy segments which are being bridged.
        self.start_segment = start
        self.end_segment = end

        # The path through the unbridged graph.
        self.graph_path = []

        # The bridge sequence, gotten from the graph path.
        self.bridge_sequence = ''

        # A score used to determine the order of bridge application. This value can only decrease
        self.quality = 1.0


        # The start segment and end segment should agree in depth. If they don't, that's very bad,
        # so depth_disagreement is applied to quality twice (squared effect).
        start_depth = graph.segments[abs(start)].depth
        end_depth = graph.segments[abs(end)].depth
        depth_disagreement = min(start_depth, end_depth) / max(start_depth, end_depth)
        self.quality *= (depth_disagreement * depth_disagreement) # has squared effect on quality

        # The loop count as determined by the repeat segment should agree with the loop count as
        # determined by the middle segment.
        mean_start_end_depth = (start_depth + end_depth) / 2
        loop_count_by_middle = graph.segments[abs(middle)].depth / mean_start_end_depth
        loop_count_by_repeat = (graph.segments[abs(repeat)].depth - mean_start_end_depth) / \
                                   mean_start_end_depth
        count_disagreement = min(loop_count_by_middle, loop_count_by_repeat) / \
                             max(loop_count_by_middle, loop_count_by_repeat)
        self.quality *= count_disagreement

        # We'll use whichever segment is longer (repeat or middle) for our loop count going
        # forward.
        if graph.segments[abs(repeat)].get_length() > graph.segments[abs(middle)].get_length():
            loop_count_float = loop_count_by_repeat
        else:
            loop_count_float = loop_count_by_middle

        # If the average loop count is near a whole number, that's better. If it's near 0.5, that's
        # very bad!
        if loop_count_float < 1.0:
            loop_count = 1
            closeness_to_whole_num = loop_count_float
        else:
            loop_count = int(round(loop_count_float))
            fractional_part = loop_count_float % 1
            distance_from_whole_num = min(fractional_part, 1.0 - fractional_part)
            closeness_to_whole_num = 1.0 - (2.0 * distance_from_whole_num)
        self.quality *= closeness_to_whole_num

        # Finally, we reduce the quality for higher loop counts, as those are harder to call.
        loop_count_penalty = (1 / loop_count) ** 0.5
        self.quality *= loop_count_penalty

        self.graph_path = [repeat]
        for _ in range(loop_count):
            self.graph_path += [middle, repeat]
        self.bridge_sequence = graph.get_path_sequence(self.graph_path)

    def __repr__(self):
        return 'loop unrolling bridge: ' + str(self.start_segment) + ' -> ' + \
               ', '.join([str(x) for x in self.graph_path]) + ' -> ' + str(self.end_segment) + \
               ' (quality = ' + float_to_str(self.quality, 2) + ')'













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

    return [SpadesContigBridge(spades_contig_path=x, graph=graph) for x in final_bridge_paths]

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

def create_loop_unrolling_bridges(graph, single_copy_segments, verbosity):
    '''
    This function creates loop unrolling bridges using the information in SPAdes paths.
    '''
    bridges = []
    simple_loops = graph.find_all_simple_loops()

    # A simple loop can either be caused by a repeat in one sequence (probably more typical) or by
    # a separate circular sequence which has some common sequence (less typical, but still very
    # possible: plasmids). We only want to unroll the former group, so we look for cases where the
    # loop's start or end is in a SPAdes contig path along with the middle. That implies that they
    # are on the same piece of DNA and can be unrolled.
    for start, end, middle, repeat in simple_loops:
        for path in graph.paths.itervalues():
            joined = False
            flipped_path = [-x for x in reversed(path)]
            if (start in path and middle in path) or \
               (end in path and middle in path) or \
               (start in flipped_path and middle in flipped_path) or \
               (end in flipped_path and middle in flipped_path):
                joined = True
                break
        if not joined:
            continue

        # If the code got here, then things look good and we'll make a loop unrolling bridge!
        bridges.append(LoopUnrollingBridge(graph, start, end, middle, repeat))

    return bridges

def create_long_read_bridges(graph, reads, single_copy_segments, verbosity):
    '''
    Makes bridges between single-copy segments using the alignments in the long reads.
    '''
    # TO DO
    # TO DO
    # TO DO
    # TO DO
    # TO DO
    # TO DO
    # TO DO
    # TO DO
    # TO DO
    # TO DO
    # TO DO
    # TO DO
    # TO DO
    # TO DO
    # TO DO
    # TO DO
    # TO DO

    return [] # TEMP


