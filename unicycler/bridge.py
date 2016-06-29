'''
This module describes bridges - links between two single-copy segments in an assembly graph.
Bridges can come from multiple sources, so there are a few different classes which all share
the important members/methods (for duck-typing purposes).

Author: Ryan Wick
email: rrwick@gmail.com
'''

from multiprocessing.dummy import Pool as ThreadPool
import time
from .misc import int_to_str, float_to_str, reverse_complement, print_progress_line, \
                  weighted_average
from .cpp_function_wrappers import multiple_sequence_alignment, fully_global_alignment

class SpadesContigBridge(object):
    '''
    This class describes a bridge created from the contigs.paths file made by SPAdes.

    Quality is affected by:
      * How well the start and end segments' depths agree.
      * The depth consistency within the path (only applies to bridges where the path segments
        exclusively lead to the start/end segments).
    '''
    def __init__(self, graph, spades_contig_path):

        # The numbers of the two single-copy segments which are being bridged.
        self.start_segment = None
        self.end_segment = None

        # The path through the unbridged graph.
        self.graph_path = []

        # The bridge sequence, gotten from the graph path.
        self.bridge_sequence = ''

        # The bridge depth, a weighted mean of the start and end depths.
        self.depth = 0.0

        # A score used to determine the order of bridge application.
        self.quality = 40.0

        # The first and last values in spades_contig_path are the start and end segments. The
        # values in between are the path.
        self.graph_path = spades_contig_path
        self.start_segment = self.graph_path.pop(0)
        self.end_segment = self.graph_path.pop()
        self.bridge_sequence = graph.get_path_sequence(self.graph_path)
        start_seg = graph.segments[abs(self.start_segment)]
        end_seg = graph.segments[abs(self.end_segment)]

        # The start segment and end segment should agree in depth. If they don't, that's very bad,
        # so depth_disagreement is applied to quality twice (squared effect).
        depth_agreement = get_num_agreement(start_seg.depth, end_seg.depth)
        self.quality *= (depth_agreement * depth_agreement) # has squared effect on quality
        self.depth = get_mean_depth(start_seg, end_seg, graph)

        # If the segments in the path exclusively lead to the start and end segments (i.e they
        # cannot lead to any another segment), then we can also scale the quality based on the
        # depth consistency of the path. E.g. if a bridge path contains a segment 3 times and that
        # segment's depth also suggests about 3 times, that's good. If they don't agree, that's
        # bad.
        if path_is_self_contained(self.graph_path, self.start_segment, self.end_segment, graph):
            graph_path_pos_nums = list(set([abs(x) for x in self.graph_path]))
            for path_segment in graph_path_pos_nums:
                actual_depth = graph.segments[path_segment].depth
                expected_depth = graph_path_pos_nums.count(path_segment) * self.depth
                agreement = get_num_agreement(actual_depth, expected_depth)
                self.quality *= agreement

    def __repr__(self):
        return 'SPAdes bridge: ' + get_bridge_str(self.start_segment, self.graph_path,
                                                  self.end_segment) + \
               ' (quality = ' + float_to_str(self.quality, 2) + ')'


class LongReadBridge(object):
    '''
    This class describes a bridge created from long read alignments.
    '''
    def __init__(self, graph, start, end):

        # The numbers of the two single-copy segments which are being bridged.
        self.start_segment = start
        self.end_segment = end

        # The individual reads contributing to the bridge. The sequences/qualities are not for the
        # entire read, just the part in the bridge. The lengths are stored separately because
        # negative bridge lengths are possible (when the start and end segment alignments overlap).
        # In these cases we don't have any read sequences but will still need to know the lengths
        # (which will be negative).
        self.full_span_reads = []

        # The bridge can also contain incomplete read sequences which don't bridge the entire span
        # between the start and end segments. These are still useful as they can contribute to the
        # consensus sequence.
        self.start_only_reads = []
        self.end_only_reads = []

        # The consensus of all read sequences. If there is only one read, then this is the same as
        # that read's sequence. If there are multiple reads, this is hopefully a more accurate
        # sequence. If this is a case where the start and end segments overlap, then there will not
        # be a consensus sequence and consensus_length will be the mean of read_lengths (a negative
        # value).
        self.consensus_sequence = ''

        # The path through the unbridged graph, if one was found.
        self.graph_path = []

        # The bridge sequence, gotten from the graph path if a good path was found. Otherwise it's
        # from the consensus sequence.
        self.bridge_sequence = ''

        # The bridge depth, a weighted mean of the start and end depths.
        self.depth = get_mean_depth(graph.segments[abs(self.start_segment)],
                                    graph.segments[abs(self.end_segment)], graph)

        # A score used to determine the order of bridge application.
        self.quality = 1.0

        # Whether or not the long read bridge got its sequence from a graph path (True) or from the
        # read consensus (False).
        self.path_support = None

        self.graph = graph

    def __repr__(self):
        return 'long read bridge: ' + get_bridge_str(self.start_segment, self.graph_path,
                                                     self.end_segment) + \
               ' (quality = ' + float_to_str(self.quality, 2) + ')'

    def finalise(self, scoring_scheme, min_alignment_length, mean_spans_per_bridge, verbosity):
        '''
        Determines the consensus sequence for the bridge, attempts to find it in the graph and
        assigns a quality score to the bridge. This is the performance-intensive step of long read
        bridging.
        '''
        start_seg = self.graph.segments[abs(self.start_segment)]
        end_seg = self.graph.segments[abs(self.end_segment)]

        output = '\n'

        # output += 'FINALISING BRIDGE\n'
        # output += '-----------------\n'
        # output += 'start: ' + str(self.start_segment) + '\n'
        # output += 'end:   ' + str(self.end_segment) + '\n'
        # output += 'start overlaps:\n'
        # for start_only_read in self.start_only_reads:
        #     output += '  ' + str(start_only_read) + '\n'
        # output += 'end overlaps:\n'
        # for end_only_read in self.end_only_reads:
        #     output += '  ' + str(end_only_read) + '\n'
        # output += 'full spans:\n'
        # for full_span_read in self.full_span_reads:
        #     output += '  ' + str(full_span_read) + '\n'

        output += str(self.start_segment) + ' to ' + str(self.end_segment) + ':\n'
        output += '  bridging reads:          ' + int_to_str(len(self.full_span_reads)) + '\n'

        # Parition the full span reads into two groups: those with negative numbers (implying that
        # the two segments overlap) and those with actual sequences.
        full_spans_without_seq = []
        full_spans_with_seq = []
        for full_span in self.full_span_reads:
            if isinstance(full_span[0], int):
                full_spans_without_seq.append(full_span)
            else:
                full_spans_with_seq.append(full_span)

        # There shouldn't usually be both full spans with sequence and without. If there are some
        # of each, we'll throw out the minority group.
        if full_spans_with_seq and full_spans_without_seq:
            if len(full_spans_without_seq) > len(full_spans_with_seq):
                full_spans_with_seq = []
            else:
                full_spans_without_seq = []

        # For full spans with sequence, we perform a MSA and get a consensus sequence.
        if full_spans_with_seq:

            # Full-span sequences are faster in the MSA and are very useful for the consensus
            # sequence, so we want to use a lot! But we still set an upper limit for cases of
            # very high read depth.
            max_full_span = 50 # TO DO: make this a parameter?
            if len(full_spans_with_seq) <= max_full_span:
                full_span_seqs = [x[0] for x in full_spans_with_seq]
                full_span_quals = [x[1] for x in full_spans_with_seq]
            else:
                sorted_full_span = sorted(full_spans_with_seq,
                                          key=lambda x: x[2].raw_score + x[3].raw_score,
                                          reverse=True)
                full_span_seqs = [x[0] for x in sorted_full_span[:max_full_span]]
                full_span_quals = [x[1] for x in sorted_full_span[:max_full_span]]


            # Start-only and end-only sequences make the MSA a lot slower, contribute less to the
            # consensus and can screw up the MSA if they are too abundant. So if there are too many
            # of them, we only include the ones with the best alignments.
            max_partial = min(5, len(full_span_seqs)) # TO DO: make this a parameter?
            if len(self.start_only_reads) <= max_partial:
                start_only_seqs = [x[0] for x in self.start_only_reads]
                start_only_quals = [x[1] for x in self.start_only_reads]
            else:
                sorted_start_only = sorted(self.start_only_reads,
                                           key=lambda x: x[2].raw_score, reverse=True)
                start_only_seqs = [x[0] for x in sorted_start_only[:max_partial]]
                start_only_quals = [x[1] for x in sorted_start_only[:max_partial]]

            if len(self.end_only_reads) <= max_partial:
                end_only_seqs = [x[0] for x in self.end_only_reads]
                end_only_quals = [x[1] for x in self.end_only_reads]
            else:
                sorted_end_only = sorted(self.end_only_reads,
                                         key=lambda x: x[2].raw_score, reverse=True)
                end_only_seqs = [x[0] for x in sorted_end_only[:max_partial]]
                end_only_quals = [x[1] for x in sorted_end_only[:max_partial]]

            consensus_start_time = time.time()
            self.consensus_sequence, _, _, _ = \
                                    multiple_sequence_alignment(full_span_seqs, full_span_quals,
                                                                start_only_seqs, start_only_quals,
                                                                end_only_seqs, end_only_quals,
                                                                scoring_scheme)
            # output += 'consensus: ' + str(self.consensus_sequence) + '\n'
            # output += 'full span consensus scores: ' + str(full_span_scores) + '\n'
            # if start_only_scores:
            #     output += 'start-only consensus scores: ' + str(start_only_scores) + '\n'
            # if end_only_scores:
            #     output += 'end-only consensus scores: ' + str(end_only_scores) + '\n'

            consensus_time = time.time() - consensus_start_time
            output += '  consensus sequence:      ' + \
                      int_to_str(len(self.consensus_sequence)) + ' bp '
            output += '(' + float_to_str(consensus_time, 2) + ' sec)\n'

            target_path_length = len(self.consensus_sequence) + (2 * self.graph.overlap)


        # For full spans without sequence, we simply need a mean distance.
        elif full_spans_without_seq:
            self.consensus_sequence = ''
            mean_overlap = int(round(sum(x[0] for x in full_spans_without_seq) / \
                                     len(full_spans_without_seq)))
            output += '  mean overlap:            ' + int_to_str(abs(mean_overlap)) + '\n'
            target_path_length = mean_overlap + (2 * self.graph.overlap)

        output += '  target path length:      ' + int_to_str(target_path_length) + ' bp\n'

        # Limit the path search to lengths near the target.
        # TO DO: adjust these or make them parameters?
        min_path_length = int(round(target_path_length * 0.5))
        max_path_length = int(round(target_path_length * 1.5))

        # output += 'min graph path length: ' + str(min_path_length) + '\n'
        # output += 'max graph path length: ' + str(max_path_length) + '\n'

        max_path_count = 100 # TO DO: make this a parameter?
        path_start_time = time.time()
        potential_paths = self.graph.all_paths(self.start_segment, self.end_segment,
                                               min_path_length, target_path_length,
                                               max_path_length, max_path_count)
        path_time = time.time() - path_start_time

        output += '  path count:              ' + int_to_str(len(potential_paths)) + ' '
        output += '(' + float_to_str(path_time, 2) + ' sec)\n'

        # output += 'potential paths:\n'
        # if not potential_paths:
        #     output += '  NONE\n'

        self.graph_path = []

        # If we have a consensus to align to, then we use that do choose the best path.
        if self.consensus_sequence:
            best_raw_score = None
            best_scaled_score = None
            alignment_start_time = time.time()
            for path in potential_paths:

                path_seq = self.graph.get_path_sequence(path)

                # output += '  ' + str(path)
                # output += ' (' + str(len(path_seq)) + '), '

                alignment_result = fully_global_alignment(self.consensus_sequence, path_seq,
                                                          scoring_scheme, True, 1000)
                if not alignment_result:
                    continue

                seqan_parts = alignment_result.split(',', 9)
                raw_score = int(seqan_parts[6])
                scaled_score = float(seqan_parts[7])
                # output += str(raw_score) + ', ' + str(scaled_score) + '\n'

                if best_raw_score is None or raw_score > best_raw_score:
                    best_raw_score = raw_score
                    best_scaled_score = scaled_score
                    self.graph_path = path

                # In case of a tie, go with the simpler path.
                elif raw_score == best_raw_score and len(path) < len(self.graph_path):
                    best_scaled_score = scaled_score
                    self.graph_path = path

            alignment_time = time.time() - alignment_start_time

        # If there isn't a consensus (i.e. the start and end overlap), then we choose the best path
        # based on its length.
        else:
            smallest_length_diff = None
            for path in potential_paths:
                path_len = self.graph.get_path_length(path)

                # output += '  ' + str(path)
                # output += ' (' + str(path_len) + ')\n'

                length_diff = abs(path_len - target_path_length)
                if smallest_length_diff is None or length_diff < smallest_length_diff:
                    smallest_length_diff = length_diff
                    self.graph_path = path
                elif smallest_length_diff == length_diff and len(path) < len(self.graph_path):
                    self.graph_path = path
            if self.graph_path:
                best_scaled_score = get_num_agreement(self.graph.get_path_length(self.graph_path),
                                                      target_path_length) * 100.0
            alignment_time = None

        # If a path was found, use its sequence for the bridge.
        if self.graph_path:
            output += '  best path:               ' + \
                      ', '.join(int_to_str(x) for x in self.graph_path) + ' '
            output += '(' + int_to_str(self.graph.get_path_length(self.graph_path)) + ' bp'
            if alignment_time:
                output += ', ' + float_to_str(alignment_time, 2) + ' sec)\n'
            else:
                output += ')\n'

            self.bridge_sequence = self.graph.get_path_sequence(self.graph_path)
            self.path_support = True

            # Now we adjust the quality. It starts on the scaled score of the path alignment (which
            # maxes out at 100.0).
            self.quality = best_scaled_score

        # If a path wasn't found, the consensus sequence is the bridge (with the overlaps added).
        else:
            output += '  best path:               none found\n'
            start_overlap = \
                self.graph.get_seq_from_signed_seg_num(self.start_segment)[-self.graph.overlap:]
            end_overlap = \
                self.graph.get_seq_from_signed_seg_num(self.end_segment)[:self.graph.overlap]
            self.bridge_sequence = start_overlap + self.consensus_sequence + end_overlap
            self.path_support = False

            # Non-graph-path-supported bridges are much lower quality than graph-path-supported
            # bridges.
            self.quality = 20.0

        if verbosity > 2:
            output += '  starting quality:        ' + float_to_str(self.quality, 2) + '\n'

        # The start segment and end segment should agree in depth. If they don't, that's very bad,
        # as it implies that they aren't actually single-copy or on the same piece of DNA.
        depth_agreement_factor = get_num_agreement(start_seg.depth, end_seg.depth)

        # The number of reads which contribute to a bridge is a big deal, so the read count factor
        # scales linearly. This is adjustes to the bridge's depth, so higher depth bridges need a
        # larger read count to get the same effect. However, the depth used here can't go below 1,
        # as that could give an unfair advantage to bridges between very low depth segments (which
        # we probably don't really care about anyways because they are low depth).
        read_count_factor = len(self.full_span_reads) / mean_spans_per_bridge / max(self.depth, 1)

        # The length of alignments to the start/end segments is positively correlated with quality
        # to reward bridges with long alignments.
        total_alignment_length = sum(x[2].get_aligned_ref_length() + x[3].get_aligned_ref_length() \
                                     for x in self.full_span_reads)
        mean_alignment_length = total_alignment_length / (2.0 * len(self.full_span_reads))
        align_length_factor = score_function(mean_alignment_length, min_alignment_length * 4)

        # The mean alignment score to the start/end segments is positively correlated with quality,
        # so bridges with high quality alignments are rewarded.
        scaled_score_total = sum(x[2].scaled_score + x[3].scaled_score \
                                 for x in self.full_span_reads)
        mean_scaled_score = scaled_score_total / (2.0 * len(self.full_span_reads))
        align_score_factor = mean_scaled_score / 100.0

        # Bridges between long start/end segments are rewarded, as they are more likely to actually
        # be single-copy.
        start_length_factor = score_function(start_seg.get_length(), min_alignment_length)
        end_length_factor = score_function(end_seg.get_length(), min_alignment_length)

        self.quality *= depth_agreement_factor * depth_agreement_factor
        self.quality *= read_count_factor
        self.quality *= align_length_factor
        self.quality *= align_score_factor
        self.quality *= start_length_factor
        self.quality *= end_length_factor

        if verbosity > 2:
            output += '  depth agreement factor:  ' + float_to_str(depth_agreement_factor, 2) + '\n'
            output += '  read count factor:       ' + float_to_str(read_count_factor, 2) + '\n'
            output += '  alignment length factor: ' + float_to_str(align_length_factor, 2) + '\n'
            output += '  alignment score factor:  ' + float_to_str(align_score_factor, 2) + '\n'
            output += '  start length factor:     ' + float_to_str(start_length_factor, 2) + '\n'
            output += '  end length factor:       ' + float_to_str(end_length_factor, 2) + '\n'
            output += '  final quality:           ' + float_to_str(self.quality, 2) + '\n'

        return output

    def contains_full_span_sequence(self):
        '''
        Some LongReadBridge objects bridge two close segments, and therefore do not actually have
        any bridging read sequence (just a negative number implying overlap). This function returns
        True if any full span read sequences exist and False if they are just negative numbers.
        '''
        for full_span_read in self.full_span_reads:
            seq_or_num = full_span_read[0]
            if not isinstance(seq_or_num, int):
                return True
        return False



class LoopUnrollingBridge(object):
    '''
    This class describes a bridge created from unrolling an assembly graph loop.

    Quality is affected by:
      * How well the start and end segments' depths agree.
      * How close the determined loop count is to a whole number.
      * The final loop count (higher counts get lower quality).
    '''
    def __init__(self, graph, start, end, middle, repeat):
        '''
        This constructor assumes the the start, end, middle and repeat segments form a simple loop
        in the graph supported by either a SPAdes contig or a long read alignment. It will use
        segment depths to determine the loop count and score the bridge's quality.
        '''
        # The numbers of the two single-copy segments which are being bridged.
        self.start_segment = start
        self.end_segment = end

        # The path through the unbridged graph.
        self.graph_path = []

        # The bridge sequence, gotten from the graph path.
        self.bridge_sequence = ''

        # The bridge depth, a weighted mean of the start and end depths.
        self.depth = 0.0

        # A score used to determine the order of bridge application. This value starts at the
        # maximum for a loop unrolling bridge and can only decrease as the constructor continues.
        self.quality = 30.0

        # Get the actual segments from the numbers. Since we are assuming they do form a simple
        # loop, we don't care about directionality.
        start_seg = graph.segments[abs(start)]
        end_seg = graph.segments[abs(end)]
        middle_seg = graph.segments[abs(middle)]
        repeat_seg = graph.segments[abs(repeat)]

        # The start segment and end segment should agree in depth. If they don't, that's very bad,
        # so depth_disagreement is applied to quality twice (squared effect).
        depth_agreement = get_num_agreement(start_seg.depth, end_seg.depth)
        self.quality *= (depth_agreement * depth_agreement) # has squared effect on quality

        # We'll use a mean loop count that's weighted by the middle and repeat segment lengths.
        self.depth = get_mean_depth(start_seg, end_seg, graph)
        loop_count_by_middle = middle_seg.depth / self.depth
        loop_count_by_repeat = (repeat_seg.depth - self.depth) / self.depth
        mean_loop_count = weighted_average(loop_count_by_middle, loop_count_by_repeat,
                                           middle_seg.get_length_no_overlap(graph.overlap),
                                           repeat_seg.get_length_no_overlap(graph.overlap))

        # If the average loop count is near a whole number, that's better. If it's near 0.5, that's
        # very bad!
        if mean_loop_count < 1.0:
            loop_count = 1
            closeness_to_whole_num = mean_loop_count
        else:
            loop_count = int(round(mean_loop_count))
            fractional_part = mean_loop_count % 1
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
        return 'loop bridge: ' + get_bridge_str(self.start_segment, self.graph_path,
                                                self.end_segment) + \
               ' (quality = ' + float_to_str(self.quality, 2) + ')'

def create_spades_contig_bridges(graph, single_copy_segments, verbosity):
    '''
    Builds graph bridges using the SPAdes contig paths.
    '''
    bridge_path_set = set()
    single_copy_numbers = [x.number for x in single_copy_segments]
    for segment in single_copy_segments:
        for path in graph.paths.values():
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
    for grouped_paths in bridge_paths_by_start.values():
        if len(grouped_paths) > 1:
            conflicting_paths += grouped_paths
    for grouped_paths in bridge_paths_by_end.values():
        if len(grouped_paths) > 1:
            conflicting_paths += grouped_paths
    conflicting_paths_no_dups = []
    for path in conflicting_paths:
        if path not in conflicting_paths_no_dups:
            conflicting_paths_no_dups.append(path)
    conflicting_paths = conflicting_paths_no_dups
    # if verbosity > 1:
    #     print('Bridge paths in conflict with single-copy segments: ', end='')
    #     if conflicting_paths:
    #         print(', '.join([str(x) for x in conflicting_paths]))
    #     else:
    #         print('none')
    #     print()

    final_bridge_paths = [x for x in bridge_path_list if x not in conflicting_paths]
    # if verbosity > 1:
    #     print('Final SPAdes contig bridge paths: ', end='')
    #     if final_bridge_paths:
    #         print(', '.join([str(x) for x in final_bridge_paths]))
    #     else:
    #         print('none')
    #     print()

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
        for path in graph.paths.values():
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

def create_long_read_bridges(graph, read_dict, read_names, single_copy_segments, verbosity,
                             existing_bridges, min_scaled_score, threads, scoring_scheme,
                             min_alignment_length):
    '''
    Makes bridges between single-copy segments using the alignments in the long reads.
    '''
    single_copy_seg_num_set = set()
    for seg in single_copy_segments:
        single_copy_seg_num_set.add(seg.number)

    # This dictionary will collect the read sequences which span between two single-copy segments.
    # These are the most useful sequences and will be used to either create a new bridge or enhance
    # an existing bridge.
    # Key = tuple of signed segment numbers (the segments being bridged)
    # Value = list of tuples containing the bridging sequence and the single-copy segment
    #         alignments.
    spanning_read_seqs = {}

    # This dictionary will collect all of the read sequences which don't manage to span between
    # two single-copy segments, but they do overlap one single-copy segment and can therefore be
    # useful for generating consensus sequences in a bridge.
    # Key = signed segment number, Value = list of tuples containing the sequence and alignment
    overlapping_read_seqs = {}

    allowed_overlap = round(int(1.1 * graph.overlap))
    for read_name in read_names:
        read = read_dict[read_name]
        alignments = get_single_copy_alignments(read, single_copy_seg_num_set, allowed_overlap,
                                                min_scaled_score)
        if not alignments:
            continue

        # print('\n')
        # print('READ:', read)
        # print('  SEQUENCE:', read.sequence)
        # print('  SINGLE-COPY SEGMENT ALIGNMENTS:')
        # for alignment in alignments:
        #     print('    ', alignment)
        # if len(alignments) > 1:
        #     print('  BRIDGING SEQUENCES:')

        # If the code got here, then we have some alignments to single-copy segments. We grab
        # neighbouring pairs of alignments, starting with the highest scoring ones and work our
        # way down. This means that we should have a pair for each neighbouring alignment, but
        # potentially also more distant pairs if the alignments are strong.
        already_added = set()
        sorted_alignments = sorted(alignments, key=lambda x: x.raw_score, reverse=True)
        available_alignments = []
        for alignment in sorted_alignments:

            # If the alignment being added is to a reference that has already been added but in the
            # opposite direction, then we don't include it. E.g. we don't add an alignment for 10
            # if we already have an alignment for -10. This is because there's no legitimate way
            # for a single-copy segment to appear in the same read in two different directions. The
            # same direction is okay, as that can happen with a circular piece of DNA, but opposite
            # directions implies multi-copy.
            opposite_num = -alignment.get_signed_ref_num()
            if opposite_num in set(x.get_signed_ref_num() for x in available_alignments):
                continue

            available_alignments.append(alignment)
            available_alignments = sorted(available_alignments,
                                          key=lambda x: x.read_start_positive_strand())
            for i in range(len(available_alignments) - 1):
                alignment_1 = available_alignments[i]
                alignment_2 = available_alignments[i+1]

                # Standardise the order so we don't end up with both directions (e.g. 5 to -6 and
                # 6 to -5) in spanning_read_seqs.
                seg_nums, flipped = flip_segment_order(alignment_1.get_signed_ref_num(),
                                                       alignment_2.get_signed_ref_num())
                if seg_nums not in already_added:
                    bridge_start = alignment_1.read_end_positive_strand()
                    bridge_end = alignment_2.read_start_positive_strand()

                    if bridge_end > bridge_start:
                        bridge_seq = read.sequence[bridge_start:bridge_end]
                        bridge_qual = read.qualities[bridge_start:bridge_end]
                        if flipped:
                            bridge_seq = reverse_complement(bridge_seq)
                            bridge_qual = bridge_qual[::-1]
                    else:
                        bridge_seq = bridge_end - bridge_start # 0 or a negative number
                        bridge_qual = ''

                    if seg_nums not in spanning_read_seqs:
                        spanning_read_seqs[seg_nums] = []

                    spanning_read_seqs[seg_nums].append((bridge_seq, bridge_qual, alignment_1,
                                                         alignment_2))
                    already_added.add(seg_nums)

                    # print('    ', seg_nums[0], seg_nums[1], bridge_seq)
                    # print('    ', seg_nums[0], seg_nums[1], bridge_qual)

        # At this point all of the alignments have been added and we are interested in the first
        # and last alignments (which may be the same if there's only one). If the read extends
        # past these alignments, then the overlapping part might be useful for consensus sequences.
        first_alignment = available_alignments[0]
        start_overlap, start_qual = first_alignment.get_start_overlapping_read_seq()
        if start_overlap:
            seg_num = -first_alignment.get_signed_ref_num()
            seq = reverse_complement(start_overlap)
            qual = start_qual[::-1]
            if seg_num not in overlapping_read_seqs:
                overlapping_read_seqs[seg_num] = []
            overlapping_read_seqs[seg_num].append((seq, qual, first_alignment))
            # print('  START OVERLAPPING SEQUENCE:')
            # print('    ', seg_num, seq)
            # print('    ', seg_num, qual)
        last_alignment = available_alignments[-1]
        end_overlap, end_qual = last_alignment.get_end_overlapping_read_seq()
        if end_overlap:
            seg_num = last_alignment.get_signed_ref_num()
            if seg_num not in overlapping_read_seqs:
                overlapping_read_seqs[seg_num] = []
            overlapping_read_seqs[seg_num].append((end_overlap, end_qual, last_alignment))
        #     print('  END OVERLAPPING SEQUENCE:')
        #     print('    ', seg_num, end_overlap)
        #     print('    ', seg_num, end_qual)
        # print('\n')

    # If a bridge already exists for a spanning sequence, we add the sequence to the bridge. If
    # not, we create a new bridge and add it.
    new_bridges = []
    for seg_nums, span in spanning_read_seqs.items():
        start, end = seg_nums
        for existing_bridge in existing_bridges:
            if isinstance(existing_bridge, LongReadBridge) and \
               existing_bridge.start_segment == start and existing_bridge.end_segment == end:
                matching_bridge = existing_bridge
                break
        else:
            new_bridge = LongReadBridge(graph, start, end)
            new_bridges.append(new_bridge)
            matching_bridge = new_bridge
        matching_bridge.full_span_reads += span
    all_bridges = existing_bridges + new_bridges
    all_bridges = sorted(all_bridges, key=lambda x: (x.start_segment, x.end_segment))

    # Add overlapping sequences to appropriate bridges, but only if they have some full span
    # sequence (if their full span reads only have overlap-indicating negative numbers, then we
    # won't be doing a consensus sequence and don't need overlapping sequences).
    for seg_num, overlaps in overlapping_read_seqs.items():
        for bridge in all_bridges:
            if not isinstance(bridge, LongReadBridge) or not bridge.contains_full_span_sequence():
                continue
            start_overlap = (bridge.start_segment == seg_num)
            end_overlap = (bridge.end_segment == -seg_num)
            if start_overlap or end_overlap:
                for overlap in overlaps:
                    if start_overlap:
                        bridge.start_only_reads.append(overlap)
                    elif end_overlap:
                        overlap_seq, overlap_qual, alignment = overlap
                        bridge.start_only_reads.append((reverse_complement(overlap_seq),
                                                        overlap_qual[::-1], alignment))

    # Figure out the average number of spanning reads per bridge, normalised to the bridge's depth.
    total_normalised_reads = 0.0
    long_read_bridge_count = 0
    for bridge in all_bridges:
        if not isinstance(bridge, LongReadBridge) or not bridge.contains_full_span_sequence():
            continue
        long_read_bridge_count += 1
        total_normalised_reads += len(bridge.full_span_reads) / bridge.depth
    mean_spans_per_bridge = total_normalised_reads / long_read_bridge_count

    # Now we need to finalise the reads. This is the intensive step, as it involves creating a
    # consensus sequence, finding graph paths and doing alignments between the consensus and the
    # graph paths. We therefore use available threads to make this faster.
    long_read_bridges = [x for x in all_bridges if isinstance(x, LongReadBridge)]
    num_long_read_bridges = len(long_read_bridges)
    completed_count = 0
    if verbosity == 1:
        print_progress_line(0, num_long_read_bridges, prefix='Bridge: ')
    completed_count = 0
    if threads == 1:
        for bridge in long_read_bridges:
            output = bridge.finalise(scoring_scheme, min_alignment_length, mean_spans_per_bridge,
                                     verbosity)
            completed_count += 1
            if verbosity == 1:
                print_progress_line(completed_count, num_long_read_bridges, prefix='Bridge: ')
            if verbosity > 1:
                print(output, end='')
    else:
        pool = ThreadPool(threads)
        arg_list = []
        for bridge in long_read_bridges:
            arg_list.append((bridge, scoring_scheme, min_alignment_length, mean_spans_per_bridge,
                             verbosity))

        # If the verbosity is 1, then the order doesn't matter, so use imap_unordered to deliver
        # the results evenly. If the verbosity is higher, deliver the results in order with imap.
        if verbosity > 1:
            imap_function = pool.imap
        else:
            imap_function = pool.imap_unordered
        for output in imap_function(finalise_bridge, arg_list):
            completed_count += 1
            if verbosity == 1:
                print_progress_line(completed_count, num_long_read_bridges, prefix='Bridge: ')
            if verbosity > 1:
                print(output, end='')

    if verbosity == 1:
        print()

    return all_bridges

def get_num_agreement(num_1, num_2):
    '''
    Returns a value between 0.0 and 1.0 describing how well the numbers agree.
    1.0 is perfect agreement and 0.0 is the worst.
    '''
    if num_1 == 0.0 and num_2 == 0.0:
        return 1.0
    if num_1 < 0.0 and num_2 < 0.0:
        num_1 = -num_1
        num_2 = -num_2
    if num_1 * num_2 < 0.0:
        return 0.0
    return min(num_1, num_2) / max(num_1, num_2)

def get_mean_depth(seg_1, seg_2, graph):
    '''
    Returns the mean depth of the two segments, weighted by their length.
    '''
    return weighted_average(seg_1.depth, seg_2.depth,
                            seg_1.get_length_no_overlap(graph.overlap),
                            seg_2.get_length_no_overlap(graph.overlap))

def path_is_self_contained(path, start, end, graph):
    '''
    Returns True if the path segments are only connected to each other and the start/end segments.
    If they are connected to anything else, it returns False.
    '''
    all_numbers_in_path = set()
    all_numbers_in_path.add(abs(start))
    all_numbers_in_path.add(abs(end))
    for segment in path:
        all_numbers_in_path.add(abs(segment))
    for segment in path:
        connected_segments = graph.get_connected_segments(segment)
        for connected_segment in connected_segments:
            if connected_segment not in all_numbers_in_path:
                return False
    return True

def flip_segment_order(seg_num_1, seg_num_2):
    '''
    Given two segment numbers, this function possibly flips them around. It returns the new numbers
    (either unchanged or flipped) and whether or not a flip took place. The decision is somewhat
    arbitrary, but it needs to be consistent so when we collect bridging read sequences they are
    always in the same direction.
    '''
    if seg_num_1 > 0 and seg_num_2 > 0:
        flip = False
    elif seg_num_1 < 0 and seg_num_2 < 0:
        flip = True
    elif seg_num_1 < 0: # only seg_num_1 is negative
        flip = abs(seg_num_1) > abs(seg_num_2)
    else: # only seg_num_2 is negative
        flip = abs(seg_num_2) > abs(seg_num_1)
    if flip:
        return (-seg_num_2, -seg_num_1), True
    else:
        return (seg_num_1, seg_num_2), False

def get_single_copy_alignments(read, single_copy_num_set, allowed_overlap, min_scaled_score):
    '''
    Returns a list of single-copy segment alignments for the read.
    '''
    sc_alignments = []
    for alignment in read.alignments:
        if alignment.ref.number in single_copy_num_set and \
           alignment.scaled_score >= min_scaled_score:
            sc_alignments.append(alignment)
    return sc_alignments

def finalise_bridge(all_args):
    '''
    Just a one-argument version of bridge.finalise, for pool.imap.
    '''
    bridge, scoring_scheme, min_alignment_length, mean_spans_per_bridge, verbosity = all_args
    return bridge.finalise(scoring_scheme, min_alignment_length, mean_spans_per_bridge, verbosity)

def score_function(val, half_score_val):
    '''
    For inputs of 0.0 and greater, this function returns a value between 0.0 and 1.0, approaching
    1.0 with large values. The half_score_val argument is the point at which the function returns
    0.5. If it's large the function approaches 1.0 more slowly, if it's small the function
    approaches 1.0 more quickly.
    '''
    return 1.0 - (half_score_val / (half_score_val + val))

def get_applicable_bridge_pieces(bridge, single_copy_nums, right_bridged, left_bridged,
                                 seg_nums_used_in_bridges):
    '''
    Using the given sets, this function returns a list of bridge pieces which can be applied to the
    graph.
    '''
    # Break the bridge into pieces based on single-copy segments in the bridge path. Each piece is
    # grouped with its index which make it easier to combine the adjacent pieces later.
    bridge_pieces = []
    current_piece = [(bridge.start_segment, 0)]
    for i, seg in enumerate(bridge.graph_path):
        if abs(seg) in single_copy_nums:
            bridge_pieces.append(current_piece + [(seg, i+1)])
            current_piece = [(seg, i+1)]
        else:
            current_piece.append((seg, i+1))
    bridge_pieces.append(current_piece + [(bridge.end_segment, len(bridge.graph_path) + 1)])

    # A piece can only be applied if its start/end haven't already been bridged or used in a
    # previous bridge.
    pieces_to_apply = []
    for piece in bridge_pieces:
        if start_end_available_to_bridge(piece[0][0], piece[-1][0], right_bridged,
                                         left_bridged, seg_nums_used_in_bridges):
            pieces_to_apply.append(piece)
    if not pieces_to_apply:
        return []

    # Merge adjacent pieces together again so they can be applied at once.
    merged_pieces_to_apply = [pieces_to_apply[0]]
    for piece in pieces_to_apply[1:]:
        if piece[0] == merged_pieces_to_apply[-1][-1]:
            merged_pieces_to_apply[-1] = merged_pieces_to_apply[-1] + piece[1:]
        else:
            merged_pieces_to_apply.append(piece)

    # We no longer need those indicies, so simplify the list to just segment numbers.
    return [[y[0] for y in x] for x in merged_pieces_to_apply]

def get_bridge_str(start, middle, end):
    '''
    Returns a bridge sequence in human-readable form.
    '''
    bridge_str = str(start) + ' -> '
    if middle:
        bridge_str += ', '.join([str(x) for x in middle]) + ' -> '
    bridge_str += str(end)
    return bridge_str

def start_end_available_to_bridge(start, end, right_bridged, left_bridged,
                                  seg_nums_used_in_bridges):
    '''
    Checks whether the start and end segments can be bridged together (i.e. that they are both
    unbridged on the relevant sides).
    '''
    if start > 0 and start in right_bridged:
        return False
    if start < 0 and -start in left_bridged:
        return False
    if end > 0 and end in left_bridged:
        return False
    if end < 0 and -end in right_bridged:
        return False
    if abs(start) in seg_nums_used_in_bridges:
        return False
    if abs(end) in seg_nums_used_in_bridges:
        return False
    return True


