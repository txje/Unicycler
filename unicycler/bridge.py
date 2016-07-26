"""
This module describes bridges - links between two single-copy segments in an assembly graph.
Bridges can come from multiple sources, so there are a few different classes which all share
the important members/methods (for duck-typing purposes).

Author: Ryan Wick
email: rrwick@gmail.com
"""

from multiprocessing.dummy import Pool as ThreadPool
import time
from collections import defaultdict
from .misc import int_to_str, float_to_str, reverse_complement, print_progress_line, \
    weighted_average, get_num_agreement, flip_number_order, score_function
from .cpp_function_wrappers import multiple_sequence_alignment


class SpadesContigBridge(object):
    """
    This class describes a bridge created from the contigs.paths file made by SPAdes.

    Quality is affected by:
      * How well the start and end segments' depths agree.
      * The depth consistency within the path (only applies to bridges where the path segments
        exclusively lead to the start/end segments).
    """

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

        # If there are segments in between the start and end (there usually will be), then they
        # give the bridge sequence. If not (i.e. if the start and end directly connect),
        # then the bridge sequence is just the overlapping sequence between them.
        if self.graph_path:
            self.bridge_sequence = graph.get_path_sequence(self.graph_path)
        else:
            overlap_seq = graph.seq_from_signed_seg_num(self.start_segment)[-graph.overlap:]
            self.bridge_sequence = overlap_seq

        # The start segment and end segment should agree in depth. If they don't, that's very bad,
        # so depth_disagreement is applied to quality twice (squared effect).
        start_seg = graph.segments[abs(self.start_segment)]
        end_seg = graph.segments[abs(self.end_segment)]
        depth_agreement = get_num_agreement(start_seg.depth, end_seg.depth)
        self.quality *= (depth_agreement * depth_agreement)  # has squared effect on quality
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

        # The quality of a SPAdes contig bridge should decline sharply if the bridging sequence
        # is too large compared to the short read insert size. E.g. if the insert size is 500 bp
        # and the SPAdes contig bridge is 2 kb long, we should not believe it!
        if self.graph_path:
            bridge_length = len(self.bridge_sequence)
            if bridge_length <= graph.insert_size_mean:
                bridge_length_factor = 1.0
            else:
                bridge_length_factor = graph.insert_size_deviation / (bridge_length -
                                                                      graph.insert_size_mean +
                                                                      graph.insert_size_deviation)
            self.quality *= bridge_length_factor

    def __repr__(self):
        return 'SPAdes bridge: ' + get_bridge_str(self.start_segment, self.graph_path,
                                                  self.end_segment) + \
               ' (quality = ' + float_to_str(self.quality, 2) + ')'

    @staticmethod
    def get_type_score():
        """
        Returns a score indicating the relative importance of the bridge types:
        LongReadBridge = 2, SpadesContigBridge = 1, LoopUnrollingBridge = 0
        """
        return 1


class LongReadBridge(object):
    """
    This class describes a bridge created from long read alignments.
    """

    def __init__(self, graph, start, end):

        # The numbers of the two single-copy segments which are being bridged.
        self.start_segment = start
        self.end_segment = end

        # The individual reads contributing to the bridge. The sequences/qualities are not for the
        # entire read, just the part in the bridge. In the case of overlapping alignments, this has
        # the overlap size instead of the sequence.
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
        self.all_best_paths = []

        # The bridge sequence, gotten from the graph path if a good path was found. Otherwise it's
        # from the consensus sequence. If there are multiple best paths, then all_bridge_sequences
        # will contain the sequence for each.
        self.bridge_sequence = ''
        self.all_bridge_sequences = []

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

    def get_total_full_span_seq_length(self):
        """
        Returns the sum of all available read sequences. This is used to prioritise the bridge
        finalisation.
        """
        total_full_span_seq_length = 0
        for full_span in self.full_span_reads:
            if not isinstance(full_span[0], int):
                total_full_span_seq_length += len(full_span[0])
        return total_full_span_seq_length

    def finalise(self, scoring_scheme, min_alignment_length, read_lengths, estimated_genome_size,
                 verbosity):
        """
        Determines the consensus sequence for the bridge, attempts to find it in the graph and
        assigns a quality score to the bridge. This is the performance-intensive step of long read
        bridging.
        """
        start_seg = self.graph.segments[abs(self.start_segment)]
        end_seg = self.graph.segments[abs(self.end_segment)]

        output = '\n'
        output += str(self.start_segment) + ' to ' + str(self.end_segment) + ':\n'
        output += '  bridging reads:          ' + int_to_str(len(self.full_span_reads)) + '\n'

        # Partition the full span reads into two groups: those with negative numbers (implying that
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
            # very high read depth, as getting the consensus sequence will be very slow if there
            # are too many reads.
            max_full_span = 25  # TO DO: make this a parameter?
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
            # Also, if we have a LOT of full span reads, there's no need to include any partial
            # sequences at all.
            max_partial = min(5, len(full_span_seqs))  # TO DO: make this a parameter?
            if len(full_span_seqs) == max_full_span:
                max_partial = 0
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
            self.consensus_sequence = multiple_sequence_alignment(full_span_seqs, full_span_quals,
                                                                  start_only_seqs, start_only_quals,
                                                                  end_only_seqs, end_only_quals,
                                                                  scoring_scheme)[0]
            consensus_time = time.time() - consensus_start_time

            output += '  consensus sequence:      ' + \
                      int_to_str(len(self.consensus_sequence)) + ' bp '
            output += '(' + float_to_str(consensus_time, 2) + ' sec)\n'
            if verbosity > 3:
                output += '                           ' + self.consensus_sequence + '\n'

            target_path_length = len(self.consensus_sequence) + (2 * self.graph.overlap)
            mean_overlap = 0

        # For full spans without sequence, we simply need a mean distance.
        else:
            self.consensus_sequence = ''
            mean_overlap = int(round(sum(x[0] for x in full_spans_without_seq) /
                                     len(full_spans_without_seq)))
            output += '  mean overlap:            ' + int_to_str(abs(mean_overlap)) + '\n'
            target_path_length = mean_overlap + (2 * self.graph.overlap)

        output += '  target path length:      ' + int_to_str(target_path_length) + ' bp\n'

        path_start_time = time.time()
        all_paths, self.all_best_paths = self.graph.get_best_paths_for_seq(self.start_segment,
                                                                           self.end_segment,
                                                                           target_path_length,
                                                                           self.consensus_sequence,
                                                                           scoring_scheme)
        path_time = time.time() - path_start_time

        output += '  path count:              ' + int_to_str(len(all_paths)) + ' '
        output += '(' + float_to_str(path_time, 2) + ' sec)\n'
        if verbosity > 2:
            for i, path in enumerate(all_paths):
                label = '  path ' + str(i + 1) + ':'
                label = label.ljust(27)
                output += label + ', '.join(str(x) for x in path[0])
                output += ' (' + int_to_str(self.graph.get_path_length(path[0])) + ' bp, '
                output += 'raw score = ' + float_to_str(path[1], 1) + ', '
                output += 'scaled score = ' + float_to_str(path[3], 2) + ', '
                output += 'length discrepancy = ' + int_to_str(path[2]) + ' bp)\n'

        # If paths were found, use a path sequence for the bridge.
        if self.all_best_paths:

            if len(self.all_best_paths) == 1:
                output += '  best path:               '
            else:
                output += '  best paths:              '
            for i, best_path in enumerate(self.all_best_paths):
                if i > 0:
                    output += '                           '
                if best_path[0]:
                    output += ', '.join(str(x) for x in best_path[0])
                else:
                    output += 'direct connection'
                output += ' (' + int_to_str(self.graph.get_path_length(best_path[0])) + ' bp, '
                output += 'raw score = ' + float_to_str(best_path[1], 1) + ', '
                output += 'scaled score = ' + float_to_str(best_path[3], 2) + ', '
                output += 'length discrepancy = ' + int_to_str(best_path[2]) + ' bp)\n'

            for best_path in self.all_best_paths:
                path = best_path[0]
                if path:
                    path_sequence = self.graph.get_path_sequence(path)
                else:  # No path implies a direction connection between the start and end segments.
                    start_overlap = \
                        self.graph.seq_from_signed_seg_num(self.start_segment)[-self.graph.overlap:]
                    end_overlap = \
                        self.graph.seq_from_signed_seg_num(self.end_segment)[:self.graph.overlap]
                    path_sequence = start_overlap + end_overlap
                self.all_bridge_sequences.append(path_sequence)

            self.graph_path = self.all_best_paths[0][0]
            self.bridge_sequence = self.all_bridge_sequences[0]
            self.path_support = True

            # Start the quality on the scaled score of the path alignment (which maxes at 100.0).
            self.quality = self.all_best_paths[0][3]

        # If a path wasn't found, the consensus sequence is the bridge (with the overlaps added).
        else:
            self.graph_path = []
            self.path_support = False
            output += '  best path:               none found\n'
            start_overlap = \
                self.graph.seq_from_signed_seg_num(self.start_segment)[-self.graph.overlap:]
            end_overlap = \
                self.graph.seq_from_signed_seg_num(self.end_segment)[:self.graph.overlap]

            # If there is a consensus sequence, we simply tack the overlaps onto its ends.
            if self.consensus_sequence:
                self.bridge_sequence = start_overlap + self.consensus_sequence + end_overlap

            # If there is the consensus sequence is exactly zero, then the two segments butt up
            # exactly.
            elif mean_overlap == 0:
                self.bridge_sequence = start_overlap + end_overlap

            # If the consensus sequence is negative, that implies the start and the end segments
            # overlap. In this case we need to find the exact overlap. If we don't find it, we
            # just glue the sequences together, like we did for mean_overlap == 0.
            else:
                larger_overlaps = list(range(abs(mean_overlap), self.graph.overlap))
                smaller_overlaps = list(range(abs(mean_overlap) - 1, 0, -1))
                test_overlaps = []
                for i in range(max(len(larger_overlaps), len(smaller_overlaps))):
                    if i < len(larger_overlaps):
                        test_overlaps.append(larger_overlaps[i])
                    if i < len(smaller_overlaps):
                        test_overlaps.append(smaller_overlaps[i])
                for test_overlap in test_overlaps:
                    if start_overlap[-test_overlap:] == end_overlap[:test_overlap]:
                        actual_overlap = test_overlap
                        break
                else:
                    actual_overlap = 0
                self.bridge_sequence = start_overlap + end_overlap[actual_overlap:]

            # The quality of non-graph-path-supported bridges depends on the dead ends. If this
            # bridge is connecting two dead ends, then it isn't penalised at all (because there
            # couldn't be a graph path). If only one of the two segments is a dead end, then it gets
            # a small penalty. If neither are a dead end, then it gets a large penalty (because in
            # this case we'd expect a graph path).
            self.quality = 100.0
            if not self.graph.ends_with_dead_end(self.start_segment):
                self.quality -= 40.0
            if not self.graph.starts_with_dead_end(self.end_segment):
                self.quality -= 40.0

        # Expected read count is determined using the read lengths and bridge size. For a given
        # read length and bridge, there are an estimable number of positions where a read of that
        # length would be able to contribute to the bridge. This is used to get the probability
        # that any read would create a bridge, and totalling those up gives us our estimated count.
        bridge_len = max(len(self.bridge_sequence), self.graph.overlap)
        min_read_len = (2 * min_alignment_length) + bridge_len - (2 * self.graph.overlap)
        total_possible_placements = 0
        for read_len, count in read_lengths.items():
            if read_len < min_read_len:
                continue
            possible_read_placements = read_len - min_read_len + 1
            possible_read_placements *= count
            possible_read_placements *= max(self.depth, 1)
            total_possible_placements += possible_read_placements
        expected_read_count = total_possible_placements / estimated_genome_size
        actual_read_count = len(self.full_span_reads)

        # Adjust the expected read count down, especially for higher values.
        # TO DO: reevaluate this step - is it necessary?
        expected_read_count = reduce_expected_count(expected_read_count, 30, 0.5)

        if verbosity > 2:
            output += '  expected bridging reads: ' + float_to_str(expected_read_count, 2) + '\n'
            output += '  starting quality:        ' + float_to_str(self.quality, 2) + '\n'

        # The start segment and end segment should agree in depth. If they don't, that's very bad,
        # as it implies that they aren't actually single-copy or on the same piece of DNA.
        depth_agreement_factor = get_num_agreement(start_seg.depth, end_seg.depth) ** 2
        self.quality *= depth_agreement_factor

        # The number of reads which contribute to a bridge is a big deal, so the read count factor
        # scales linearly. This value is capped at 1, which means that bridges with too few reads
        # are punished but bridges with excess reads are not rewarded.
        read_count_factor = min(1.0, actual_read_count / expected_read_count)
        self.quality *= read_count_factor

        # The length of alignments to the start/end segments is positively correlated with quality
        # to reward bridges with long alignments. Specifically, we want there to be at least one
        # spanning read with a long alignment to the start segment and at least one spanning read
        # with a long alignment to the end segment.
        longest_start_alignment = max(x[2].get_aligned_ref_length() for x in self.full_span_reads)
        longest_end_alignment = max(x[3].get_aligned_ref_length() for x in self.full_span_reads)
        alignment_length = min(longest_start_alignment, longest_end_alignment)
        align_length_factor = score_function(alignment_length, min_alignment_length * 4)
        self.quality *= align_length_factor

        # The mean alignment score to the start/end segments is positively correlated with quality,
        # so bridges with high quality alignments are rewarded. Specifically, we want there to be at
        # least one spanning read with a high quality alignment to the start segment and at least
        # one spanning read with a high quality alignment to the end segment.
        best_start_alignment = max(x[2].scaled_score for x in self.full_span_reads)
        best_end_alignment = max(x[3].scaled_score for x in self.full_span_reads)
        alignment_quality = min(best_start_alignment, best_end_alignment)
        align_score_factor = alignment_quality / 100.0
        self.quality *= align_score_factor

        # Bridges between long start/end segments are rewarded, as they are more likely to actually
        # be single-copy. We apply a length factor for both the start and the end segments,
        # and then apply the smaller of two again. This is to punish cases where both segments
        # are not long.
        start_length_factor = score_function(start_seg.get_length(), min_alignment_length * 4)
        self.quality *= start_length_factor
        end_length_factor = score_function(end_seg.get_length(), min_alignment_length * 4)
        self.quality *= end_length_factor
        smaller_length_factor = min(start_length_factor, end_length_factor)
        self.quality *= smaller_length_factor

        if verbosity > 2:
            output += '  depth agreement factor:  ' + float_to_str(depth_agreement_factor, 2) + '\n'
            output += '  read count factor:       ' + float_to_str(read_count_factor, 2) + '\n'
            output += '  alignment length factor: ' + float_to_str(align_length_factor, 2) + '\n'
            output += '  alignment score factor:  ' + float_to_str(align_score_factor, 2) + '\n'
            output += '  start length factor:     ' + float_to_str(start_length_factor, 2) + '\n'
            output += '  end length factor:       ' + float_to_str(end_length_factor, 2) + '\n'
            output += '  smaller_length_factor:   ' + float_to_str(smaller_length_factor, 2) + '\n'
            output += '  final quality:           ' + float_to_str(self.quality, 2) + '\n'
        if verbosity == 2:
            output += '  quality:                 ' + float_to_str(self.quality, 2) + '\n'
        return output

    def set_path_based_on_availability(self, graph):
        """
        If this bridge has multiple equally-good paths, this will change its path to the one most
        available in the graph.
        """
        best_path = self.all_best_paths[0][0]
        best_availability = graph.get_path_availability(best_path)
        best_sequence = self.all_bridge_sequences[0]
        for i in range(1, len(self.all_best_paths)):
            potential_path = self.all_best_paths[i][0]
            potential_path_availability = graph.get_path_availability(potential_path)
            if potential_path_availability > best_availability:
                best_availability = potential_path_availability
                best_path = potential_path
                best_sequence = self.all_bridge_sequences[i]
        self.graph_path = best_path
        self.bridge_sequence = best_sequence

    def contains_full_span_sequence(self):
        """
        Some LongReadBridge objects bridge two close segments, and therefore do not actually have
        any bridging read sequence (just a negative number implying overlap). This function returns
        True if any full span read sequences exist and False if they are just negative numbers.
        """
        for full_span_read in self.full_span_reads:
            seq_or_num = full_span_read[0]
            if not isinstance(seq_or_num, int):
                return True
        return False

    @staticmethod
    def get_type_score():
        """
        Returns a score indicating the relative importance of the bridge types:
        LongReadBridge = 2, SpadesContigBridge = 1, LoopUnrollingBridge = 0
        """
        return 2


class LoopUnrollingBridge(object):
    """
    This class describes a bridge created from unrolling an assembly graph loop.

    Quality is affected by:
      * How well the start and end segments' depths agree.
      * How close the determined loop count is to a whole number.
      * The final loop count (higher counts get lower quality).
    """

    def __init__(self, graph, start, end, middle, repeat):
        """
        This constructor assumes the the start, end, middle and repeat segments form a simple loop
        in the graph supported by either a SPAdes contig or a long read alignment. It will use
        segment depths to determine the loop count and score the bridge's quality.
        """
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
        self.quality *= (depth_agreement * depth_agreement)  # has squared effect on quality

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

    @staticmethod
    def get_type_score():
        """
        Returns a score indicating the relative importance of the bridge types:
        LongReadBridge = 2, SpadesContigBridge = 1, LoopUnrollingBridge = 0
        """
        return 0


def create_spades_contig_bridges(graph, single_copy_segments):
    """
    Builds graph bridges using the SPAdes contig paths.
    """
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
    conflicting_paths_no_dupes = []
    for path in conflicting_paths:
        if path not in conflicting_paths_no_dupes:
            conflicting_paths_no_dupes.append(path)
    conflicting_paths = conflicting_paths_no_dupes
    final_bridge_paths = [x for x in bridge_path_list if x not in conflicting_paths]

    return [SpadesContigBridge(spades_contig_path=x, graph=graph) for x in final_bridge_paths]


def find_contig_bridges(segment_num, path, single_copy_numbers):
    """
    This function returns a list of lists: every part of the path which starts on the segment_num
    and ends on any of the single_copy_numbers.
    """
    bridge_paths = []
    indices = [i for i, x in enumerate(path) if abs(x) == segment_num]
    for index in indices:
        bridge_path = [path[index]]
        for i in range(index + 1, len(path)):
            bridge_path.append(path[i])
            if path[i] in single_copy_numbers or -path[i] in single_copy_numbers:
                break
        else:
            bridge_path = []
        if bridge_path:
            bridge_paths.append(bridge_path)
    return bridge_paths


def create_loop_unrolling_bridges(graph):
    """
    This function creates loop unrolling bridges using the information in SPAdes paths.
    """
    bridges = []
    simple_loops = graph.find_all_simple_loops()

    # A simple loop can either be caused by a repeat in one sequence (probably more typical) or by
    # a separate circular sequence which has some common sequence (less typical, but still very
    # possible: plasmids). We only want to unroll the former group, so we look for cases where the
    # loop's start or end is in a SPAdes contig path along with the middle. That implies that they
    # are on the same piece of DNA and can be unrolled.
    for start, end, middle, repeat in simple_loops:
        joined = False
        for path in graph.paths.values():
            flipped_path = [-x for x in reversed(path)]
            if (start in path and middle in path) or \
                    (end in path and middle in path) or \
                    (start in flipped_path and middle in flipped_path) or \
                    (end in flipped_path and middle in flipped_path):
                joined = True
                break

        # If we've found evidence the simply loop is a single piece of DNA, then we'll make a loop
        # unrolling bridge!
        if joined:
            bridges.append(LoopUnrollingBridge(graph, start, end, middle, repeat))

    return bridges


def create_long_read_bridges(graph, read_dict, read_names, single_copy_segments, verbosity,
                             existing_bridges, min_scaled_score, threads, scoring_scheme,
                             min_alignment_length):
    """
    Makes bridges between single-copy segments using the alignments in the long reads.
    """
    single_copy_seg_num_set = set()
    for seg in single_copy_segments:
        single_copy_seg_num_set.add(seg.number)

    # This dictionary will collect the read sequences which span between two single-copy segments.
    # These are the most useful sequences and will be used to either create a new bridge or enhance
    # an existing bridge.
    # Key = tuple of signed segment numbers (the segments being bridged)
    # Value = list of tuples containing the bridging sequence and the single-copy segment
    #         alignments.
    spanning_read_seqs = defaultdict(list)

    # This dictionary will collect all of the read sequences which don't manage to span between
    # two single-copy segments, but they do overlap one single-copy segment and can therefore be
    # useful for generating consensus sequences in a bridge.
    # Key = signed segment number, Value = list of tuples containing the sequence and alignment
    overlapping_read_seqs = defaultdict(list)

    for read_name in read_names:
        read = read_dict[read_name]
        alignments = get_single_copy_alignments(read, single_copy_seg_num_set, min_scaled_score)
        if not alignments:
            continue

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
            if len(available_alignments) < 2:
                continue

            for i in range(len(available_alignments)):
                if i < len(available_alignments) - 1:
                    alignment_1 = available_alignments[i]
                    alignment_2 = available_alignments[i + 1]

                # Special case: when the first and last alignments are to the same graph segment,
                # make a bridge for them, even if they aren't a particularly high scoring pair of
                # alignments. This can help to circularise plasmids which are very tied up with
                # other, similar plasmids.
                elif available_alignments[0].ref.name == available_alignments[-1].ref.name:
                    alignment_1 = available_alignments[0]
                    alignment_2 = available_alignments[-1]
                else:
                    continue

                # Standardise the order so we don't end up with both directions (e.g. 5 to -6 and
                # 6 to -5) in spanning_read_seqs.
                seg_nums, flipped = flip_number_order(alignment_1.get_signed_ref_num(),
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
                        bridge_seq = bridge_end - bridge_start  # 0 or a negative number
                        bridge_qual = ''

                    spanning_read_seqs[seg_nums].append((bridge_seq, bridge_qual, alignment_1,
                                                         alignment_2))
                    already_added.add(seg_nums)

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
        last_alignment = available_alignments[-1]
        end_overlap, end_qual = last_alignment.get_end_overlapping_read_seq()
        if end_overlap:
            seg_num = last_alignment.get_signed_ref_num()
            overlapping_read_seqs[seg_num].append((end_overlap, end_qual, last_alignment))

    # If a bridge already exists for a spanning sequence, we add the sequence to the bridge. If
    # not, we create a new bridge and add it.
    new_bridges = []
    for seg_nums, span in spanning_read_seqs.items():
        start, end = seg_nums

        # If the start and end are the same and already exclusively connect (i.e. if this segment
        # is already circular), then skip - there's no need to bridge.
        if start == end and graph.get_downstream_seg_nums(start) == [start] and \
                graph.get_upstream_seg_nums(start) == [start]:
            continue

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

    # During finalisation, we will compare the expected read count to the actual read count for
    # each bridge. To do this, we'll need the lengths of all reads (excluding those with no
    # alignments). We also need an estimate of the genome size.
    read_lengths = defaultdict(int)
    for read_name in read_names:
        read = read_dict[read_name]
        if read.alignments:
            read_lengths[read.get_length()] += 1
    estimated_genome_size = graph.get_estimated_sequence_len()

    # Now we need to finalise the bridges. This is the intensive step, as it involves creating a
    # consensus sequence, finding graph paths and doing alignments between the consensus and the
    # graph paths. We can therefore use threads to make this faster.
    long_read_bridges = [x for x in all_bridges if isinstance(x, LongReadBridge)]
    num_long_read_bridges = len(long_read_bridges)
    if verbosity == 1:
        print_progress_line(0, num_long_read_bridges, prefix='Bridge: ')
    completed_count = 0
    if threads == 1:
        for bridge in long_read_bridges:
            output = bridge.finalise(scoring_scheme, min_alignment_length, read_lengths,
                                     estimated_genome_size, verbosity)
            completed_count += 1
            if verbosity == 1:
                print_progress_line(completed_count, num_long_read_bridges, prefix='Bridge: ')
            if verbosity > 1:
                print(output, end='', flush=True)
    else:
        pool = ThreadPool(threads)
        arg_list = []

        # Sort the bridges based on their size, with longer spanning sequences coming first. This
        # will make the big ones runs first which helps to more efficiently use the CPU cores.
        # E.g. if the biggest bridge was at the end, we'd be left waiting for it to finish with
        # only one core (bad), but if it was at the start, other work could be done in parallel.
        long_read_bridges = sorted(long_read_bridges, reverse=True,
                                   key=lambda x: x.get_total_full_span_seq_length())

        for bridge in long_read_bridges:
            arg_list.append((bridge, scoring_scheme, min_alignment_length, read_lengths,
                             estimated_genome_size, verbosity))

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
                print(output, end='', flush=True)

    if verbosity == 1:
        print()

    return all_bridges


def get_mean_depth(seg_1, seg_2, graph):
    """
    Returns the mean depth of the two segments, weighted by their length.
    """
    return weighted_average(seg_1.depth, seg_2.depth,
                            seg_1.get_length_no_overlap(graph.overlap),
                            seg_2.get_length_no_overlap(graph.overlap))


def path_is_self_contained(path, start, end, graph):
    """
    Returns True if the path segments are only connected to each other and the start/end segments.
    If they are connected to anything else, it returns False.
    """
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


def get_single_copy_alignments(read, single_copy_num_set, min_scaled_score):
    """
    Returns a list of single-copy segment alignments for the read.
    """
    sc_alignments = []
    for alignment in read.alignments:
        if alignment.ref.number in single_copy_num_set and \
                        alignment.scaled_score >= min_scaled_score:
            sc_alignments.append(alignment)
    return sc_alignments


def finalise_bridge(all_args):
    """
    Just a one-argument version of bridge.finalise, for pool.imap.
    """
    bridge, scoring_scheme, min_alignment_length, read_lengths, estimated_genome_size, verbosity = \
        all_args
    return bridge.finalise(scoring_scheme, min_alignment_length, read_lengths,
                           estimated_genome_size, verbosity)


def get_bridge_str(start, middle, end):
    """
    Returns a bridge sequence in human-readable form.
    """
    bridge_str = str(start) + ' -> '
    if middle:
        bridge_str += ', '.join([str(x) for x in middle]) + ' -> '
    bridge_str += str(end)
    return bridge_str


def reduce_expected_count(expected_count, a, b):
    """
    This function reduces the expected read count. It reduces by a factor which is a function of
    the read count, so low expected values aren't reduced much, but high expected values are
    reduced more. This is to help with high read depth cases where expected counts get quite high.

    https://www.desmos.com/calculator
    y=x\cdot \left(\left(\frac{a}{a+x}\right)\cdot \left(1-b\right)+b\right)
    """
    return expected_count * ((a / (a + expected_count)) * (1.0 - b) + b)
