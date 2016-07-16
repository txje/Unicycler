"""
This module describes an assembly graph and many related functions.

Author: Ryan Wick
email: rrwick@gmail.com
"""

import math
from collections import deque, defaultdict
from .misc import int_to_str, float_to_str, weighted_average, weighted_average_list, \
    print_section_header, get_num_agreement, reverse_complement
from .bridge import SpadesContigBridge, LoopUnrollingBridge, LongReadBridge, \
    get_applicable_bridge_pieces, get_bridge_str
from .cpp_function_wrappers import fully_global_alignment, path_alignment


class TooManyPaths(Exception):
    pass


class CannotTrimOverlaps(Exception):
    pass


class AssemblyGraph(object):
    """
    This class holds an assembly graph with segments and links.
    """

    def __init__(self, filename, overlap, paths_file=None,
                 insert_size_mean=250, insert_size_deviation=50):
        self.segments = {}  # Dict of unsigned segment number -> segment
        self.forward_links = {}  # Dict of signed segment number -> list of signed segment numbers
        self.reverse_links = {}  # Dict of signed segment number <- list of signed segment numbers
        self.copy_depths = {}  # Dict of unsigned segment number -> list of copy depths
        self.paths = {}  # Dict of path name -> list of signed segment numbers
        self.overlap = overlap
        self.insert_size_mean = insert_size_mean
        self.insert_size_deviation = insert_size_deviation

        if filename.endswith('.fastg'):
            self.load_from_fastg(filename)
        else:
            self.load_from_gfa(filename)
            if not overlap:
                self.overlap = get_overlap_from_gfa_link(filename)

        if paths_file:
            self.load_spades_paths(paths_file)

    def load_from_fastg(self, filename):
        """
        Loads a Graph from a SPAdes-style FASTG file.
        """
        # Load in the graph segments.
        headers, sequences = get_headers_and_sequences(filename)
        for i, header in enumerate(headers):
            num = get_unsigned_number_from_header(header)
            sequence = sequences[i]
            positive = is_header_positive(header)

            # If the segment already exists, then add this sequence.
            if num in self.segments:
                self.segments[num].add_sequence(sequence, positive)

            # If the segment does not exist, make it.
            else:
                depth = get_depth_from_header(header)
                segment = Segment(num, depth, sequence, positive)
                self.segments[num] = segment

        # Make sure that every segment has both a forward and reverse sequence.
        for segment in self.segments.values():
            segment.build_other_sequence_if_necessary()

        # Load in the links.
        for header in headers:
            start, end_list = get_links_from_header(header)
            if end_list:
                self.forward_links[start] = end_list
        self.forward_links = build_rc_links_if_necessary(self.forward_links)
        self.reverse_links = build_reverse_links(self.forward_links)
        self.sort_link_order()

    def load_from_gfa(self, filename):
        """
        Loads a Graph from a GFA file. It does not load any GFA file, but makes some restrictions:
        1) The segment names must be integers.
        2) The depths should be stored in a DP tag.
        3) All link overlaps are the same (equal to the graph overlap value).
        """
        # Load in the segments.
        gfa_file = open(filename, 'rt')
        for line in gfa_file:
            if line.startswith('S'):
                line_parts = line.strip().split('\t')
                num = int(line_parts[1])
                depth = 1.0
                for part in line_parts:
                    if part.startswith('DP:'):
                        depth = float(part[5:])
                sequence = line_parts[2]
                self.segments[num] = Segment(num, depth, sequence, True)
                self.segments[num].build_other_sequence_if_necessary()
            if line.startswith('i'):
                line_parts = line.strip().split('\t')
                try:
                    self.insert_size_mean = float(line_parts[1])
                    self.insert_size_deviation = float(line_parts[2])
                except ValueError:
                    pass
        gfa_file.close()

        # Load in the links.
        gfa_file = open(filename, 'rt')
        for line in gfa_file:
            if line.startswith('L'):
                line_parts = line.strip().split('\t')
                start = signed_string_to_int(line_parts[1] + line_parts[2])
                end = signed_string_to_int(line_parts[3] + line_parts[4])
                if start not in self.forward_links:
                    self.forward_links[start] = [end]
                else:
                    self.forward_links[start].append(end)
        self.forward_links = build_rc_links_if_necessary(self.forward_links)
        self.reverse_links = build_reverse_links(self.forward_links)
        gfa_file.close()
        self.sort_link_order()

        # Load in the paths
        gfa_file = open(filename, 'rt')
        for line in gfa_file:
            if line.startswith('P'):
                line_parts = line.strip().split('\t')
                path_name = line_parts[1]
                segments = [signed_string_to_int(x) for x in line_parts[2].split(',')]
                self.paths[path_name] = segments
        gfa_file.close()

    def load_spades_paths(self, filename):
        """
        Loads in SPAdes contig paths from file.
        It only saves the positive paths and does not save paths with only one segment.
        If a SPAdes path has a gap (semicolon), then it treats each component part as a separate
        path (i.e. paths do not span gaps).
        """
        names = []
        segment_strings = []
        name = ''
        segment_string = ''

        paths_file = open(filename, 'rt')
        for line in paths_file:
            line = line.strip()
            if not line:
                continue
            if line.startswith('NODE'):
                if name:
                    names.append(name)
                    segment_strings.append(segment_string)
                    segment_string = ''
                name = line
            else:
                segment_string += line
        paths_file.close()
        if name:
            names.append(name)
            segment_strings.append(segment_string)

        for i, name in enumerate(names):
            if name.endswith("'"):
                continue
            name_parts = name.split('_')
            if len(name_parts) < 2:
                continue
            name = '_'.join(name_parts[:2])
            segment_string = segment_strings[i]
            if not segment_string:
                continue
            segment_string_parts = segment_string.split(';')
            segment_string_parts = [x for x in segment_string_parts if len(x.split(',')) > 1]
            for j, segment_string_part in enumerate(segment_string_parts):
                path_name = name
                if len(segment_string_parts) > 1:
                    path_name += '_' + str(j + 1)
                segments = [signed_string_to_int(x) for x in segment_string_part.split(',')]
                self.paths[path_name] = segments

    def get_median_read_depth(self, segment_list=None):
        """
        Returns the assembly graph's median read depth (by base).  Optionally, a list of segments
        can be given, in which case only those segments are used for the calculation.
        """
        if not segment_list:
            segment_list = self.segments.values()
        sorted_segments = sorted(segment_list, key=lambda x: x.depth)
        total_length = 0
        for segment in sorted_segments:
            total_length += segment.get_length_no_overlap(self.overlap)
        halfway_length = total_length // 2
        length_so_far = 0
        for segment in sorted_segments:
            length_so_far += segment.get_length_no_overlap(self.overlap)
            if length_so_far >= halfway_length:
                return segment.depth
        return 0.0

    def normalise_read_depths(self):
        """
        For every segment in the graph, divide its depth by the graph's median.
        This makes segments with the median depth have a depth of 1, segments with more than the
        median a depth of greater than 1 and segments with less than the median a depth of less
        than 1.
        """
        median_depth = self.get_median_read_depth()
        for segment in self.segments.values():
            segment.divide_depth(median_depth)

    def get_total_length(self):
        """
        Returns the sum of all segment sequence lengths.
        """
        return sum([x.get_length() for x in self.segments.values()])

    def get_total_length_no_overlaps(self):
        """
        Returns the sum of all segment sequence lengths, subtracting the overlap size from each
        segment.
        """
        return sum([x.get_length_no_overlap(self.overlap) for x in self.segments.values()])

    def save_to_fasta(self, filename, verbosity=0):
        """
        Saves whole graph (only forward sequences) to a FASTA file.
        """
        fasta = open(filename, 'w')
        if verbosity > 0:
            print('\nSaving', filename)
        sorted_segments = sorted(self.segments.values(), key=lambda x: x.number)
        for segment in sorted_segments:
            fasta.write('>' + str(segment.number) + '\n')
            fasta.write(add_line_breaks_to_sequence(segment.forward_sequence, 60))

    def save_to_fastg(self, filename, verbosity=0):
        """
        Saves whole graph to a SPAdes-style FASTG file.
        """
        fastg = open(filename, 'w')
        if verbosity > 0:
            print('\nSaving', filename)
        sorted_segments = sorted(self.segments.values(), key=lambda x: x.number)
        for segment in sorted_segments:
            fastg.write(self.get_fastg_header_with_links(segment, True))
            fastg.write(add_line_breaks_to_sequence(segment.forward_sequence, 60))
            fastg.write(self.get_fastg_header_with_links(segment, False))
            fastg.write(add_line_breaks_to_sequence(segment.reverse_sequence, 60))

    def save_to_gfa(self, filename, verbosity, save_copy_depth_info=False,
                    save_seg_type_info=False, single_copy_segments=None):
        """
        Saves whole graph to a GFA file.
        """
        gfa = open(filename, 'w')
        if verbosity > 0:
            print('\nSaving', filename)
        sorted_segments = sorted(self.segments.values(), key=lambda x: x.number)
        if single_copy_segments is not None:
            single_copy_segments = set([x.number for x in single_copy_segments])
        for segment in sorted_segments:
            segment_line = segment.gfa_segment_line()
            if save_copy_depth_info and segment.number in self.copy_depths:
                segment_line = segment_line[:-1]  # Remove newline
                segment_line += '\tLB:z:' + self.get_depth_string(segment)
                segment_line += '\tCL:z:' + self.get_copy_number_colour(segment)
                segment_line += '\n'
            if save_seg_type_info:
                segment_line = segment_line[:-1]  # Remove newline
                segment_line += '\tLB:z:' + segment.get_seg_type_label()
                if segment.number in single_copy_segments:
                    colour = 'forestgreen'
                elif segment.bridge is None:
                    colour = 'grey'
                else:
                    colour = 'pink'
                segment_line += '\tCL:z:' + colour
                segment_line += '\n'

            gfa.write(segment_line)
        gfa.write(self.get_all_gfa_link_lines())
        paths = sorted(self.paths.items())
        overlap_cigar = str(self.overlap) + 'M'
        for path_name, segment_list in paths:
            gfa.write('P\t' + path_name + '\t')
            gfa.write(','.join([int_to_signed_string(x) for x in segment_list]))
            gfa.write('\t')
            gfa.write(','.join([overlap_cigar] * (len(segment_list) - 1)))
            gfa.write('\n')
        if self.insert_size_mean is not None and self.insert_size_deviation is not None:
            gfa.write('i\t')
            gfa.write(str(self.insert_size_mean))
            gfa.write('\t')
            gfa.write(str(self.insert_size_deviation))
            gfa.write('\n')
        gfa.close()

    def get_all_gfa_link_lines(self):
        """
        Returns a string of the link component of the GFA file for this graph.
        """
        gfa_link_lines = ''
        for start, ends in self.forward_links.items():
            for end in ends:
                if is_link_positive(start, end):
                    gfa_link_lines += self.gfa_link_line(start, end)
        return gfa_link_lines

    def get_fastg_header_with_links(self, segment, positive):
        """
        Returns a full SPAdes-style FASTG header for a segment, including the leading '>', all of
        the links, the trailing ';' and a newline.
        """
        number = segment.number
        if not positive:
            number *= -1
        header = '>' + segment.get_fastg_header(positive)
        if number in self.forward_links:
            header += ':'
            next_segment_headers = []
            for next_num in self.forward_links[number]:
                if next_num < 0:
                    next_positive = False
                    next_num *= -1
                else:
                    next_positive = True
                next_segment = self.segments[next_num]
                next_segment_headers.append(next_segment.get_fastg_header(next_positive))
            header += ','.join(next_segment_headers)
        header += ';\n'
        return header

    def total_dead_end_count(self):
        """
        Returns the total number of dead ends in the assembly graph.
        """
        dead_ends = 0
        for seg_num in self.segments:
            dead_ends += self.dead_end_count(seg_num)
        return dead_ends

    def dead_end_count(self, seg_num):
        """
        Returns the number of dead ends for one segment: 0, 1 or 2.
        """
        dead_ends = 0
        if seg_num not in self.forward_links or not self.forward_links[seg_num]:
            dead_ends += 1
        if seg_num not in self.reverse_links or not self.reverse_links[seg_num]:
            dead_ends += 1
        return dead_ends

    def filter_by_read_depth(self, relative_depth_cutoff):
        """
        This function removes segments from the graph based on a relative depth cutoff. Segments
        are considered below the cutoff if they are less than the cutoff for the entire graph or
        less than the cutoff for their connected component.
        To be removed, one of the following must also be true:
          1) the segment has at least one dead end
          2) the segment is part of a connected component where all of the segments are below the
             whole graph cutoff
          3) deleting the segment would not create any dead ends
        """
        segment_nums_to_remove = []
        whole_graph_cutoff = self.get_median_read_depth() * relative_depth_cutoff
        connected_components = self.get_connected_components()
        for component in connected_components:
            component_segs = [self.segments[x] for x in component]
            component_cutoff = self.get_median_read_depth(component_segs) * relative_depth_cutoff
            for seg_num in component:
                segment = self.segments[seg_num]
                if segment.depth < whole_graph_cutoff or segment.depth < component_cutoff:
                    if self.dead_end_count(seg_num) > 0 or \
                            self.all_segments_below_depth(component, whole_graph_cutoff) or \
                            self.dead_end_change_if_deleted(seg_num) <= 0:
                        segment_nums_to_remove.append(seg_num)
        self.remove_segments(segment_nums_to_remove)

    def filter_homopolymer_loops(self):
        """
        A common feature in SPAdes graphs is a small piece of the graph (often just one segment)
        which has nothing but one base.  Filter these out.
        """
        segment_nums_to_remove = []
        connected_components = self.get_connected_components()
        for component_nums in connected_components:
            component_segments = [self.segments[x] for x in component_nums]
            if all_segments_are_one_base(component_segments):
                segment_nums_to_remove += component_nums
        self.remove_segments(segment_nums_to_remove)

    def remove_segments(self, nums_to_remove):
        """
        Given a list of segment numbers to remove, this function rebuilds the graph's segments
        and links, excluding those segments. It also deletes any paths which contain those
        segments.
        """
        new_segments = {}
        for num, segment in self.segments.items():
            if num not in nums_to_remove:
                new_segments[num] = segment
        self.segments = new_segments

        for num in nums_to_remove:
            if num in self.copy_depths:
                del self.copy_depths[num]

        self.forward_links = remove_nums_from_links(self.forward_links, nums_to_remove)
        self.reverse_links = remove_nums_from_links(self.reverse_links, nums_to_remove)

        paths_to_delete = set()
        neg_nums_to_remove = [-x for x in nums_to_remove]
        for path_name, path_nums in self.paths.items():
            if len(list(set(nums_to_remove) & set(path_nums))) > 0:
                paths_to_delete.add(path_name)
            if len(list(set(neg_nums_to_remove) & set(path_nums))) > 0:
                paths_to_delete.add(path_name)
        for path_to_delete in paths_to_delete:
            del self.paths[path_to_delete]

    def remove_small_components(self, min_component_size, verbosity):
        """
        Remove small graph components, but only if they do not contain any bridges. The idea is
        to clean up parts of the graph that were orphaned by the bridging process. But if they
        contain a bridge, then they are more likely to be genuine and we keep them.
        """
        segment_nums_to_remove = []
        connected_components = self.get_connected_components()
        for component_nums in connected_components:
            component_segments = [self.segments[x] for x in component_nums]
            component_length = sum(x.get_length() for x in component_segments)
            if component_length >= min_component_size:
                continue
            if any(x.bridge is not None for x in component_segments):
                continue
            segment_nums_to_remove += component_nums
        self.remove_segments(segment_nums_to_remove)
        if verbosity > 1 and segment_nums_to_remove:
            print('\nRemoved small components:', ', '.join(str(x) for x in segment_nums_to_remove))

    def remove_small_dead_ends(self, min_dead_end_size, verbosity):
        """
        Remove small segments which are graph dead-ends. This is just to tidy things up a bit
        before the final merge.
        """
        removed_segments = []
        while True:
            for seg_num, segment in self.segments.items():
                if segment.get_length() >= min_dead_end_size:
                    continue
                if self.dead_end_change_if_deleted(seg_num) < 0:
                    self.remove_segments([seg_num])
                    removed_segments.append(seg_num)
                    break
            else:
                break
        if verbosity > 1 and removed_segments:
            print('\nRemoved small dead ends: ', ', '.join(str(x) for x in removed_segments))

    def merge_all_possible(self):
        """
        This function merges segments which are in a simple, unbranching path.
        """
        while True:
            # Sort the segment numbers first so we apply the merging in a consistent order.
            seg_nums = sorted(list(self.segments.keys()))
            for num in seg_nums:
                path = self.get_simple_path(num)
                if len(path) <= 1:
                    continue
                else:
                    self.merge_simple_path(path)
                    break
            else:
                break
        self.renumber_segments()

    def merge_simple_path(self, merge_path):
        """
        Merges the path into a single segment and adjusts any graph paths as necessary. Assumes
        that the path is a simple, unbranching path and can be merged.
        """
        start = merge_path[0]
        end = merge_path[-1]

        # The merged sequence depth is the weighted mean of the components.
        depths = [self.segments[abs(x)].depth for x in merge_path]
        lengths = [self.segments[abs(x)].get_length() - self.overlap for x in merge_path]
        if sum(lengths) > 0.0:
            mean_depth = weighted_average_list(depths, lengths)
        else:
            mean_depth = 1.0

        new_seg_num = self.get_next_available_seg_number()
        merged_forward_seq = self.get_path_sequence(merge_path)
        new_seg = Segment(new_seg_num, mean_depth, merged_forward_seq, True)
        new_seg.build_other_sequence_if_necessary()

        # Save some info that we'll need, and then delete the old segments.
        paths_copy = self.paths.copy()
        outgoing_links = []
        if end in self.forward_links:
            outgoing_links = self.forward_links[end]
        incoming_links = []
        if start in self.reverse_links:
            incoming_links = self.reverse_links[start]
        outgoing_links = find_replace_one_val_in_list(outgoing_links, start, new_seg_num)
        outgoing_links = find_replace_one_val_in_list(outgoing_links, -end, -new_seg_num)
        incoming_links = find_replace_one_val_in_list(incoming_links, end, new_seg_num)
        incoming_links = find_replace_one_val_in_list(incoming_links, -start, -new_seg_num)
        self.remove_segments([abs(x) for x in merge_path])

        # Add the new segment to the graph and give it the links from its source segments.
        self.segments[new_seg_num] = new_seg
        for link in outgoing_links:
            self.add_link(new_seg_num, link)
        for link in incoming_links:
            self.add_link(link, new_seg_num)

        # Merge the segments in any paths.
        flipped_merge_path = [-x for x in reversed(merge_path)]
        for path_name in paths_copy:
            paths_copy[path_name] = find_replace_in_list(paths_copy[path_name], merge_path,
                                                         [new_seg_num])
            paths_copy[path_name] = find_replace_in_list(paths_copy[path_name], flipped_merge_path,
                                                         [-new_seg_num])

        # If any paths still contain the original segments, then split those paths into pieces,
        # removing the original segments.
        new_paths = {}
        for path_name, path_segments in paths_copy.items():
            split_paths = split_path_multiple(path_segments, merge_path + flipped_merge_path)
            if len(split_paths) == 1:
                new_paths[path_name] = split_paths[0]
            elif len(split_paths) > 1:
                for i, path in enumerate(split_paths):
                    new_paths[path_name + '_' + str(i + 1)] = path
        self.paths = new_paths

    def add_link(self, start, end):
        """
        Adds a link to the graph in all necessary ways: forward and reverse, and for reverse
        complements too.
        """
        if start not in self.forward_links:
            self.forward_links[start] = []
        if end not in self.forward_links[start]:
            self.forward_links[start].append(end)

        if end not in self.reverse_links:
            self.reverse_links[end] = []
        if start not in self.reverse_links[end]:
            self.reverse_links[end].append(start)

        if -start not in self.reverse_links:
            self.reverse_links[-start] = []
        if -end not in self.reverse_links[-start]:
            self.reverse_links[-start].append(-end)

        if -end not in self.forward_links:
            self.forward_links[-end] = []
        if -start not in self.forward_links[-end]:
            self.forward_links[-end].append(-start)

    def remove_link(self, start, end):
        """
        Removes a link from the graph in all necessary ways: forward and reverse, and for reverse
        complements too.
        """
        if start in self.forward_links:
            self.forward_links[start].remove(end)
        if -end in self.forward_links:
            self.forward_links[-end].remove(-start)
        if end in self.reverse_links:
            self.reverse_links[end].remove(start)
        if -start in self.reverse_links:
            self.reverse_links[-start].remove(-end)

    def get_seq_from_signed_seg_num(self, signed_num):
        """
        Returns the forwards or reverse sequence of a segment, if the number is next_positive or
        negative, respectively. Assumes the segment number is in the graph.
        """
        if signed_num > 0:
            return self.segments[signed_num].forward_sequence
        else:
            return self.segments[-signed_num].reverse_sequence

    def get_connected_components(self):
        """
        Returns a list of lists, where each inner list is the segment numbers of one connected
        component of the graph.
        E.g. [[1, 2], [3, 4, 5]] would mean that segments 1 and 2 are in a connected component
        and segments 3, 4 and 5 are in another connected component.
        """
        visited = set()
        components = []
        for v in self.segments:
            if v not in visited:
                component = []
                q = deque()
                q.append(v)
                visited.add(v)
                while q:
                    w = q.popleft()
                    component.append(w)
                    connected_segments = self.get_connected_segments(w)
                    for k in connected_segments:
                        if k not in visited:
                            visited.add(k)
                            q.append(k)
                components.append(sorted(component))

        # Sort (just for consistency from one run to the next)
        return sorted(components)

    def get_connected_segments(self, segment_num):
        """
        Given a segment number, this function returns a list of all other segment numbers for
        segments that are directly connected.
        It only returns positive numbers (i.e. is not strand-specific).
        """
        connected_segments = set()
        if segment_num in self.forward_links:
            downstream_segments = self.forward_links[segment_num]
            for segment in downstream_segments:
                connected_segments.add(abs(segment))
        if segment_num in self.reverse_links:
            upstream_segments = self.reverse_links[segment_num]
            for segment in upstream_segments:
                connected_segments.add(abs(segment))
        return list(connected_segments)

    def all_segments_below_depth(self, segment_nums, cutoff):
        """
        Returns true if all segments in the list are below the depth cutoff.
        """
        for num in segment_nums:
            if self.segments[num].depth >= cutoff:
                return False
        return True

    def get_n_segment_length(self, n_percent):
        """
        Returns the length for which segments that length and longer make up >= n% of the total
        bases.  E.g. if n = 50, this function returns the N50.  n must be from 0 to 100.
        """
        total_length = self.get_total_length_no_overlaps()
        target_length = total_length * (n_percent / 100.0)
        sorted_segments = sorted(self.segments.values(),
                                 key=lambda x: x.get_length_no_overlap(self.overlap),
                                 reverse=True)
        length_so_far = 0
        for segment in sorted_segments:
            seg_length = segment.get_length_no_overlap(self.overlap)
            length_so_far += seg_length
            if length_so_far >= target_length:
                return seg_length
        return 0

    def gfa_link_line(self, start, end):
        """
        Returns an entire L line for GFA output, including the newline.
        """
        l_line = 'L\t'
        l_line += str(abs(start)) + '\t'
        l_line += get_sign_string(start) + '\t'
        l_line += str(abs(end)) + '\t'
        l_line += get_sign_string(end) + '\t'
        l_line += str(self.overlap) + 'M\n'
        return l_line

    def get_all_outputs(self, segment):
        """
        Returns a list of segments which lead out from the given segment.
        """
        if segment.number in self.reverse_links:
            return [self.segments[x] for x in self.forward_links[segment.number]]
        else:
            return []

    def get_exclusive_inputs(self, segment_number):
        """
        This function finds all segments which lead into the given segment.  If those segments
        do not lead into any other segments, then this function returns them in a list.  If they
        do lead into other segments, then this function returns None.
        Specifically, this function returns a list of unsigned numbers.
        """
        if segment_number not in self.reverse_links:
            return []
        return [abs(x) for x in self.reverse_links[segment_number] if
                self.lead_exclusively_to(x, segment_number)]

    def get_exclusive_outputs(self, segment_number):
        """
        Does the same thing as get_exclusive_inputs, but in the other direction.
        """
        if segment_number not in self.forward_links:
            return []
        return [abs(x) for x in self.forward_links[segment_number] if
                self.lead_exclusively_from(x, segment_number)]

    def lead_exclusively_to(self, segment_num_1, segment_num_2):
        """
        Returns whether or not the first segment leads to and only to the second segment.
        """
        if segment_num_1 not in self.forward_links:
            return False
        return self.forward_links[segment_num_1] == [segment_num_2]

    def lead_exclusively_from(self, segment_num_1, segment_num_2):
        """
        Does the same thing as lead_exclusively_to, but follows links in the opposite direction.
        """
        if segment_num_1 not in self.reverse_links:
            return False
        return self.reverse_links[segment_num_1] == [segment_num_2]

    def dead_end_change_if_deleted(self, seg_num):
        """
        Returns the change in graph dead end count if this segment was deleted. 0 means no change,
        positive values mean an increase in dead ends, negative values mean a decrease.
        """
        potential_dead_ends = 0
        if seg_num in self.forward_links:
            downstream_segments = self.forward_links[seg_num]
        else:
            downstream_segments = []
        for downstream_segment in downstream_segments:
            if len(self.reverse_links[downstream_segment]) == 1:
                potential_dead_ends += 1

        if seg_num in self.reverse_links:
            upstream_segments = self.reverse_links[seg_num]
        else:
            upstream_segments = []
        for upstream_segment in upstream_segments:
            if len(self.forward_links[upstream_segment]) == 1:
                potential_dead_ends += 1

        return potential_dead_ends - self.dead_end_count(seg_num)

    def dead_end_change_if_path_deleted(self, path_segments):
        """
        Like the above function, but considered the whole path at once. It assumes that the path is
        simple and unbranching (i.e. could be merged into a single segment).
        This function does not check whether the path start and end both connect to the same
        segment. So if they form a hairpin loop, this function will return 0 even though the
        deletion of the path would create a dead end. This behaviour is intentionally left, as it
        helps to clean up such loops from the graph when they have been entirely used in bridges.
        """
        start = path_segments[0]
        end = path_segments[-1]

        potential_dead_ends = 0
        if end in self.forward_links:
            downstream_segments = self.forward_links[end]
        else:
            downstream_segments = []
        for downstream_segment in downstream_segments:
            if len(self.reverse_links[downstream_segment]) == 1:
                potential_dead_ends += 1

        if start in self.reverse_links:
            upstream_segments = self.reverse_links[start]
        else:
            upstream_segments = []
        for upstream_segment in upstream_segments:
            if len(self.forward_links[upstream_segment]) == 1:
                potential_dead_ends += 1

        dead_ends = 0
        if downstream_segments == 0:
            dead_ends += 1
        if upstream_segments == 0:
            dead_ends += 1
        return potential_dead_ends - dead_ends

    def clean(self, read_depth_filter):
        """
        This function does various graph repairs, filters and normalisations to make it a bit
        nicer.
        """
        self.repair_multi_way_junctions()
        self.filter_by_read_depth(read_depth_filter)
        self.filter_homopolymer_loops()
        self.merge_all_possible()
        self.normalise_read_depths()
        self.remove_zero_length_segs(0)
        self.sort_link_order()

    def final_clean(self, verbosity):
        """
        This function cleans up the final assembled graph, in preparation for saving.
        """
        print_section_header('Finalising graph', verbosity, last_newline=(verbosity > 2))
        if self.overlap:
            try:
                self.remove_all_overlaps(verbosity)
                if verbosity > 0:
                    print('\nSuccessfully removed all graph overlaps')
            except CannotTrimOverlaps:
                if verbosity > 0:
                    print('\nUnable to remove graph overlaps')
        self.remove_zero_length_segs(verbosity)
        self.merge_small_segments(verbosity, 5)
        self.normalise_read_depths()
        self.renumber_segments()
        self.sort_link_order()

    def repair_multi_way_junctions(self):
        """
        This function finds and fixes multi-way junctions in the graph, as these can mess up copy
        number determination. It fixes them by creating a new segment with no length (i.e with the
        overlap size) to bridge the connection.
        For example: A->B,C and D->B,C becomes A->E and D->E and E->B and E->C
        """
        while True:
            seg_nums = list(self.segments) + [-x for x in self.segments]
            for seg_num in seg_nums:

                # For the segment, get all of its downstream segments.
                if seg_num not in self.forward_links:
                    continue
                ending_segs = set(self.forward_links[seg_num])
                if len(ending_segs) < 2:
                    continue

                # Now for all of the downstream segments, get their upstream segments.
                starting_segs = set()
                for ending_seg in ending_segs:
                    if ending_seg in self.reverse_links and self.reverse_links[ending_seg]:
                        starting_segs.update(self.reverse_links[ending_seg])
                if len(starting_segs) < 2:
                    continue

                # Now for all of the upstream (starting) segments, get their downstream segments.
                # If this set is the same as the downstream segments of the first segment, then
                # we have ourselves a multi-way junction!
                ending_segs_2 = set()
                for starting_seg in starting_segs:
                    if starting_seg in self.forward_links and self.forward_links[starting_seg]:
                        ending_segs_2.update(self.forward_links[starting_seg])
                if ending_segs_2 != ending_segs:
                    continue

                # Double-check that all of the overlaps agree.
                starting_segs = list(starting_segs)
                ending_segs = list(ending_segs)
                bridge_seq = self.get_seq_from_signed_seg_num(ending_segs[0])[:self.overlap]
                for start_seg in starting_segs:
                    assert bridge_seq == self.get_seq_from_signed_seg_num(start_seg)[-self.overlap:]
                for end_seg in ending_segs:
                    assert bridge_seq == self.get_seq_from_signed_seg_num(end_seg)[:self.overlap]

                # Create a new segment to bridge the starting and ending segments.
                bridge_num = self.get_next_available_seg_number()
                start_seg_depth_sum = sum(self.segments[abs(x)].depth for x in starting_segs)
                end_seg_depth_sum = sum(self.segments[abs(x)].depth for x in ending_segs)
                bridge_depth = (start_seg_depth_sum + end_seg_depth_sum) / 2.0
                bridge_seg = Segment(bridge_num, bridge_depth, bridge_seq, True)
                bridge_seg.build_other_sequence_if_necessary()
                self.segments[bridge_num] = bridge_seg

                # Now rebuild the links around the junction.
                for start_seg in starting_segs:
                    self.forward_links[start_seg] = [bridge_num]
                    self.reverse_links[-start_seg] = [-bridge_num]
                for end_seg in ending_segs:
                    self.reverse_links[end_seg] = [bridge_num]
                    self.forward_links[-end_seg] = [-bridge_num]
                self.forward_links[bridge_num] = ending_segs
                self.reverse_links[bridge_num] = starting_segs
                self.reverse_links[-bridge_num] = [-x for x in ending_segs]
                self.forward_links[-bridge_num] = [-x for x in starting_segs]

                # Finally, we need to check to see if there were any paths through the junction. If
                # so, they need to be adjusted to contain the new segment.
                for name in self.paths:
                    for start_num in starting_segs:
                        for end_num in ending_segs:
                            self.paths[name] = insert_num_in_list(self.paths[name], start_num,
                                                                  end_num, bridge_num)
                            self.paths[name] = insert_num_in_list(self.paths[name], -end_num,
                                                                  -start_num, -bridge_num)
                break
            else:
                break

    def get_next_available_seg_number(self):
        """
        This function finds the largest used segment number and returns the next
        """
        current_largest = max(self.segments)
        return current_largest + 1

    def get_depth_string(self, segment):
        """
        Given a particular segment, this function returns a string with the segment's copy depths
        (if it has any).
        """
        if segment.number not in self.copy_depths:
            return ''
        return ', '.join(['%.3f' % x for x in self.copy_depths[segment.number]])

    def get_copy_number_colour(self, segment):
        """
        Given a particular segment, this function returns a colour string based on the copy number.
        """
        if segment.number not in self.copy_depths:
            return 'black'
        copy_number = len(self.copy_depths[segment.number])
        if copy_number == 1:
            return 'forestgreen'
        if copy_number == 2:
            return 'gold'
        if copy_number == 3:
            return 'darkorange'
        else:  # 4+
            return 'red'

    def determine_copy_depth(self, verbosity):
        """
        Assigns a copy depth to each segment in the graph.
        """
        # Reset any existing copy depths.
        self.copy_depths = {}

        # TO DO: These should be parameters, after I have them sorted out.
        initial_tolerance = 0.1
        propagation_tolerance = 0.5
        min_half_median_for_diploid = 0.1
        min_single_copy_length = 1000

        # Determine the single-copy read depth for the graph. In haploid and some diploid cases,
        # this will be the median depth. But in some diploid cases, the single-copy depth may at
        # about half the median (because the median depth holds the sequences shared between sister
        # chromosomes). To catch these cases, we look to see whether the graph peaks more strongly
        # at half the median or double the median. In the former case, we move the single-copy
        # depth down to half the median.
        median_depth = self.get_median_read_depth()
        if verbosity > 1:
            print('Median graph depth:', float_to_str(median_depth, 3))
        bases_near_half_median = self.get_base_count_in_depth_range(median_depth * 0.4,
                                                                    median_depth * 0.6)
        bases_near_double_median = self.get_base_count_in_depth_range(median_depth * 1.6,
                                                                      median_depth * 2.4)
        total_graph_bases = self.get_total_length()
        half_median_frac = bases_near_half_median / total_graph_bases
        double_median_frac = bases_near_double_median / total_graph_bases
        if half_median_frac > double_median_frac and \
                half_median_frac >= min_half_median_for_diploid:
            single_copy_depth = median_depth / 2.0
        else:
            single_copy_depth = median_depth
        if verbosity > 1:
            print('Single-copy depth:', float_to_str(median_depth, 3))

        # Assign single-copy status to segments within the tolerance of the single-copy depth.
        max_depth = single_copy_depth + initial_tolerance
        initial_single_copy_segments = []
        for seg_num, segment in self.segments.items():
            if segment.depth <= max_depth and \
                    self.at_most_one_link_per_end(segment):
                self.copy_depths[segment.number] = [segment.depth]
                initial_single_copy_segments.append(seg_num)
        if verbosity > 1:
            if initial_single_copy_segments:
                print()
                print('Initial single copy segments:',
                      ', '.join([str(x) for x in initial_single_copy_segments]))
            else:
                print('Initial single copy segments: none')
            print_section_header('Propagating copy numbers', verbosity)

        # Propagate copy depth as possible using those initial assignments.
        self.determine_copy_depth_part_2(propagation_tolerance, verbosity)

        # Assign single-copy to the largest available segment, propagate and repeat.
        while True:
            assignments = self.assign_single_copy_depth(verbosity, min_single_copy_length)
            self.determine_copy_depth_part_2(propagation_tolerance, verbosity)
            if not assignments:
                break

        # Now propagate with no tolerance threshold to complete the remaining segments.
        self.determine_copy_depth_part_2(1.0, verbosity)

    def determine_copy_depth_part_2(self, tolerance, verbosity):
        """
        Propagates copy depth repeatedly until assignments stop.
        """
        while self.merge_copy_depths(tolerance, verbosity):
            pass
        if self.redistribute_copy_depths(tolerance, verbosity):
            self.determine_copy_depth_part_2(tolerance, verbosity)

    def assign_single_copy_depth(self, verbosity, min_single_copy_length):
        """
        This function assigns a single copy to the longest available segment.
        """
        segments = sorted(self.get_segments_without_copies(),
                          key=lambda x: x.get_length(), reverse=True)
        for segment in segments:
            if segment.get_length() < min_single_copy_length:
                continue
            if self.exactly_one_link_per_end(segment):
                self.copy_depths[segment.number] = [segment.depth]
                if verbosity > 1:
                    print('New single copy:', segment.number,
                          '(' + float_to_str(segment.depth, 2) + 'x)')
                return 1
        return 0

    def merge_copy_depths(self, error_margin, verbosity):
        """
        This function looks for segments where they have input on one end where:
          1) All input segments have copy depth assigned.
          2) All input segments exclusively input to this segment.
        All such cases are evaluated, and the segment with the lowest error (if that error is below
        the allowed error margin) is assigned copy depths, scaling the inputs so their sum
        exactly matches the segment's depth.
        """
        segments = self.get_segments_without_copies()
        if not segments:
            return 0

        best_segment_num = None
        best_source_nums = None
        best_new_depths = []
        lowest_error = float('inf')

        for segment in segments:
            num = segment.number
            exclusive_inputs = self.get_exclusive_inputs(num)
            exclusive_outputs = self.get_exclusive_outputs(num)
            in_depth_possible = exclusive_inputs and self.all_have_copy_depths(exclusive_inputs)
            out_depth_possible = exclusive_outputs and self.all_have_copy_depths(exclusive_outputs)
            if in_depth_possible:
                depths, error = self.scale_copy_depths_from_source_segments(num, exclusive_inputs)
                if error < lowest_error:
                    lowest_error = error
                    best_segment_num = num
                    best_source_nums = exclusive_inputs
                    best_new_depths = depths
            if out_depth_possible:
                depths, error = self.scale_copy_depths_from_source_segments(num, exclusive_outputs)
                if error < lowest_error:
                    lowest_error = error
                    best_segment_num = num
                    best_source_nums = exclusive_outputs
                    best_new_depths = depths
        if best_segment_num and lowest_error < error_margin:
            self.copy_depths[best_segment_num] = best_new_depths
            if verbosity > 1:
                print('Merged copies:  ',
                      ' + '.join([str(x) + ' (' + float_to_str(self.segments[x].depth, 2) + 'x)'
                                  for x in best_source_nums]), '->',
                      best_segment_num,
                      '(' + float_to_str(self.segments[best_segment_num].depth, 2) + 'x)')
            return 1
        else:
            return 0

    def redistribute_copy_depths(self, error_margin, verbosity):
        """
        This function deals with the easier case of copy depth redistribution: where one segments
        with copy depth leads exclusively to multiple segments without copy depth.
        We will then try to redistribute the source segment's copy depths among the destination
        segments.  If it can be done within the allowed error margin, the destination segments will
        get their copy depths.
        """
        segments = self.get_segments_with_two_or_more_copies()
        if not segments:
            return 0
        for segment in segments:
            num = segment.number
            connections = self.get_exclusive_inputs(num)
            if not connections or self.all_have_copy_depths(connections):
                connections = self.get_exclusive_outputs(num)
            if not connections or self.all_have_copy_depths(connections):
                continue

            # If we got here, then we can try to redistribute the segment's copy depths to its
            # connections which are lacking copy depth.
            copy_depths = self.copy_depths[num]
            bins = [[]] * len(connections)
            targets = [None if x not in self.copy_depths else len(self.copy_depths[x])
                       for x in connections]

            # For cases where there are many copy depths being distributed to many segments, there
            # will be too many combinations, so we don't bother trying.
            arrangement_count = len(bins) ** len(copy_depths)
            if arrangement_count > 10000:  # TO DO: choose this cutoff more carefully?
                continue

            arrangements = shuffle_into_bins(copy_depths, bins, targets)
            if not arrangements:
                continue

            lowest_error = float('inf')
            for i, arrangement in enumerate(arrangements):
                error = self.get_error_for_multiple_segments_and_depths(connections, arrangement)
                if i == 0 or error < lowest_error:
                    lowest_error = error
                    best_arrangement = arrangement
            if lowest_error < error_margin:
                if self.assign_copy_depths_where_needed(connections, best_arrangement,
                                                        error_margin):
                    if verbosity > 1:
                        print('Split copies:   ', num,
                              '(' + float_to_str(self.segments[num].depth, 2) + 'x) ->',
                              ' + '.join([str(x) + ' (' +
                                          float_to_str(self.segments[x].depth, 2) + 'x)'
                                          for x in connections]))
                    return 1
        return 0

    def at_most_one_link_per_end(self, segment):
        """
        Returns True if the given segment has no more than one link on either end.
        """
        num = segment.number
        if num in self.forward_links and len(self.forward_links[num]) > 1:
            return False
        if num in self.reverse_links and len(self.reverse_links[num]) > 1:
            return False
        return True

    def exactly_one_link_per_end(self, segment):
        """
        Returns True if the given segment has exactly one link on either end.
        """
        num = segment.number
        if num in self.forward_links and len(self.forward_links[num]) != 1:
            return False
        if num in self.reverse_links and len(self.reverse_links[num]) != 1:
            return False
        return True

    def all_have_copy_depths(self, segment_numbers):
        """
        Takes a list of segment numbers and returns whether every segment in the list has copy
        depths assigned.
        """
        for num in segment_numbers:
            if num not in self.copy_depths:
                return False
        return True

    def scale_copy_depths_from_source_segments(self, segment_number, source_segment_numbers):
        """
        Using a list of segments which are the source of copy depth, this function scales them so
        that their sum matches the depth of the given segment.
        It returns:
          1) a list of depth numbers
          2) the error (i.e. the degree of scaling which had to occur)
        It assumes that all of the source segments definitely have copy depths.
        """
        source_depths = []
        for num in source_segment_numbers:
            source_depths += self.copy_depths[num]
        target_depth = self.segments[segment_number].depth
        return self.scale_copy_depths(target_depth, source_depths)

    @staticmethod
    def scale_copy_depths(target_depth, source_depths):
        """
        This function takes the source depths and scales them so their sum matches the target
        depth.  It returns the scaled depths and the error.
        """
        source_depth_sum = sum(source_depths)
        scaling_factor = target_depth / source_depth_sum
        scaled_depths = sorted([scaling_factor * x for x in source_depths], reverse=True)
        error = get_error(source_depth_sum, target_depth)
        return scaled_depths, error

    def get_segments_without_copies(self):
        """
        Returns a list of the graph segments lacking copy depth information.
        """
        return [x for x in self.segments.values() if x.number not in self.copy_depths]

    def get_segments_with_two_or_more_copies(self):
        return [x for x in self.segments.values() if
                x.number in self.copy_depths and len(self.copy_depths[x.number]) > 1]

    def get_error_for_multiple_segments_and_depths(self, segment_numbers, copy_depths):
        """
        For the given segments, this function assesses how well the given copy depths match up.
        The maximum error for any segment is what's returned at the end.
        """
        max_error = 0.0
        for i, num in enumerate(segment_numbers):
            segment_depth = self.segments[num].depth
            depth_sum = sum(copy_depths[i])
            max_error = max(max_error, get_error(depth_sum, segment_depth))
        return max_error

    def assign_copy_depths_where_needed(self, segment_numbers, new_depths, error_margin):
        """
        For the given segments, this function assigns the corresponding copy depths, scaled to fit
        the segment.  If a segment already has copy depths, it is skipped (i.e. this function only
        write new copy depths, doesn't overwrite existing ones).
        It will only create copy depths if doing so is within the allowed error margin.
        """
        success = False
        for i, num in enumerate(segment_numbers):
            if num not in self.copy_depths:
                new_copy_depths, error = self.scale_copy_depths(self.segments[num].depth,
                                                                new_depths[i])
                if error <= error_margin:
                    self.copy_depths[num] = new_copy_depths
                    success = True
        return success

    def remove_segment_depth(self, seg_num, depth_to_remove):
        """
        This function removes the depth from the specified segment. If the segment has copy depths,
        it will also removes the copy depth closest to the given depth.
        """
        seg_num = abs(seg_num)
        if seg_num not in self.segments:
            return
        seg = self.segments[seg_num]
        seg.depth = max(0.0, seg.depth - depth_to_remove)
        if seg_num not in self.copy_depths or not self.copy_depths[seg_num]:
            return
        closest_depth = min(self.copy_depths[seg_num], key=lambda x: abs(x - depth_to_remove))
        del self.copy_depths[seg_num][self.copy_depths[seg_num].index(closest_depth)]

    def get_base_count_in_depth_range(self, min_depth, max_depth):
        """
        Returns the total number of bases in the graph in the given depth range.
        """
        total_bases = 0
        for segment in self.segments.values():
            if min_depth <= segment.depth <= max_depth:
                total_bases += segment.get_length()
        return total_bases

    def get_single_copy_segments(self):
        """
        Returns a list of the graph segments with a copy number of 1.
        """
        single_copy_segments = []
        for num, segment in self.segments.items():
            if num in self.copy_depths and len(self.copy_depths[num]) == 1:
                single_copy_segments.append(segment)
        return single_copy_segments

    def is_path_valid(self, path_segments):
        """
        Returns whether or not a path exists in the graph, i.e. an edge connects each segment in
        the path to the next segment in the path.
        """
        for i, seg_num in enumerate(path_segments):
            if abs(seg_num) not in self.segments:
                return False
            if i > 0:
                prev_seg = path_segments[i - 1]
                if seg_num not in self.get_downstream_seg_nums(prev_seg):
                    return False
        return True

    def get_path_sequence(self, path_segments):
        """
        Gets a linear (i.e. not circular) path sequence from the graph.
        """
        path_sequence = ''
        prev_segment_number = None
        for i, seg_num in enumerate(path_segments):
            segment = self.segments[abs(seg_num)]
            if seg_num > 0:
                seg_sequence = segment.forward_sequence
            else:
                seg_sequence = segment.reverse_sequence
            if i == 0:
                path_sequence = seg_sequence
            else:
                assert seg_num in self.forward_links[prev_segment_number]
                if self.overlap > 0:
                    assert path_sequence[-self.overlap:] == seg_sequence[:self.overlap]
                path_sequence += seg_sequence[self.overlap:]
            prev_segment_number = seg_num
        return path_sequence

    def apply_bridges(self, bridges, verbosity, min_bridge_qual, single_copy_segments):
        """
        Uses the supplied bridges to simplify the graph.
        """
        # Each segment can have only one bridge per side, so we will track which segments have had
        # a bridge applied off one side or the other.
        right_bridged = set()
        left_bridged = set()

        # Sort bridges first by type: LongReadBridge, SpadesContigBridge and then
        # LoopUnrollingBridge. Then sort by  quality so within each type we apply the best bridges
        # first.
        sorted_bridges = sorted(bridges, key=lambda x: (x.get_type_score(), x.quality),
                                reverse=True)
        bridge_segs = []
        seg_nums_used_in_bridges = []
        single_copy_nums = set(x.number for x in single_copy_segments)

        for bridge in sorted_bridges:

            # If a bridge has multiple equally-good graph paths, then we can choose which one to
            # use based on the availability of the path.
            if isinstance(bridge, LongReadBridge) and len(bridge.all_best_paths) > 1:
                best_path = bridge.all_best_paths[0][0]
                best_availability = self.get_path_availability(best_path)
                for i in range(1, len(bridge.all_best_paths)):
                    potential_path = bridge.all_best_paths[i][0]
                    potential_path_availability = self.get_path_availability(potential_path)
                    if potential_path_availability > best_availability:
                        best_availability = potential_path_availability
                        best_path = potential_path
                bridge.graph_path = best_path
                if self.is_path_valid(best_path):
                    bridge.bridge_sequence = self.get_path_sequence(best_path)

            # Get the pieces of the bridge which can be applied.
            pieces = get_applicable_bridge_pieces(bridge, single_copy_nums, right_bridged,
                                                  left_bridged, seg_nums_used_in_bridges)

            # If no parts of the bridge can be applied, then there's nothing to bridge.
            if not pieces:
                if verbosity > 1:
                    print('Unused', bridge)
                continue

            # If the bridge's quality is too low, we don't use it.
            if bridge.quality < min_bridge_qual:
                if verbosity > 1:
                    print('Rejected', bridge)
                continue

            # If the bridge can be applied as a whole, go ahead and apply it.
            if len(pieces) == 1 and \
                    pieces[0] == [bridge.start_segment] + bridge.graph_path + [bridge.end_segment]:
                bridge = self.apply_entire_bridge(bridge, verbosity, right_bridged, left_bridged,
                                                  seg_nums_used_in_bridges, single_copy_nums)
                bridge_segs.append(bridge)

            else:  # Either multiple pieces or a single piece which isn't the whole bridge.
                bridges = self.apply_bridge_in_pieces(bridge, pieces, verbosity, right_bridged,
                                                      left_bridged, seg_nums_used_in_bridges,
                                                      single_copy_nums)
                bridge_segs += bridges

            # Remove duplicates, while preserving order.
            seen = set()
            seg_nums_used_in_bridges = [x for x in seg_nums_used_in_bridges
                                        if not (x in seen or seen.add(x))]
        return seg_nums_used_in_bridges

    def clean_up_after_bridging(self, single_copy_segments, seg_nums_used_in_bridges,
                                min_component_size, min_dead_end_size, verbosity):
        """
        Cleans up unnecessary segments to produce a clean graph. Used after bridge application.
        """
        if verbosity > 1:
            print_section_header('Cleaning up leftover segments', verbosity, last_newline=False)

        # Graph components which contain no single-copy segments are deleted.
        single_copy_seg_nums = set(x.number for x in single_copy_segments)
        segment_nums_to_remove = []
        connected_components = self.get_connected_components()
        for component_nums in connected_components:
            for seg_num in component_nums:
                if abs(seg_num) in single_copy_seg_nums:
                    break
            else:
                segment_nums_to_remove += component_nums
        if verbosity > 1 and segment_nums_to_remove:
            print('\nRemoved components with no single-copy segments:',
                  ', '.join(str(x) for x in segment_nums_to_remove))
        self.remove_segments(segment_nums_to_remove)

        # Graph components which have been entirely used in bridges are deleted.
        segment_nums_to_remove = []
        connected_components = self.get_connected_components()
        for component_nums in connected_components:
            for seg_num in component_nums:
                if abs(seg_num) not in seg_nums_used_in_bridges:
                    break
            else:
                segment_nums_to_remove += component_nums
        if verbosity > 1 and segment_nums_to_remove:
            print('\nRemoved components used in bridges:',
                  ', '.join(str(x) for x in segment_nums_to_remove))
        self.remove_segments(segment_nums_to_remove)

        # Any multi-copy segments which cannot possibly connect two single-copy segments are
        # deleted.
        segment_nums_to_remove = []
        for seg_num in self.segments:
            if seg_num in single_copy_seg_nums:
                continue
            if not (self.search(seg_num, single_copy_seg_nums) and
                    self.search(-seg_num, single_copy_seg_nums)):
                segment_nums_to_remove.append(seg_num)
        if verbosity > 1 and segment_nums_to_remove:
            print('\nRemoved unbridging segments:',
                  ', '.join(str(x) for x in segment_nums_to_remove))
        self.remove_segments(segment_nums_to_remove)

        # Remove segments used in bridges, if doing do would not break up the graph.
        removed_segments = []
        while True:
            for seg_num in seg_nums_used_in_bridges:
                if seg_num not in self.segments:
                    continue

                # If the segment already has a dead end, just go ahead and delete it.
                if self.dead_end_count(seg_num) > 0:
                    self.remove_segments([seg_num])
                    removed_segments.append(seg_num)
                    break

                # If deleting the segment would not lead to an increase in dead ends, go ahead and
                # delete it.
                elif self.dead_end_change_if_deleted(seg_num) <= 0:
                    self.remove_segments([seg_num])
                    removed_segments.append(seg_num)
                    break

                # It's possible that a group of segments have all been used and they form a simple
                # path. Deleting the whole path would not create a dead end (but deleting any one
                # segment would), so we look for these and delete them all together.
                else:
                    path = self.get_simple_path(seg_num)
                    unsigned_path = [abs(x) for x in path]
                    if len(path) > 1 and \
                            all(x in seg_nums_used_in_bridges for x in unsigned_path) and \
                            self.dead_end_change_if_path_deleted(path) <= 0:
                        self.remove_segments(unsigned_path)
                        removed_segments += unsigned_path
                        break
            else:
                break

        # Clean up connected components which have been entirely used in bridges.
        segment_nums_to_remove = []
        connected_components = self.get_connected_components()
        for component in connected_components:
            component_seg_nums = [self.segments[x].number for x in component]
            for component_seg_num in component_seg_nums:
                if component_seg_num not in seg_nums_used_in_bridges:
                    break
            else:
                segment_nums_to_remove += component_seg_nums
        self.remove_segments(segment_nums_to_remove)
        removed_segments += segment_nums_to_remove
        if verbosity > 1 and segment_nums_to_remove:
            print('\nRemoved segments used in bridges', ', '.join(str(x) for x in removed_segments))

        self.remove_small_components(min_component_size, verbosity)
        self.remove_small_dead_ends(min_dead_end_size, verbosity)

    def apply_entire_bridge(self, bridge, verbosity, right_bridged, left_bridged,
                            seg_nums_used_in_bridges, single_copy_nums):
        """
        Applies a whole bridge, start to end.
        """
        if verbosity > 1:
            print('Applying', bridge)
        new_seg = self.apply_bridge(bridge, bridge.start_segment, bridge.end_segment,
                                    bridge.bridge_sequence, bridge.graph_path, verbosity)
        for seg_num in bridge.graph_path:
            single_copy_nums.discard(abs(seg_num))
        add_to_bridged_sets(bridge.start_segment, bridge.end_segment, right_bridged, left_bridged)
        seg_nums_used_in_bridges.extend([abs(x) for x in bridge.graph_path])
        return new_seg

    def apply_bridge_in_pieces(self, bridge, pieces, verbosity, right_bridged, left_bridged,
                               seg_nums_used_in_bridges, single_copy_nums):
        """
        Applies a bridge in chunks, as appropriate. Used for bridges which can't be applied in
        their entirety.
        """
        new_segs = []
        if verbosity > 1:
            if len(pieces) == 1:
                print('Applying 1 piece of', bridge)
            else:
                print('Applying ' + str(len(pieces)) + ' pieces of', bridge)

        for piece in pieces:
            piece_start = piece[0]
            piece_middle = piece[1:-1]
            piece_end = piece[-1]
            piece_seq = self.get_path_sequence(piece_middle)
            if verbosity > 1:
                print('        ' + get_bridge_str(piece_start, piece_middle, piece_end))
            new_seg = self.apply_bridge(bridge, piece_start, piece_end, piece_seq, piece_middle,
                                        verbosity)
            new_segs.append(new_seg)
            for seg_num in piece_middle:
                single_copy_nums.discard(abs(seg_num))
            add_to_bridged_sets(piece_start, piece_end, right_bridged, left_bridged)
            seg_nums_used_in_bridges.extend([abs(x) for x in piece_middle])
        return new_segs

    def apply_bridge(self, bridge, start, end, sequence, graph_path, verbosity):
        """
        Applies one bridge to the graph, from the start segment to the end and with the given
        sequence. This may be the entire bridge, or possibly just a piece of the bridge.
        """
        # Remove all existing links for the segments being bridged.
        if start in self.forward_links:
            for link in self.forward_links[start]:
                self.remove_link(start, link)
        if end in self.reverse_links:
            for link in self.reverse_links[end]:
                self.remove_link(link, end)

        # Create a new bridge segment.
        new_seg_num = self.get_next_available_seg_number()
        new_seg = Segment(new_seg_num, bridge.depth, sequence, True, bridge, graph_path)
        new_seg.build_other_sequence_if_necessary()
        self.segments[new_seg_num] = new_seg

        # Link the bridge segment in to the start/end segments.
        self.add_link(start, new_seg_num)
        self.add_link(new_seg_num, end)

        # Subtract the bridge depth from the segments in the bridge.
        for seg in graph_path:
            self.remove_segment_depth(seg, bridge.depth)

        return new_seg

    def find_all_simple_loops(self):
        """
        This function finds all cases of a simple loop in the graph: A->B->C->B->D.
        It returns them as a list of 4-tuples of segment numbers in this order:
        (start, end, middle, repeat).
        """
        simple_loops = []

        # We'll search specifically for the middle segments as they should be easy to spot.
        for middle in self.segments:

            # A middle segment will always have exactly one connection on each end which connect
            # to the same segment (the repeat segment).
            if middle not in self.forward_links or middle not in self.reverse_links:
                continue
            if len(self.forward_links[middle]) != 1 or len(self.reverse_links[middle]) != 1:
                continue
            if self.forward_links[middle][0] != self.reverse_links[middle][0]:
                continue
            repeat = self.forward_links[middle][0]

            # The repeat segment should have exactly two connections on each end. If less, then we
            # have a simple path which can be merged. If more, it's a more complex loop.
            if len(self.forward_links[repeat]) != 2 or len(self.reverse_links[repeat]) != 2:
                continue

            # Find the start and end segment numbers. It's okay if the start and the end are the
            # same, but we exclude any other screwy cases where the start or end is the middle or
            # repeat segment.
            start = self.reverse_links[repeat][0]
            if abs(start) == abs(middle):
                start = self.reverse_links[repeat][1]
            if abs(start) == abs(middle) or abs(start) == abs(repeat):
                continue

            end = self.forward_links[repeat][0]
            if abs(end) == abs(middle):
                end = self.forward_links[repeat][1]
            if abs(end) == abs(middle) or abs(end) == abs(repeat):
                continue

            simple_loops.append((start, end, middle, repeat))
        return simple_loops

    def all_paths(self, start, end, min_length, target_length, max_length, max_path_count):
        """
        Returns a list of all paths which connect the starting segment to the ending segment and
        are within the length bounds. The start and end segments are not themselves included in the
        paths. Returns an empty list if no paths exist.
        Loops in the graph (especially loops of short segments which don't add much to the path
        length) can result in very large numbers of potential paths in complex areas. To somewhat
        manage this, we exclude paths which include too many copies of a segment. 'Too many copies'
        is defined as double the copy depth count or the depth to start/end count.
        """
        if start not in self.forward_links:
            return []

        start_seg = self.segments[abs(start)]
        end_seg = self.segments[abs(end)]
        start_end_depth = weighted_average(start_seg.depth, end_seg.depth,
                                           start_seg.get_length_no_overlap(self.overlap),
                                           end_seg.get_length_no_overlap(self.overlap))
        max_working_paths = 1000  # TO DO: make this a parameter?
        working_paths = [[x] for x in self.forward_links[start]]
        final_paths = []
        while working_paths:
            new_working_paths = []
            for working_path in working_paths:
                last_seg = working_path[-1]
                if last_seg == end:
                    potential_result = working_path[:-1]
                    if self.get_path_length(potential_result) >= min_length:
                        final_paths.append(potential_result)
                        if len(final_paths) > max_path_count:
                            raise TooManyPaths
                elif self.get_path_length(working_path) <= max_length and \
                        last_seg in self.forward_links:
                    for next_seg in self.forward_links[last_seg]:
                        max_allowed_count = self.max_path_segment_count(next_seg, start_end_depth)
                        count_so_far = working_path.count(next_seg) + working_path.count(-next_seg)
                        if count_so_far < max_allowed_count:
                            new_working_paths.append(working_path + [next_seg])

            # If the number of working paths is too high, we give up.
            if len(working_paths) > max_working_paths:
                raise TooManyPaths
            working_paths = new_working_paths

        # Sort by length discrepancy from the target so the closest length matches come first.
        final_paths = sorted(final_paths,
                             key=lambda x: abs(target_length - self.get_path_length(x)))

        return final_paths

    def progressive_path_find(self, start, end, min_length, target_length, max_length, sequence,
                              scoring_scheme):
        """
        Finds paths from the start to the end that best match the sequence. It works outward and
        culls paths when they grow too numerous. This function is called when all_paths fails due to
        too many paths. The start and end segments are not themselves included in the paths.
        """
        if start not in self.forward_links:
            return []

        start_seg = self.segments[abs(start)]
        end_seg = self.segments[abs(end)]
        start_end_depth = weighted_average(start_seg.depth, end_seg.depth,
                                           start_seg.get_length_no_overlap(self.overlap),
                                           end_seg.get_length_no_overlap(self.overlap))

        # When the number of paths exceeds max_working_paths, they are all evaluated and the best
        # keep_count paths will be kept.
        max_working_paths = 100
        keep_count = 10

        final_paths = []
        best_final_path_score = 0.0
        working_paths = [[x] for x in self.forward_links[start]]

        # print('\n\n\n\n\n\n\n')
        # print(start, 'to', end, 'progressive path finding:')
        # print('  target length:', target_length)

        while working_paths:

            # Find the length of the shortest working path.
            shortest_len = min(self.get_path_length(x) for x in working_paths)

            # Extend the shortest working path(s) by adding downstream segments. Check to see if
            # this finishes the path or makes it excessively long.
            new_working_paths = []
            for path in working_paths:
                path_len = self.get_path_length(path)
                if path_len == shortest_len:
                    if path[-1] in self.forward_links:
                        downstream_segments = self.forward_links[path[-1]]
                        for next_seg in downstream_segments:
                            if next_seg == end:
                                if min_length <= path_len <= max_length:
                                    final_paths.append(path)
                                    path_seq = self.get_path_sequence(path)
                                    alignment_result = fully_global_alignment(path_seq, sequence,
                                                                              scoring_scheme, True,
                                                                              1000)
                                    if alignment_result:
                                        seqan_parts = alignment_result.split(',', 9)
                                        scaled_score = float(seqan_parts[7])
                                        best_final_path_score = max(best_final_path_score,
                                                                    scaled_score)
                                        # print('FINAL PATH:')
                                        # print('  ' + ','.join(str(x) for x in path) + ' (' +
                                        #       int_to_str(self.get_path_length(
                                        #           path)) + ' bp, score = ' +
                                        #       float_to_str(scaled_score, 2) + ')')
                            else:
                                max_allowed_count = self.max_path_segment_count(next_seg,
                                                                                start_end_depth)
                                count_so_far = path.count(next_seg) + path.count(-next_seg)
                                if count_so_far < max_allowed_count:
                                    extended_path = path + [next_seg]
                                    if self.get_path_length(extended_path) <= max_length:
                                        new_working_paths.append(extended_path)
                else:
                    new_working_paths.append(path)

            # print('\n  Working paths:', len(new_working_paths))
            # for path in new_working_paths:
            #     print('  ' + ','.join(str(x) for x in path) + ' (' +
            #           int_to_str(self.get_path_length(path)) + ' bp)')

            # If we've acquired too many working paths, cull out the worst ones.
            if len(new_working_paths) > max_working_paths:
                # print('\n  CULL\n')

                # It's possible at this point that all of the working paths share quite a bit in
                # common at their start. We can therefore find the common starting sequence and
                # align to that once, and then only do separate alignments for the remainder of
                # the paths. This can save a lot of time!
                common_start = []
                smallest_seg_count = min(len(x) for x in new_working_paths)
                for i in range(smallest_seg_count):
                    potential_common_seg = new_working_paths[0][i]
                    for path in new_working_paths:
                        if path[i] != potential_common_seg:
                            break
                    else:
                        common_start.append(potential_common_seg)
                        continue
                    break
                # print('Common path: ', common_start)
                common_path_seq = self.get_path_sequence(common_start)[:-100]
                path_seq_align_start = len(common_path_seq)
                consensus_seq_align_start = 0
                if common_path_seq:
                    alignment_result = path_alignment(common_path_seq, sequence, scoring_scheme,
                                                      True, 1000)
                    # print(alignment_result)
                    consensus_seq_align_start = int(alignment_result.split(',')[5])
                    if consensus_seq_align_start >= len(sequence):
                        consensus_seq_align_start = 0
                    # print('Alignment start pos in paths:    ', path_seq_align_start)
                    # print('Alignment start pos in consensus:', consensus_seq_align_start)

                scored_paths = []
                shortest_len = min(self.get_path_length(x) for x in new_working_paths)
                for path in new_working_paths:
                    path_seq = self.get_path_sequence(path)[path_seq_align_start:shortest_len]
                    consensus_seq = sequence[consensus_seq_align_start:]
                    alignment_result = path_alignment(path_seq, consensus_seq, scoring_scheme, True,
                                                      500)
                    if not alignment_result:
                        continue
                    seqan_parts = alignment_result.split(',', 9)
                    scaled_score = float(seqan_parts[7])
                    # print('  ' + ','.join(str(x) for x in path) + ' (score = ' +
                    #       float_to_str(scaled_score, 2) + ')')
                    scored_paths.append((path, scaled_score))
                scored_paths = sorted(scored_paths, key=lambda x: x[1], reverse=True)

                # Keep the keep count number of paths, but go over the keep count if there's a tie.
                # E.g. if the keep count is 12 but the top 16 paths all have the same score, keep
                # all 16.
                surviving_paths = []
                if scored_paths:
                    best_score = scored_paths[0][1]
                    for i, scored_path in enumerate(scored_paths):
                        if i >= keep_count and scored_paths[i][1] < best_score:
                            break
                        if best_final_path_score == 0.0 or \
                                scored_paths[i][1] >= 0.75 * best_final_path_score:
                            surviving_paths.append(scored_paths[i])

                # If any of the surviving paths ends in the same segment but have different
                # scores, only keep one with the top score. This is because the segments with a
                # lower score will have the same future paths, and so will always be lower. And
                # segments with the same score will have the same future paths, and so will always
                # be the same.
                surviving_paths_by_terminal_seg = {}
                for surviving_path in surviving_paths:
                    terminal_seg = surviving_path[0][-1]
                    if terminal_seg not in surviving_paths_by_terminal_seg or \
                            surviving_path[1] > surviving_paths_by_terminal_seg[terminal_seg][1]:
                        surviving_paths_by_terminal_seg[terminal_seg] = surviving_path
                working_paths = list(x[0] for x in surviving_paths_by_terminal_seg.values())

                # print('\n  SURVIVING PATHS:', len(working_paths))
                # for i, path in enumerate(working_paths):
                #     print('  ' + ','.join(str(x) for x in path) + ' (' +
                #           int_to_str(self.get_path_length(path)) + ' bp, score = ' +
                #           int_to_str(scored_paths[i][1]) + ')')
            else:
                working_paths = new_working_paths

        # print('\n  FINAL PATHS:')
        # for path in final_paths:
        #     print('  ' + ','.join(str(x) for x in path))
        # print()

        # Sort by length discrepancy from the target so the closest length matches come first.
        final_paths = sorted(final_paths,
                             key=lambda x: abs(target_length - self.get_path_length(x)))

        return final_paths

    def max_path_segment_count(self, seg_num, start_end_depth):
        """
        This function returns the maximum allowed number of times a segment can be in a bridge
        path. It uses both the segment's copy depth (if it has one) and the relative depth of the
        segment as compared to the start/end of the bridge.
        """
        if abs(seg_num) in self.copy_depths:
            count_by_copies = len(self.copy_depths[abs(seg_num)])
        else:
            count_by_copies = 1
        depth = self.segments[abs(seg_num)].depth
        count_by_depth = max(1, int(round(depth / start_end_depth)))
        return 2 * max(count_by_copies, count_by_depth)

    def get_path_length(self, path):
        """
        Returns the length of the given path.
        """
        try:
            path_length = 0
            for seg in path:
                path_length += self.segments[abs(seg)].get_length()
            overlap_count = len(path) - 1
            path_length -= overlap_count * self.overlap
            return path_length
        except KeyError:
            return 0

    def renumber_segments(self):
        """
        This function gives the longest segment the number 1, the second-longest the number 2, etc.
        """
        old_nums = [x.number for x in sorted(self.segments.values(), reverse=True,
                                             key=lambda x: x.get_length())]
        new_nums = list(range(1, len(old_nums) + 1))
        old_nums += [-x for x in old_nums]
        new_nums += [-x for x in new_nums]
        changes = dict(zip(old_nums, new_nums))

        new_segments = {}
        for seg_num, seg in self.segments.items():
            new_num = changes[seg_num]
            seg.number = new_num
            new_segments[new_num] = seg
        self.segments = new_segments

        new_forward_links = {}
        for seg_num, link_nums in self.forward_links.items():
            new_forward_links[changes[seg_num]] = [changes[x] for x in link_nums]
        self.forward_links = new_forward_links

        new_reverse_links = {}
        for seg_num, link_nums in self.reverse_links.items():
            new_reverse_links[changes[seg_num]] = [changes[x] for x in link_nums]
        self.reverse_links = new_reverse_links

        self.copy_depths = {changes[x]: y for x, y in self.copy_depths.items()}

        new_paths = {}
        for name, path_nums in self.paths.items():
            new_paths[name] = [changes[x] for x in path_nums]
        self.paths = new_paths

    def get_summary(self, verbosity, file=None, score=None, show_components=False):
        """
        Returns a nice table describing the graph.
        """
        total_length = self.get_total_length()
        max_v = max(total_length, 1000000)
        n50, shortest, lower_quartile, median, upper_quartile, longest = self.get_contig_stats()
        summary = ''
        if file:
            summary += file + '\n'
        summary += 'segments:              ' + int_to_str(len(self.segments), max_v) + '\n'
        summary += 'links:                 ' + int_to_str(self.get_total_link_count(), max_v) + '\n'
        summary += 'total length (bp):     ' + int_to_str(total_length, max_v) + '\n'
        summary += 'N50:                   ' + int_to_str(n50, max_v) + '\n'
        if verbosity > 2:
            summary += 'shortest segment (bp): ' + int_to_str(shortest, max_v) + '\n'
            summary += 'lower quartile (bp):   ' + int_to_str(lower_quartile, max_v) + '\n'
            summary += 'median segment (bp):   ' + int_to_str(median, max_v) + '\n'
            summary += 'upper quartile (bp):   ' + int_to_str(upper_quartile, max_v) + '\n'
        summary += 'longest segment (bp):  ' + int_to_str(longest, max_v) + '\n'
        summary += 'dead ends:             ' + int_to_str(self.total_dead_end_count(), max_v) + '\n'
        if show_components:
            components = self.get_connected_components()
            for i, component in enumerate(components):
                summary += 'Component ' + str(i + 1) + ' ('
                if self.is_component_complete(component):
                    summary += 'complete):\n'
                else:
                    summary += 'incomplete):\n'
                summary += '    segments:          ' + int_to_str(len(component), max_v) + '\n'
                component_len = sum(self.segments[x].get_length() for x in component)
                summary += '    length (bp):       ' + int_to_str(component_len, max_v) + '\n'
        if score:
            pad_size = len(int_to_str(max_v))
            summary += 'score:                 ' + '{:.2e}'.format(score).rjust(pad_size) + '\n'
        return summary

    def get_total_link_count(self):
        """
        Returns the total number of forward links in the graph, not counting rev comp duplicates.
        """
        links = set()
        for start, ends in self.forward_links.items():
            for end in ends:
                if (start, end) not in links and (-end, -start) not in links:
                    links.add((start, end))
        return len(links)

    def get_contig_stats(self):
        """
        Returns various contig length metrics.
        """
        segment_lengths = sorted([x.get_length() for x in self.segments.values()])
        if not segment_lengths:
            return 0, 0, 0, 0, 0, 0

        shortest = segment_lengths[0]
        longest = segment_lengths[-1]

        first_quartile_index = (len(segment_lengths) - 1) / 4
        median_index = (len(segment_lengths) - 1) / 2
        third_quartile_index = (len(segment_lengths) - 1) * 3 / 4

        first_quartile = int(round(value_from_fractional_index(segment_lengths,
                                                               first_quartile_index)))
        median = int(round(value_from_fractional_index(segment_lengths, median_index)))
        third_quartile = int(round(value_from_fractional_index(segment_lengths,
                                                               third_quartile_index)))

        half_total_length = sum(segment_lengths) / 2
        total_so_far = 0
        segment_lengths = segment_lengths[::-1]
        for length in segment_lengths:
            total_so_far += length
            if total_so_far >= half_total_length:
                n50 = length
                break
        else:
            n50 = 0

        return n50, shortest, first_quartile, median, third_quartile, longest

    def completed_circular_replicons(self):
        """
        Returns a list of graph components which are simple loops: one segment connected to itself
        to make a circular piece of DNA.
        """
        completed_components = []
        single_segment_components = [x for x in self.get_connected_components() if len(x) == 1]
        for component in single_segment_components:
            only_segment = component[0]
            if only_segment in self.forward_links and \
                    self.forward_links[only_segment] == [only_segment] and \
                    only_segment in self.reverse_links and \
                    self.reverse_links[only_segment] == [only_segment]:
                completed_components.append(only_segment)
        return completed_components

    def is_component_complete(self, component):
        """
        Given a list of unsigned segment numbers, this function returns whether the component is
        a completed circular replicon.
        """
        if len(component) != 1:
            return False
        seg = component[0]
        if self.get_downstream_seg_nums(seg) != [seg]:
            return False
        return self.get_upstream_seg_nums(seg) == [seg]

    def get_simple_path(self, starting_seg):
        """
        Starting with the given segment, this function tries to expand outward as far as possible
        while maintaining a simple (i.e. mergeable) path. If it can't expand at all, it will just
        return a list of the starting segment.
        """
        simple_path = [starting_seg]

        # Expand forward as much as possible.
        while True:
            if simple_path[-1] not in self.forward_links or \
                            len(self.forward_links[simple_path[-1]]) != 1:
                break
            potential = self.forward_links[simple_path[-1]][0]
            if potential in simple_path or -potential in simple_path:
                break
            if len(self.reverse_links[potential]) == 1 and \
                    self.reverse_links[potential][0] == simple_path[-1]:
                simple_path.append(potential)
            else:
                break

        # Expand backward as much as possible.
        while True:
            if simple_path[0] not in self.reverse_links or \
                            len(self.reverse_links[simple_path[0]]) != 1:
                break
            potential = self.reverse_links[simple_path[0]][0]
            if potential in simple_path or -potential in simple_path:
                break
            if len(self.forward_links[potential]) == 1 and \
                    self.forward_links[potential][0] == simple_path[0]:
                simple_path.insert(0, potential)
            else:
                break

        return simple_path

    def sort_link_order(self):
        """
        This function sorts the lists in links so path finding can be consistent from one run to
        the next.
        """
        for seg_num in self.forward_links:
            self.forward_links[seg_num].sort()
        for seg_num in self.reverse_links:
            self.reverse_links[seg_num].sort()

    def search(self, start, ends):
        """
        Conducts a DFS from the start segment to see if it leads to any of the end segments.
        The start segment is signed, i.e. positive start and negative start will conduct the
        search in different directions. The end segments are not signed, i.e. the search is
        successful if it reaches either orientation of an end.
        """
        end_set = set(ends)
        end_set.update(-x for x in ends)
        visited, stack = set(), [start]
        while stack:
            seg = stack.pop()
            if seg not in visited:
                visited.add(seg)
                if seg in self.forward_links:
                    for next_seg in self.forward_links[seg]:
                        if next_seg in end_set:
                            return True
                        if next_seg not in visited:
                            stack.append(next_seg)
        return False

    def get_best_paths_for_seq(self, start_seg, end_seg, target_length, sequence, scoring_scheme):
        """
        Given a sequence and target length, this function finds the best path from the start
        segment to the end segment.
        """
        # Limit the path search to lengths near the target.
        min_length = int(round(target_length * 0.8))  # TO DO: adjust or make these parameters?
        max_length = int(round(target_length * 1.2))

        # The overlap isn't present in the consensus sequence, so we need to add it on.
        if sequence:
            sequence = self.get_seq_from_signed_seg_num(start_seg)[-self.overlap:] + sequence + \
                       self.get_seq_from_signed_seg_num(end_seg)[:self.overlap]

        # If there are few enough possible paths, we just try aligning to them all.
        max_path_count = 100  # TO DO: make this a parameter?
        try:
            paths = self.all_paths(start_seg, end_seg, min_length, target_length, max_length,
                                   max_path_count)

        # If there are too many paths to try exhaustively, we use a progressive approach to find
        # the best path. The path search is done in both directions, in case one direction proves
        # more difficult than the other.
        except TooManyPaths:
            paths = self.progressive_path_find(start_seg, end_seg, min_length, target_length,
                                               max_length, sequence, scoring_scheme)
            rev_paths = self.progressive_path_find(-end_seg, -start_seg, min_length, target_length,
                                                   max_length, reverse_complement(sequence),
                                                   scoring_scheme)
            for rev_path in rev_paths:
                flipped_rev_path = [-x for x in rev_path[::-1]]
                if flipped_rev_path not in paths:
                    paths.append(flipped_rev_path)

        # We now have some paths and want to see how well the consensus sequence matches each of
        # them.
        paths_and_scores = []
        for path in paths:
            path_len = self.get_path_length(path)
            length_discrepancy = abs(path_len - target_length)

            # If there is a consensus sequence, then we actually do an alignment against the path.
            if sequence:
                path_seq = self.get_path_sequence(path)
                alignment_result = fully_global_alignment(sequence, path_seq, scoring_scheme,
                                                          True, 1000)
                if not alignment_result:
                    continue

                seqan_parts = alignment_result.split(',', 9)
                raw_score = int(seqan_parts[6])
                scaled_score = float(seqan_parts[7])

            # If there isn't a consensus sequence (i.e. the start and end overlap), then each
            # path is only scored on how well its length agrees with the target length.
            else:
                raw_score = get_num_agreement(path_len, target_length) * 100.0
                scaled_score = 100.0

            paths_and_scores.append((path, raw_score, length_discrepancy, scaled_score))

        # Sort the paths from highest to lowest quality.
        paths_and_scores = sorted(paths_and_scores, key=lambda x: (-x[1], x[2], -x[3]))

        # Only keep the paths which tie for the top slot. This will usually mean only one path,
        # but could be multiple paths that tie.
        if not paths_and_scores:
            return [], []
        best_path = paths_and_scores[0]
        best_paths = [best_path]
        for i in range(1, len(paths_and_scores)):
            potential_best = paths_and_scores[i]
            if potential_best[1:4] == best_path[1:4]:
                best_paths.append(potential_best)
            else:
                break

        return paths_and_scores, best_paths

    def get_path_availability(self, path):
        """
        Given a path, this function returns the fraction that is available. A segment is considered
        unavailable if it has no more copy depth or has a depth of below 0.5.
        """
        total_bases = 0
        available_bases = 0
        for seg_num in path:
            abs_seg_num = abs(seg_num)
            seg = self.segments[abs_seg_num]
            if abs_seg_num in self.copy_depths and len(self.copy_depths[abs_seg_num]) == 0:
                available = False
            elif seg.depth < 0.5:
                available = False
            else:
                available = True
            seg_len = seg.get_length() - self.overlap
            total_bases += seg_len
            if available:
                available_bases += seg_len
        if total_bases == 0:
            return 1.0
        else:
            return available_bases / total_bases

    def get_estimated_sequence_len(self):
        """
        Returns an estimate sequence length, based on the depth. E.g. a segment of depth 1.0
        contributes it own length (minus overlaps), a segment of depth 2.5 contributes 2.5 times
        its own length, etc.
        """
        total_seq_len = 0.0
        for seg_num, seg in self.segments.items():
            seg_len = seg.get_length()
            if seg_num in self.forward_links:
                seg_len -= self.overlap / 2
            if seg_num in self.reverse_links:
                seg_len -= self.overlap / 2
            seg_len *= seg.depth
            total_seq_len += seg_len
        return total_seq_len

    def remove_all_overlaps(self, verbosity):
        """
        This function removes all overlaps from the graph by shortening segments. It assumes that
        all overlaps in the graph are the same size.
        """
        # First we create a set of all graph edges, in both directions.
        all_edges = set()
        for start, ends in self.forward_links.items():
            for end in ends:
                all_edges.add((start, end))
                all_edges.add((-end, -start))

        # The overlap to be removed is an odd number, as SPAdes only uses odd k-mers. We'll split
        # this value approximately in half.
        large_half = int(math.ceil(self.overlap / 2))
        small_half = int(math.floor(self.overlap / 2))

        # This function will strive to put each edge into one of two groups:
        #   1) Trim more sequence from the end of the starting segment
        #   2) Trim more sequence from start of the ending segment

        # To do this, we determine which edges which must be in the same group as each other and
        # edges which must be in different groups.
        must_match = defaultdict(set)
        must_differ = defaultdict(set)

        # Firstly, each edge must be in the opposite group as its complement edge. E.g. if we
        # have an edge (5, -4) and trim more from the first segment, then we need to trim more from
        # the second segment of its complement edge (4, -5).
        for edge in all_edges:
            rev_edge = (-edge[1], -edge[0])
            must_differ[edge].add(rev_edge)
            must_differ[rev_edge].add(edge)

        # Edges which connect to the same side of a segment must be grouped together.
        pos_and_neg_seg_nums = list(self.segments) + [-x for x in self.segments]
        for seg in pos_and_neg_seg_nums:
            downstream_segs = self.get_downstream_seg_nums(seg)
            if len(downstream_segs) > 1:
                edge_1_for = (seg, downstream_segs[0])
                edge_1_rev = (-downstream_segs[0], -seg)
                for downstream_seg in downstream_segs[1:]:
                    edge_2_for = (seg, downstream_seg)
                    edge_2_rev = (-downstream_seg, -seg)
                    must_match[edge_1_for].add(edge_2_for)
                    must_match[edge_2_for].add(edge_1_for)
                    must_match[edge_1_rev].add(edge_2_rev)
                    must_match[edge_2_rev].add(edge_1_rev)
            upstream_segs = self.get_upstream_seg_nums(seg)
            if len(upstream_segs) > 1:
                edge_1_for = (upstream_segs[0], seg)
                edge_1_rev = (-seg, -upstream_segs[0])
                for upstream_seg in upstream_segs[1:]:
                    edge_2_for = (upstream_seg, seg)
                    edge_2_rev = (-seg, -upstream_seg)
                    must_match[edge_1_for].add(edge_2_for)
                    must_match[edge_2_for].add(edge_1_for)
                    must_match[edge_1_rev].add(edge_2_rev)
                    must_match[edge_2_rev].add(edge_1_rev)

        # Segments which are equal to the overlap size cannot have the larger trim applied to
        # both sides, so we require that edges on opposite sides of these segments to be grouped
        # together.
        small_seg_nums = [x for x in pos_and_neg_seg_nums
                          if self.segments[abs(x)].get_length() == self.overlap]
        for seg in small_seg_nums:
            downstream_segs = self.get_downstream_seg_nums(seg)
            upstream_segs = self.get_upstream_seg_nums(seg)
            if downstream_segs and upstream_segs:
                for downstream_seg in downstream_segs:
                    edge_1_for = (seg, downstream_seg)
                    edge_1_rev = (-downstream_seg, -seg)
                    for upstream_seg in upstream_segs:
                        edge_2_for = (upstream_seg, seg)
                        edge_2_rev = (-seg, -seg)
                        must_match[edge_1_for].add(edge_2_for)
                        must_match[edge_2_for].add(edge_1_for)
                        must_match[edge_1_rev].add(edge_2_rev)
                        must_match[edge_2_rev].add(edge_1_rev)

        # Now for the grouping! We order the edges for consistency from one run to the next.
        ordered_edges = list(all_edges)
        group_1 = set()
        group_2 = set()
        for edge in ordered_edges:
            if edge in group_1 or edge in group_2:
                continue

            # Put an edge in the first group (arbitrary decision).
            group_1.add(edge)

            while True:
                new_group_1 = set()
                new_group_2 = set()

                for group_1_edge in group_1:
                    for must_match_edge in must_match[group_1_edge]:
                        if must_match_edge in group_2 or must_match_edge in new_group_2:
                            print(must_match_edge)
                            raise CannotTrimOverlaps
                        else:
                            new_group_1.add(must_match_edge)
                    for must_differ_edge in must_differ[group_1_edge]:
                        if must_differ_edge in group_1 or must_differ_edge in new_group_1:
                            print(must_differ_edge)
                            raise CannotTrimOverlaps
                        else:
                            new_group_2.add(must_differ_edge)

                for group_2_edge in group_2:
                    for must_match_edge in must_match[group_2_edge]:
                        if must_match_edge in group_1 or must_match_edge in new_group_1:
                            print(must_match_edge)
                            raise CannotTrimOverlaps
                        else:
                            new_group_2.add(must_match_edge)
                    for must_differ_edge in must_differ[group_2_edge]:
                        if must_differ_edge in group_2 or must_differ_edge in new_group_2:
                            print(must_differ_edge)
                            raise CannotTrimOverlaps
                        else:
                            new_group_1.add(must_differ_edge)

                # We continue to loop until the groups stop growing.
                group_1_size_before = len(group_1)
                group_2_size_before = len(group_2)
                group_1.update(new_group_1)
                group_2.update(new_group_2)
                if len(group_1) == group_1_size_before and len(group_2) == group_2_size_before:
                    break

        # If the code got here, that means that all edges have been grouped according to the
        # rules, so now we produce sets of what to do for each segment. Segments in the
        # large_trim_end set will have the larger overlap trimmed from their end. Segments not in
        # that set will have the shorter overlap trimmed from their end. And similarly for the
        # large_trim_start set. These two sets are not exclusive - a segment may be in neither,
        # just one or both.
        large_trim_end = set()
        large_trim_start = set()
        for edge in group_1:
            start_seg = edge[0]
            if start_seg > 0:
                large_trim_end.add(start_seg)
            else:
                large_trim_start.add(-start_seg)
        for edge in group_2:
            end_seg = edge[1]
            if end_seg > 0:
                large_trim_start.add(end_seg)
            else:
                large_trim_end.add(-end_seg)

        # Now we finally do the segment trimming!
        if verbosity > 2:
            print('             Bases     Bases')
            print('           trimmed   trimmed')
            print(' Segment      from      from')
            print('  number     start       end')
        for seg_num, segment in self.segments.items():
            start_trim = large_half if seg_num in large_trim_start else small_half
            end_trim = large_half if seg_num in large_trim_end else small_half
            segment.trim_from_start(start_trim)
            segment.trim_from_end(end_trim)
            if verbosity > 2:
                print(str(seg_num).rjust(8) + str(start_trim).rjust(10) + str(end_trim).rjust(10))

        self.overlap = 0

    def get_downstream_seg_nums(self, seg_num):
        """
        Returns the list of downstream segments for the given segment (or an empty list if there
        are none).
        """
        if seg_num not in self.forward_links:
            return []
        else:
            return self.forward_links[seg_num]

    def get_upstream_seg_nums(self, seg_num):
        """
        Returns the list of upstream segments for the given segment (or an empty list if there
        are none).
        """
        if seg_num not in self.reverse_links:
            return []
        else:
            return self.reverse_links[seg_num]

    def remove_zero_length_segs(self, verbosity):
        """
        This function removes zero-length segments from the graph (segments with a length equal
        to the graph overlap), but only if they they aren't serving a purpose (such as in a
        multi-way junction).
        """
        segs_to_remove = []
        for seg_num, seg in self.segments.items():
            if seg.get_length() != self.overlap:
                continue
            if seg_num in self.forward_links:
                forward_connections = len(self.forward_links[seg_num])
            else:
                forward_connections = 0
            if seg_num in self.reverse_links:
                reverse_connections = len(self.reverse_links[seg_num])
            else:
                reverse_connections = 0

            if forward_connections == 1 and reverse_connections > 0:
                downstream_seg = self.forward_links[seg_num][0]
                for upstream_seg in self.reverse_links[seg_num]:
                    self.add_link(upstream_seg, downstream_seg)
                segs_to_remove.append(seg_num)

            elif reverse_connections == 1 and forward_connections > 0:
                upstream_seg = self.reverse_links[seg_num][0]
                for downstream_seg in self.forward_links[seg_num]:
                    self.add_link(upstream_seg, downstream_seg)
                segs_to_remove.append(seg_num)

        self.remove_segments(segs_to_remove)
        if verbosity > 0 and segs_to_remove:
            print('\nRemoved zero-length segments:', ', '.join(str(x) for x in segs_to_remove))

    def merge_small_segments(self, verbosity, max_merge_size):
        """
        In some cases, small segments can be merged into neighbouring segments to simplify the
        graph somewhat. This function does just that!
        """
        # This function assumes an overlap-free graph.
        if self.overlap > 0:
            return

        merged_seg_nums = []
        while True:
            for seg_num, segment in self.segments.items():
                if segment.get_length() > max_merge_size or segment.get_length() == 0:
                    continue
                downstream_segs = self.get_downstream_seg_nums(seg_num)
                upstream_segs = self.get_upstream_seg_nums(seg_num)

                # If the segment has one downstream segment and multiple upstream segments,
                # then we can merge it into the upstream segments.
                if len(downstream_segs) == 1 and len(upstream_segs) > 1 and \
                        all(self.lead_exclusively_to(x, seg_num) for x in upstream_segs):
                    for upstream_seg_num in upstream_segs:
                        upstream_seg = self.segments[abs(upstream_seg_num)]
                        if upstream_seg_num > 0:
                            upstream_seg.append_to_forward_sequence(segment.forward_sequence)
                        else:
                            upstream_seg.append_to_reverse_sequence(segment.forward_sequence)
                    segment.remove_sequence()
                    merged_seg_nums.append(seg_num)
                    break

                # If the segment has one upstream segment and multiple downstream segments,
                # then we can merge it into the downstream segments.
                if len(upstream_segs) == 1 and len(downstream_segs) > 1 and \
                        all(self.lead_exclusively_from(x, seg_num) for x in downstream_segs):
                    for downstream_seg_num in downstream_segs:
                        downstream_seg = self.segments[abs(downstream_seg_num)]
                        if downstream_seg_num > 0:
                            downstream_seg.prepend_to_forward_sequence(segment.forward_sequence)
                        else:
                            downstream_seg.prepend_to_reverse_sequence(segment.forward_sequence)
                    segment.remove_sequence()
                    merged_seg_nums.append(seg_num)
                    break
            else:
                break
            self.remove_zero_length_segs(0)

        if verbosity > 0 and merged_seg_nums:
            print('\nMerged small segments:', ', '.join(str(x) for x in merged_seg_nums))

    def starts_with_dead_end(self, signed_seg_num):
        """
        Returns whether or not the given segment (in the direction implied by the sign) is not
        lead into by any other segments.
        """
        if signed_seg_num not in self.reverse_links:
            return True
        return not self.reverse_links[signed_seg_num]

    def ends_with_dead_end(self, signed_seg_num):
        """
        Returns whether or not the given segment (in the direction implied by the sign) does not
        lead to any other segments.
        """
        if signed_seg_num not in self.forward_links:
            return True
        return not self.forward_links[signed_seg_num]


class Segment(object):
    """
    This hold a graph segment with a number, depth, direction and sequence.
    """

    def __init__(self, number, depth, sequence, positive, bridge=None, graph_path=None):
        self.number = number
        self.depth = depth
        self.forward_sequence = ''
        self.reverse_sequence = ''
        self.bridge = bridge
        self.graph_path = graph_path
        if positive:
            self.forward_sequence = sequence
        else:
            self.reverse_sequence = sequence

    def __repr__(self):
        if len(self.forward_sequence) > 6:
            seq_string = self.forward_sequence[:3] + '...' + self.forward_sequence[-3:]
        else:
            seq_string = self.forward_sequence
        return str(self.number) + ' (' + seq_string + ')'

    def add_sequence(self, sequence, positive):
        if positive:
            self.forward_sequence = sequence
        else:
            self.reverse_sequence = sequence

    def build_other_sequence_if_necessary(self):
        if not self.forward_sequence:
            self.forward_sequence = reverse_complement(self.reverse_sequence)
        if not self.reverse_sequence:
            self.reverse_sequence = reverse_complement(self.forward_sequence)

    def divide_depth(self, divisor):
        self.depth /= divisor

    def get_fastg_header(self, positive):
        """
        Returns a SPAdes-style FASTG header, without the leading '>' or ending ';'.
        """
        header = 'EDGE_' + str(self.number) + '_length_' + str(
            len(self.forward_sequence)) + '_cov_' + str(self.depth)
        if not positive:
            header += "'"
        return header

    def get_length(self):
        return len(self.forward_sequence)

    def get_length_no_overlap(self, overlap):
        return len(self.forward_sequence) - overlap

    def is_homopolymer(self):
        """
        Returns True if the segment's sequence is made up of only one base.
        """
        if len(self.forward_sequence) == 0:
            return False
        first_base = self.forward_sequence[0].lower()
        for base in self.forward_sequence[1:]:
            if base.lower() != first_base:
                return False
        return True

    def gfa_segment_line(self):
        """
        Returns an entire S line for GFA output, including the newline.
        """
        s_line = 'S\t'
        s_line += str(self.number) + '\t'
        s_line += self.forward_sequence + '\t'
        s_line += 'LN:i:' + str(self.get_length()) + '\t'
        s_line += 'DP:f:' + str(self.depth) + '\n'
        return s_line

    def save_to_fasta(self, fasta_filename):
        """
        Saves the segment's sequence to FASTA file.
        """
        fasta = open(fasta_filename, 'w')
        fasta.write('>' + self.get_fastg_header(True) + '\n')
        fasta.write(add_line_breaks_to_sequence(self.forward_sequence, 60))
        fasta.close()

    def get_seg_type_label(self):
        """
        Given a particular segment, this function returns a label string based its type.
        """
        if self.bridge is None:
            return ''
        if isinstance(self.bridge, SpadesContigBridge):
            label = 'SPAdes contig bridge'
        elif isinstance(self.bridge, LoopUnrollingBridge):
            label = 'Loop unrolling bridge'
        else:  # LongReadBridge
            label = 'Long read bridge'
        if self.graph_path:
            graph_path_str = ', '.join([str(x) for x in self.graph_path])
            label += ': ' + graph_path_str
        return label

    def trim_from_end(self, amount):
        """
        Removes the specified number of bases from the end of the segment sequence.
        """
        self.forward_sequence = self.forward_sequence[:-amount]
        self.reverse_sequence = self.reverse_sequence[amount:]

    def trim_from_start(self, amount):
        """
        Removes the specified number of bases from the end of the segment sequence.
        """
        self.forward_sequence = self.forward_sequence[amount:]
        self.reverse_sequence = self.reverse_sequence[:-amount]

    def append_to_forward_sequence(self, additional_seq):
        """
        Adds the given sequence to the end of the forward sequence (and updates the reverse
        sequence accordingly).
        """
        self.forward_sequence = self.forward_sequence + additional_seq
        self.reverse_sequence = reverse_complement(self.forward_sequence)

    def append_to_reverse_sequence(self, additional_seq):
        """
        Adds the given sequence to the end of the reverse sequence (and updates the forward
        sequence accordingly).
        """
        self.reverse_sequence = self.reverse_sequence + additional_seq
        self.forward_sequence = reverse_complement(self.reverse_sequence)

    def prepend_to_forward_sequence(self, additional_seq):
        """
        Adds the given sequence to the end of the forward sequence (and updates the reverse
        sequence accordingly).
        """
        self.forward_sequence = additional_seq + self.forward_sequence
        self.reverse_sequence = reverse_complement(self.forward_sequence)

    def prepend_to_reverse_sequence(self, additional_seq):
        """
        Adds the given sequence to the end of the reverse sequence (and updates the forward
        sequence accordingly).
        """
        self.reverse_sequence = additional_seq + self.reverse_sequence
        self.forward_sequence = reverse_complement(self.reverse_sequence)

    def remove_sequence(self):
        """
        Gets rid of the segment sequence entirely, turning it into a zero-length segment.
        """
        self.forward_sequence = ''
        self.reverse_sequence = ''

    def rotate_sequence(self, start_pos, flip, overlap):
        """
        Rotates the sequence so it begins at start_pos. If flip is True, it also switches the
        forward and reverse strands. This function assumes that the segment is a circular
        completed replicon.
        """
        unrotated_seq = self.forward_sequence
        if overlap > 0:
            unrotated_seq = unrotated_seq[:-overlap]
        rotated_seq = unrotated_seq[start_pos:] + unrotated_seq[:start_pos]
        if overlap > 0:
            rotated_seq += rotated_seq[:overlap]
        rev_comp_rotated_seq = reverse_complement(rotated_seq)

        if flip:
            self.forward_sequence = rev_comp_rotated_seq
            self.reverse_sequence = rotated_seq
        else:
            self.forward_sequence = rotated_seq
            self.reverse_sequence = rev_comp_rotated_seq


def get_error(source, target):
    """
    Returns the relative error from trying to assign the source value to the target value.
    E.g. if source = 1.6 and target = 2.0, the error is 0.2
    """
    if target > 0.0:
        return abs(source - target) / target
    else:
        return float('inf')


def within_error_margin(val_1, val_2, error_margin):
    """
    Returns whether val_1 is within the error margin of val_2.
    I.e. val_2 * (1 - em) <= val_1 <= val_2 * (1 + em)
    E.g. if val_2 is 100 and the error margin is 0.3, then val_1 must be in the range of 70 to 130
         (inclusive) for this function to return true.
    """
    return val_2 * (1 - error_margin) <= val_1 <= val_2 * (1 + error_margin)


def shuffle_into_bins(items, bins, targets):
    """
    Shuffle items into bins in all possible arrangements that satisfy these conditions:
      1) All bins must have at least one item.
      2) Any bins with a specified target must have exactly that number of items.
    """
    arrangements = []

    # If there are items not yet in a bin, place the first item in each possible bin and call this
    # function recursively.
    if items:

        # If there are only enough items to fill the empty bins, then we will only put the next
        # item in an empty bin (because putting it in a non-empty bin would prevent us from filling
        # all bins).
        empty_bin_count = sum(1 for x in bins if not x)
        only_put_in_empty = len(items) <= empty_bin_count

        for i, _ in enumerate(bins):

            # Don't put an item in a bin if that bin is already at capacity.
            if targets[i] and len(bins[i]) >= targets[i]:
                continue

            if only_put_in_empty and bins[i]:
                continue

            bins_copy = [list(x) for x in bins]
            bins_copy[i].append(items[0])
            arrangements += shuffle_into_bins(items[1:], bins_copy, targets)

    # If all items are in a bin, all bins have at least one item and any bins with a target have
    # the appropriate amount, then add the arrangement to the results.
    elif all(x for x in bins) and \
            all([not target or target == len(bins[i]) for i, target in enumerate(targets)]):
        arrangements.append(bins)
    return arrangements


def get_headers_and_sequences(filename):
    """
    Reads through a SPAdes assembly graph file and returns two lists:
    1) the headers for each segment (without the leading '>')
    2) the sequences for each segment
    """
    headers = []
    sequences = []
    header = ''
    sequence = ''
    graph_file = open(filename, 'rt')
    for line in graph_file:
        line = line.strip()
        if not line:
            continue
        if line[0] == '>':
            if header:
                headers.append(header)
                sequences.append(sequence)
                sequence = ''
            header = line[1:]
        else:
            sequence += line
    if header:
        headers.append(header)
        sequences.append(sequence)
    return headers, sequences


def get_unsigned_number_from_header(header):
    """
    Input: a SPAdes FASTG header line
    Output: an int for the segment number (always positive)
    """
    return int(header.split('_')[1])


def get_signed_number_from_header(header):
    """
    Input: a SPAdes FASTG header line
    Output: an int for the segment number (always positive)
    """
    number = get_unsigned_number_from_header(header)
    if not is_header_positive(header):
        number *= -1
    return number


def is_header_positive(header):
    """
    Input: a SPAdes FASTG header line
    Output: True if the header is for a positive segment, False for a negative segment.
    """
    if header[-1] == ';':
        header = header[:-1]
    return header.split(':')[0][-1] != "'"


def get_depth_from_header(header):
    """
    Input: a SPAdes FASTG header line
    Output: The segment's depth
    """
    header = header.split(':')[0]
    if header[-1] == "'":
        header = header[:-1]
    parts = header.split('_')
    depth_str = parts[5]
    if depth_str.endswith(';'):
        depth_str = depth_str[:-1]
    if depth_str.endswith("'"):
        depth_str = depth_str[:-1]
    return float(depth_str)


def get_links_from_header(header):
    """
    Input: a SPAdes FASTG header line
    Output: a tuple of starting segment and a list of ending segments
    """
    if header[-1] == ';':
        header = header[:-1]
    start = get_signed_number_from_header(header)
    end_list = []
    pieces = header.split(':')
    if len(pieces) > 1:
        ends = pieces[1].split(',')
        for end in ends:
            end_list.append(get_signed_number_from_header(end))
    return start, end_list


def build_rc_links_if_necessary(links):
    """
    This function makes sure that every link also has a reverse complement.  E.g. if there is a
    link from 5+ to 7-, there should also be a link from 7+ to 5-.
    """
    new_links = links.copy()
    for start, ends in links.items():
        rc_start = -start
        for end in ends:
            rc_end = -end
            if rc_end not in new_links:
                new_links[rc_end] = []
            if rc_start not in new_links[rc_end]:
                new_links[rc_end].append(rc_start)
    return new_links


def build_reverse_links(links):
    """
    This function builds a dictionary of links going the other way.  I.e. if given a dictionary
    of start to end links, it will return a dictionary of end to start links.
    """
    reverse_links = {}
    for start, ends in links.items():
        for end in ends:
            if end not in reverse_links:
                reverse_links[end] = []
            reverse_links[end].append(start)
    return reverse_links


def add_line_breaks_to_sequence(sequence, length):
    """
    Wraps sequences to the defined length.  All resulting sequences end in a line break.
    """
    seq_with_breaks = ''
    while len(sequence) > length:
        seq_with_breaks += sequence[:length] + '\n'
        sequence = sequence[length:]
    if len(sequence) > 0:
        seq_with_breaks += sequence
        seq_with_breaks += '\n'
    return seq_with_breaks


def remove_nums_from_links(links, nums_to_remove):
    """
    This function rebuilds a link dictionary excluding the given numbers.
    nums_to_remove is expected to be a list of positive (unsigned) segment numbers.
    """
    new_links = {}
    for n_1, n_2 in links.items():
        if abs(n_1) not in nums_to_remove:
            new_links[n_1] = [x for x in n_2 if abs(x) not in nums_to_remove]
            if not new_links[n_1]:
                del new_links[n_1]
    return new_links


def all_segments_are_one_base(segments):
    """
    This function returns true if all given segments have nothing but one base.
    """
    non_empty_segments = [x for x in segments if x.get_length() > 0]
    if not non_empty_segments:
        return False
    base = non_empty_segments[0].forward_sequence[0].lower()
    for segment in non_empty_segments:
        if not segment.is_homopolymer():
            return False
        forward_base = segment.forward_sequence[0].lower()
        reverse_base = segment.reverse_sequence[0].lower()
        if forward_base != base and reverse_base != base:
            return False
    return True


def is_link_positive(start, end):
    """
    Returns True if the link is 'positive'.  This is a somewhat arbitrary call that allows us to
    only get one link per RC pair.
    A link is positive if:
      1) Both segments are positive
      2) It has no RC link (i.e. is its own RC)
      3) The starting segment has a higher absolute value than the ending segment.
    """
    if start > 0 and end > 0:
        return True
    if start < 0 and end < 0:
        return False
    if start == -end:
        return True
    return abs(start) > abs(end)


def get_sign_string(num):
    """
    Returns '+' for positive numbers (and zero) and '-' for negative numbers.
    """
    if num >= 0:
        return '+'
    else:
        return '-'


def int_to_signed_string(num):
    """
    Takes an integer and returns a string with the sign at the end.
    Examples:
      5 -> 5+
      -6 -> 6-
    """
    return str(abs(num)) + get_sign_string(num)


def signed_string_to_int(signed_str):
    """
    Takes a string with the sign at the end and returns an integer.
    """
    sign = signed_str[-1]
    num = int(signed_str[:-1])
    if sign == '+':
        return num
    else:
        return -num


def insert_num_in_list(lst, val_1, val_2, insert_val):
    """
    If the list lst contains val_1 immediately followed by val_2, the function returns a new list
    with insert_val between them. If the list does not contain that sequence of values, this
    function just returns the original list.
    """
    if len(lst) < 2:
        return lst
    new_list = []
    for i, val in enumerate(lst[:-1]):
        next_val = lst[i + 1]
        new_list.append(val)
        if val == val_1 and next_val == val_2:
            new_list.append(insert_val)
    new_list.append(lst[-1])
    return new_list


def find_replace_in_list(lst, pattern, replacement):
    """
    This function looks for the given pattern in the list and if found, replaces it.
    Example: find_replace_in_list([1,5,8,3], [5,8], 7) -> [1,7,3]
    If there are multiple occurrences, it will replace them all.
    """
    replacement_made = True
    while replacement_made:
        replacement_made = False
        for i, _ in enumerate(lst):
            if lst[i] == pattern[0] and lst[i:i + len(pattern)] == pattern:
                replacement_made = True
                lst = lst[:i] + replacement + lst[i + len(pattern):]
                break
    return lst


def find_replace_one_val_in_list(lst, val, replacement):
    """
    This function looks for the given value in the list and if found, replaces it.
    Like the above function, but simpler.
    """
    if val not in lst:
        return lst
    return [replacement if x == val else x for x in lst]


def split_path(path, seg):
    """
    If val is in the list, it returns multiple lists split at that point, excluding val.
    Sort of like the string split function, but it throws out lists of 1 (because they aren't
    useful as paths).
    """
    return_paths = []
    while seg in path:
        seg_i = path.index(seg)
        return_paths.append(path[:seg_i])
        path = path[seg_i + 1:]
    return_paths.append(path)
    return_paths = [x for x in return_paths if len(x) > 1]
    return return_paths


def split_path_multiple(path, segs):
    """
    Like split_path, but vals is a list of vals, all of which split the list.
    """
    path_parts = [path]
    for seg in segs:
        new_path_parts = []
        for part in path_parts:
            new_path_parts += split_path(part, seg)
        path_parts = new_path_parts
    return path_parts


def value_from_fractional_index(lst, index):
    """
    Given a list of numbers and a fractional index, this function will interpolate between the
    values.
    """
    if not lst:
        return 0
    if len(lst) == 1:
        return lst[0]

    whole_part = int(index)
    if whole_part < 0:
        return lst[0]
    if whole_part >= len(lst) - 1:
        return lst[-1]

    fractional_part = index - float(whole_part)
    piece_1 = lst[whole_part]
    piece_2 = lst[whole_part + 1]
    return piece_1 * (1.0 - fractional_part) + piece_2 * fractional_part


def add_to_bridged_sets(start, end, right_bridged, left_bridged):
    """
    Adds the start and end segments to the sets which track bridging direction,
    based on their sign.
    """
    if start > 0:
        right_bridged.add(start)
    else:
        left_bridged.add(-start)
    if end > 0:
        left_bridged.add(end)
    else:
        right_bridged.add(-end)


def get_overlap_from_gfa_link(filename):
    """
    Looks for the first link line and gets the overlap. Assumes that all overlaps in the graph are
    the same.
    """
    gfa_file = open(filename, 'rt')
    for line in gfa_file:
        if line.startswith('L'):
            line_parts = line.strip().split('\t')
            if len(line_parts) > 5:
                cigar = line_parts[5]
                return int(cigar[:-1])
    return 0
