from __future__ import print_function
from __future__ import division
from collections import deque


class AssemblyGraph(object):
    '''
    This class holds an assembly graph with segments and links.
    '''
    def __init__(self, filename, overlap):
        self.segments = {} # Dictionary of segment number -> segment
        self.forward_links = {} # Dictionary of segment number -> segment number
        self.reverse_links = {} # Dictionary of segment number <- segment number
        self.overlap = overlap
        self.load_graph(filename)

    def load_graph(self, filename):
        # Load in the graph segments.
        headers, sequences = get_headers_and_sequences(filename)
        for i, header in enumerate(headers):
            sequence = sequences[i]
            num = get_unsigned_number_from_header(header)

            # If the segment already exists, then add this sequence.
            if num in self.segments:
                self.segments[num].add_sequence(header, sequence)

            # If the segment does not exist, make it.
            else:
                segment = Segment(header, sequence)
                self.segments[num] = segment

        # Make sure that every segment has both a forward and reverse sequence.
        for segment in self.segments.itervalues():
            segment.build_other_sequence_if_necessary()

        # Load in the links.
        for header in headers:
            start, end_list = get_links_from_header(header)
            if end_list:
                self.forward_links[start] = end_list
        self.forward_links = build_rc_links_if_necessary(self.forward_links)
        self.reverse_links = build_reverse_links(self.forward_links)

    def get_median_read_depth(self, segment_list=None):
        '''
        Returns the assembly graph's median read depth (by base).  Optionally, a list of segments
        can be given, in which case only those segments are used for the calculation.
        '''
        if not segment_list:
            segment_list = self.segments.values()
        sorted_segments = sorted(segment_list, key=lambda x: x.depth)
        total_length = 0
        for segment in sorted_segments:
            total_length += segment.get_length() - self.overlap
        halfway_length = total_length // 2
        length_so_far = 0
        for segment in sorted_segments:
            length_so_far += segment.get_length() - self.overlap
            if length_so_far >= halfway_length:
                return segment.depth
        return 0.0

    def normalise_read_depths(self):
        '''
        For every segment in the graph, divide its depth by the graph's median.
        This makes segments with the median depth have a depth of 1, segments with more than the
        median a depth of greater than 1 and segments with less than the median a depth of less
        than 1.
        '''
        median_depth = self.get_median_read_depth()
        for segment in self.segments.itervalues():
            segment.divide_depth(median_depth)

    def get_total_length(self):
        '''
        Returns the sum of all segment sequence lengths.
        '''
        total_length = 0
        for segment in self.segments.itervalues():
            total_length += segment.get_length()
        return total_length

    def save_to_fastg(self, filename):
        fastg = open(filename, 'w')
        sorted_segments = sorted(self.segments.values(), key=lambda x: x.number)
        for segment in sorted_segments:
            fastg.write(self.get_fastg_header_with_links(segment, True))
            fastg.write(add_line_breaks_to_sequence(segment.forward_sequence, 60))
            fastg.write(self.get_fastg_header_with_links(segment, False))
            fastg.write(add_line_breaks_to_sequence(segment.reverse_sequence, 60))

    def save_to_gfa(self, filename):
        gfa = open(filename, 'w')
        sorted_segments = sorted(self.segments.values(), key=lambda x: x.number)
        for segment in sorted_segments:
            gfa.write(segment.gfa_segment_line())
        gfa.write(self.get_all_gfa_link_lines())

    def get_all_gfa_link_lines(self):
        gfa_link_lines = ''
        for start, ends in self.forward_links.iteritems():
            for end in ends:
                if is_link_positive(start, end):
                    gfa_link_lines += self.gfa_link_line(start, end)
        return gfa_link_lines


    def get_fastg_header_with_links(self, segment, positive):
        '''
        Returns a full SPAdes-style FASTG header for a segment, including the leading '>', all of
        the links, the trailing ';' and a newline.
        '''
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
        '''
        Returns the total number of dead ends in the assembly graph.
        '''
        dead_ends = 0
        for segment in self.segments.itervalues():
            dead_ends += self.dead_end_count(segment)
        return dead_ends

    def dead_end_count(self, segment):
        '''
        Returns the number of dead ends for one segment: 0, 1 or 2.
        '''
        segment_num = segment.number
        dead_ends = 0
        if segment_num not in self.forward_links:
            dead_ends += 1
        if segment_num not in self.reverse_links:
            dead_ends += 1
        return dead_ends

    def filter_by_read_depth(self, cutoff):
        '''
        This function removes segments from the graph with a read depth less than the given cutoff,
        if one of the following is also true:
          1) the segment has at least one dead end
          2) the segment is part of a connected component where all of the segments are below the
             depth cutoff
          3) deleting the segment would not create any dead ends
        '''
        segment_nums_to_remove = []
        connected_components = self.get_connected_components()
        for num, segment in self.segments.iteritems():
            if segment.depth < cutoff:
                if self.dead_end_count(segment) > 0 or \
                   self.all_segments_below_depth(get_list_with_num(num, connected_components), cutoff) or \
                   not self.deleting_would_create_dead_end(segment):
                    segment_nums_to_remove.append(num)

        self.remove_segments(segment_nums_to_remove)

    def filter_homopolymer_loops(self):
        '''
        A common feature in SPAdes graphs is a small piece of the graph (often just one node) which
        has nothing but one base.  Filter these out.
        '''
        segment_nums_to_remove = []
        connected_components = self.get_connected_components()
        for component_nums in connected_components:
            component_segments = [self.segments[x] for x in component_nums]
            if all_segments_are_one_base(component_segments):
                segment_nums_to_remove += component_nums
        self.remove_segments(segment_nums_to_remove)

    def get_connected_components(self):
        '''
        Returns a list of lists, where each inner list is the segment numbers of one connected
        component of the graph.
        E.g. [[1, 2], [3, 4, 5]] would mean that segments 1 and 2 are in a connected component
        and segments 3, 4 and 5 are in another connected component. 
        '''
        visited = set()
        components = []
        for v in self.segments.iterkeys():
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
                components.append(component)
        return components

    def get_connected_segments(self, segment_num):
        '''
        Given a segment number, this function returns a list of all other segment numbers for
        segments that are directly connected.
        It only returns positive numbers (i.e. is not strand-specific).
        '''
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
        '''
        Returns true if all segments in the list are below the depth cutoff.
        '''
        for num in segment_nums:
            if self.segments[num].depth >= cutoff:
                return False
        return True

    def remove_segments(self, nums_to_remove):
        '''
        Given a list of segment numbers to remove, this function rebuilds the graph's segments
        and links, exclude those segments.
        '''
        new_segments = {}
        for num, segment in self.segments.iteritems():
            if num not in nums_to_remove:
                new_segments[num] = segment
        self.segments = new_segments
        self.forward_links = remove_nums_from_links(self.forward_links, nums_to_remove)
        self.reverse_links = remove_nums_from_links(self.reverse_links, nums_to_remove)

    def get_n_segment_length(self, n_percent):
        '''
        Returns the length for which segments that length and longer make up >= n% of the total
        bases.  E.g. if n = 50, this function returns the N50.  n must be from 0 to 100.
        '''
        total_length = self.get_total_length()
        target_length = total_length * (n_percent / 100.0)
        sorted_segments = sorted(self.segments.values(), key=lambda x: x.get_length(), reverse=True)
        length_so_far = 0
        for segment in sorted_segments:
            seg_length = segment.get_length()
            length_so_far += seg_length
            if length_so_far >= target_length:
                return seg_length
        return 0

    def gfa_link_line(self, start, end):
        '''
        Returns an entire L line for GFA output, including the newline.
        '''
        l_line = 'L\t'
        l_line += str(abs(start)) + '\t'
        l_line += get_sign_string(start) + '\t'
        l_line += str(abs(end)) + '\t'
        l_line += get_sign_string(end) + '\t'
        l_line += str(self.overlap) + 'M\n'
        return l_line

    def get_all_outputs(self, segment):
        '''
        Returns a list of segments which lead out from the given segment.
        '''
        if segment.number in self.reverse_links:
            return [self.segments[x] for x in self.forward_links[segment.number]]
        else:
            return []

    def get_exclusive_inputs(self, segment_number):
        '''
        This function finds all segments which lead into the given segment.  If those segments
        do not lead into any other segments, then this function returns them in a list.  If they
        do lead into other segments, then this function returns None.
        Specifically, this function returns a list of unsigned numbers.
        '''
        if segment_number not in self.reverse_links:
            return []
        return [abs(x) for x in self.reverse_links[segment_number] if self.lead_exclusively_to(x, segment_number)]

    def get_exclusive_outputs(self, segment_number):
        '''
        Does the same thing as get_exclusive_inputs, but in the other direction.
        '''
        if segment_number not in self.forward_links:
            return []
        return [abs(x) for x in self.forward_links[segment_number] if self.lead_exclusively_from(x, segment_number)]

    def lead_exclusively_to(self, segment_num_1, segment_num_2):
        '''
        Returns whether or not the first segment leads to and only to the second segment.
        '''
        if segment_num_1 not in self.forward_links:
            return False
        return self.forward_links[segment_num_1] == [segment_num_2]

    def lead_exclusively_from(self, segment_num_1, segment_num_2):
        '''
        Does the same thing as lead_exclusively_to, but follows links in the opposite direction.
        '''
        if segment_num_1 not in self.reverse_links:
            return False
        return self.reverse_links[segment_num_1] == [segment_num_2]

    def deleting_would_create_dead_end(self, segment):
        '''
        If deleting the given segment would create a dead end, this function returns True.
        '''
        downstream_segments = self.forward_links[segment.number]
        for downstream_segment in downstream_segments:
            if len(self.reverse_links[downstream_segment]) == 1:
                return True
        upstream_segments = self.reverse_links[segment.number]
        for upstream_segment in upstream_segments:
            if len(self.forward_links[upstream_segment]) == 1:
                return True
        return False

    def repair_four_way_junctions(self):
        '''
        This function finds and fixes four-way junctions in the graph, as these can mess up copy
        number determination. It fixes them by creating a new node with no length (i.e with the
        overlap size) to bridge the connection.
        For example: A->B,C and D->B,C becomes A->E and D->E and E->B and E->C
        '''
        seg_nums = self.segments.keys()
        seg_nums += [-x for x in self.segments.keys()]
        for seg_num in seg_nums:
            ending_segs = self.forward_links[seg_num]
            if len(ending_segs) != 2:
                continue
            end_num_1 = ending_segs[0]
            end_num_2 = ending_segs[1]
            if len(self.reverse_links[end_num_1]) != 2 or \
               len(self.reverse_links[end_num_2]) != 2:
                continue
            starting_segs = set(self.reverse_links[end_num_1])
            starting_segs.union(set(self.reverse_links[end_num_2]))
            starting_segs = list(starting_segs)
            if len(starting_segs) != 2:
                continue

            # If the code got here, then we've found a four-way junction! Create a new segment
            # to bridge the starting and ending segments.
            start_num_1 = starting_segs[0]
            start_num_2 = starting_segs[1]
            start_1 = self.segments[abs(start_num_1)]
            start_2 = self.segments[abs(start_num_2)]
            end_1 = self.segments[abs(end_num_1)]
            end_2 = self.segments[abs(end_num_2)]
            if end_1 > 0:
                bridge_seq = end_1.forward_sequence[:self.overlap]
            else:
                bridge_seq = end_1.reverse_sequence[:self.overlap]
            bridge_depth = (start_1.depth + start_2.depth + end_1.depth + end_2.depth) / 2.0
            bridge_num = self.get_next_available_seg_number()
            bridge_header = 'Node_' + str(bridge_num) + '_length_' + str(len(bridge_seq)) + \
                            '_cov_' + str(bridge_depth)
            bridge_seg = Segment(bridge_header, bridge_seq)
            self.segments[bridge_num] = bridge_seg
            self.forward_links[start_num_1] = [bridge_num]
            self.forward_links[start_num_2] = [bridge_num]
            self.forward_links[bridge_num] = [end_num_1, end_num_2]
            self.reverse_links[bridge_num] = [start_num_1, start_num_2]
            self.reverse_links[end_num_1] = [bridge_num]
            self.reverse_links[end_num_2] = [bridge_num]
            self.reverse_links[-start_num_1] = [-bridge_num]
            self.reverse_links[-start_num_2] = [-bridge_num]
            self.reverse_links[-bridge_num] = [-end_num_1, -end_num_2]
            self.forward_links[-bridge_num] = [-start_num_1, -start_num_2]
            self.forward_links[-end_num_1] = [-bridge_num]
            self.forward_links[-end_num_2] = [-bridge_num]

    def get_next_available_seg_number(self):
        '''
        This function finds the largest used segment number and returns the next 
        '''
        current_largest = max(self.segments.iterkeys())
        return current_largest + 1



class Segment(object):
    '''
    This hold a graph segment with a number, depth, direction and sequence.
    '''
    def __init__(self, header, sequence):
        self.number = 0
        self.depth = 0.0
        self.forward_sequence = ''
        self.reverse_sequence = ''
        self.parse_header(header)
        if is_header_positive(header):
            self.forward_sequence = sequence
        else:
            self.reverse_sequence = sequence

    def __repr__(self):
        if len(self.forward_sequence) > 6:
            seq_string = self.forward_sequence[:3] + '...' + self.forward_sequence[-3:]
        else:
            seq_string = self.forward_sequence
        return str(self.number) + ' (' + seq_string + ')'


    def parse_header(self, header):
        header = header[:-1]
        header = header.split(':')[0]
        if header[-1] == "'":
            header = header[:-1]
        parts = header.split('_')
        self.number = int(parts[1])
        self.depth = float(parts[5])

    def add_sequence(self, header, sequence):
        if is_header_positive(header):
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
        '''
        Returns a SPAdes-style FASTG header, without the leading '>' or ending ';'.
        '''
        header = 'EDGE_' + str(self.number) + '_length_' + str(len(self.forward_sequence)) + '_cov_' + str(self.depth)
        if not positive:
            header += "'"
        return header

    def get_length(self):
        return len(self.forward_sequence)

    def is_homopolymer(self):
        '''
        Returns True if the segment's sequence is made up of only one base.
        '''
        if len(self.forward_sequence) == 0:
            return False
        first_base = self.forward_sequence[0].lower()
        for base in self.forward_sequence[1:]:
            if base.lower() != first_base:
                return False
        return True

    def gfa_segment_line(self):
        '''
        Returns an entire S line for GFA output, including the newline.
        '''
        s_line = 'S\t'
        s_line += str(self.number) + '\t'
        s_line += self.forward_sequence + '\t'
        s_line += 'LN:i:' + str(self.get_length()) + '\t'
        s_line += 'DP:f:' + str(self.depth) + '\n'
        return s_line










def get_headers_and_sequences(filename):
    '''
    Reads through a SPAdes assembly graph file and returns two lists:
    1) the headers for each segment (without the leading '>')
    2) the sequences for each segment
    '''
    headers = []
    sequences = []
    header = ''
    sequence = ''
    graph_file = open(filename, 'r')
    for line in graph_file:
        line = line.strip()
        if not line:
            continue
        if line[0] == '>':
            if header:
                headers.append(header)
                sequences.append(sequence)
                header = ''
                sequence = ''
            header = line[1:]
        else:
            sequence += line
    if header:
        headers.append(header)
        sequences.append(sequence)
    return headers, sequences

def reverse_complement(seq):
    '''
    Given a DNA sequences, this function returns the reverse complement sequence.
    '''
    rev_comp = ''
    for i in reversed(range(len(seq))):
        rev_comp += complement_base(seq[i])
    return rev_comp

def complement_base(base):
    '''
    Given a DNA base, this returns the complement.
    '''
    forward = 'ATGCatgcRYSWKMryswkmBDHVbdhvNn.-?'
    reverse = 'TACGtacgYRSWMKyrswmkVHDBvhdbNn.-?N'
    return reverse[forward.find(base)]

def get_unsigned_number_from_header(header):
    '''
    Input: a SPAdes FASTG header line
    Output: an int for the segment number (always positive)
    '''
    return int(header.split('_')[1])

def get_signed_number_from_header(header):
    '''
    Input: a SPAdes FASTG header line
    Output: an int for the segment number (always positive)
    '''
    number = get_unsigned_number_from_header(header)
    if not is_header_positive(header):
        number *= -1
    return number

def is_header_positive(header):
    '''
    Input: a SPAdes FASTG header line
    Output: True if the header is for a positive segment, False for a negative segment.
    '''
    if header[-1] == ';':
        header = header[:-1]
    return header.split(':')[0][-1] != "'"
        
def get_links_from_header(header):
    '''
    Input: a SPAdes FASTG header line
    Output: a tuple of starting segment and a list of ending segments
    '''
    if header[-1] == ';':
        header = header[:-1]
    start = get_signed_number_from_header(header)
    end_list = []
    pieces = header.split(':')
    if len(pieces) > 1:
        ends = pieces[1].split(',')
        for end in ends:
            end_list.append(get_signed_number_from_header(end))
    return (start, end_list)

def build_rc_links_if_necessary(links):
    '''
    This function makes sure that every link also has a reverse complement.  E.g. if there is a
    link from 5+ to 7-, there should also be a link from 7+ to 5-.
    '''
    new_links = links.copy()
    for start, ends in links.iteritems():
        rc_start = -start
        for end in ends:
            rc_end = -end
            if rc_end not in new_links:
                new_links[rc_end] = []
            if rc_start not in new_links[rc_end]:
                new_links[rc_end].append(rc_start)
    return new_links

def build_reverse_links(links):
    '''
    This function builds a dictionary of links going the other way.  I.e. if given a dictionary
    of start to end links, it will return a dictionary of end to start links.
    '''
    reverse_links = {}
    for start, ends in links.iteritems():
        for end in ends:
            if end not in reverse_links:
                reverse_links[end] = []
            reverse_links[end].append(start)
    return reverse_links

def add_line_breaks_to_sequence(sequence, length):
    '''
    Wraps sequences to the defined length.  All resulting sequences end in a line break.
    '''
    seq_with_breaks = ''
    while len(sequence) > length:
        seq_with_breaks += sequence[:length] + '\n'
        sequence = sequence[length:]
    if len(sequence) > 0:
        seq_with_breaks += sequence
        seq_with_breaks += '\n'
    return seq_with_breaks

def get_list_with_num(num, num_lists):
    '''
    Given a number and a list of lists of numbers, this function returns the list of numbers
    which contains the specified number.
    E.g. given 5 and [[1,2], [3, 4, 5]], this will return [3, 4, 5]
    '''
    for num_list in num_lists:
        if num in num_list:
            return num_list
    return []

def remove_nums_from_links(links, nums_to_remove):
    '''
    This function rebuilds a link dictionary excluding the given numbers.
    nums_to_remove is expected to be a list of positive (unsigned) segment numbers.
    '''
    new_links = {}
    for n_1, n_2 in links.iteritems():
        if abs(n_1) not in nums_to_remove:
            new_links[n_1] = [x for x in n_2 if abs(x) not in nums_to_remove]
            if new_links[n_1] == []:
                del new_links[n_1]
    return new_links

def all_segments_are_one_base(segments):
    '''
    This function returns true if all given segments have nothing but one base.
    '''
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
    '''
    Returns True if the link is 'positive'.  This is a somewhat arbitrary call that allows us to
    only get one link per RC pair.
    A link is positive if:
      1) Both segments are positive
      2) It has no RC link (i.e. is its own RC)
      3) The starting segment has a higher absolute value than the ending segment.
    '''
    if start > 0 and end > 0:
        return True
    if start < 0 and end < 0:
        return False
    if start == -end:
        return True
    return abs(start) > abs(end)

def get_sign_string(num):
    '''
    Returns '+' for positive numbers (and zero) and '-' for negative numbers.
    '''
    if num >= 0:
        return '+'
    else:
        return '-'

