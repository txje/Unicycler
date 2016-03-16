from __future__ import print_function
from __future__ import division
from collections import deque


class AssemblyGraph:
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

    def get_median_read_depth(self):
        '''
        Returns the assembly graph's median read depth (by base).
        '''
        sorted_segments = sorted(self.segments.values(), key=lambda x: x.depth)
        total_length = self.get_total_length()
        halfway_length = total_length // 2
        length_so_far = 0
        for segment in sorted_segments:
            length_so_far += segment.get_length()
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
        '''
        segment_nums_to_remove = []
        connected_components = self.get_connected_components()
        for num, segment in self.segments.iteritems():
            if segment.depth < cutoff:
                if self.dead_end_count(segment) > 0:
                    segment_nums_to_remove.append(num)
                else:
                    component = get_list_with_num(num, connected_components)
                    if self.all_segments_below_depth(component, cutoff):
                        segment_nums_to_remove.append(num)
        self.remove_segments(segment_nums_to_remove)

    def filter_homopolymer_loops(self):
        '''
        A common feature in SPAdes graphs is a k+1 segment that is all one base and loops to itself
        (and doesn't have any other connections).  This function removes these from the graph.
        '''
        segment_nums_to_remove = []
        for num, segment in self.segments.iteritems():
            if segment.get_length() == self.overlap + 1 and \
               segment.is_homopolymer() and \
               self.forward_links[num] == [num] and \
               self.forward_links[-num] == [-num]:
                segment_nums_to_remove.append(num)
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

    def get_n_segment_length(self, n):
        '''
        Returns the length for which segments that length and longer make up >= n% of the total
        bases.  E.g. if n = 50, this function returns the N50.
        n must be from 0 to 100.
        '''
        total_length = self.get_total_length()
        target_length = total_length * (n / 100.0)
        sorted_segments = sorted(self.segments.values(), key=lambda x: x.get_length(), reverse=True)
        length_so_far = 0
        for segment in sorted_segments:
            seg_length = segment.get_length()
            length_so_far += seg_length
            if length_so_far >= target_length:
                return seg_length
        return 0




class Segment:
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
        base = seq[i]
        if base == 'A': rev_comp += 'T'
        elif base == 'T': rev_comp += 'A'
        elif base == 'G': rev_comp += 'C'
        elif base == 'C': rev_comp += 'G'
        elif base == 'a': rev_comp += 't'
        elif base == 't': rev_comp += 'a'
        elif base == 'g': rev_comp += 'c'
        elif base == 'c': rev_comp += 'g'
        elif base == 'R': rev_comp += 'Y'
        elif base == 'Y': rev_comp += 'R'
        elif base == 'S': rev_comp += 'S'
        elif base == 'W': rev_comp += 'W'
        elif base == 'K': rev_comp += 'M'
        elif base == 'M': rev_comp += 'K'
        elif base == 'r': rev_comp += 'y'
        elif base == 'y': rev_comp += 'r'
        elif base == 's': rev_comp += 's'
        elif base == 'w': rev_comp += 'w'
        elif base == 'k': rev_comp += 'm'
        elif base == 'm': rev_comp += 'k'
        elif base == 'B': rev_comp += 'V'
        elif base == 'D': rev_comp += 'H'
        elif base == 'H': rev_comp += 'D'
        elif base == 'V': rev_comp += 'B'
        elif base == 'b': rev_comp += 'v'
        elif base == 'd': rev_comp += 'h'
        elif base == 'h': rev_comp += 'd'
        elif base == 'v': rev_comp += 'b'
        elif base == 'N': rev_comp += 'N'
        elif base == 'n': rev_comp += 'n'
        elif base == '.': rev_comp += '.'
        elif base == '-': rev_comp += '-'
        elif base == '?': rev_comp += '?'
        else: rev_comp += 'N'
    return rev_comp

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
    return new_links

