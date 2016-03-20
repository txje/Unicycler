from __future__ import print_function
from __future__ import division
from assembly_graph import AssemblyGraph


class CopyDepthGraph(AssemblyGraph):
    '''
    This class add copy depth tracking to the basic assembly graph.
    '''
    def __init__(self, filename, overlap, error_margin):
        '''
        The only argument in addition to those used in AssemblyGraph is error_margin.  This
        specifies the allowed relative depth used when making copy depth assignments.  For example,
        if error_margin = 0.2, then 20% error is acceptable.
        '''
        AssemblyGraph.__init__(self, filename, overlap)
        self.copy_depths = {} # Dictionary of segment number -> list of copy depths
        self.error_margin = error_margin

    def save_to_gfa(self, filename):
        gfa = open(filename, 'w')
        sorted_segments = sorted(self.segments.values(), key=lambda x: x.number)
        for segment in sorted_segments:
            segment_line = segment.gfa_segment_line()
            segment_line = segment_line[:-1] # Remove newline
            if segment.number in self.copy_depths:
                segment_line += '\tLB:z:' + self.get_depth_string(segment)
                segment_line += '\tCL:z:' + self.get_copy_number_colour(segment)
            segment_line += '\n'
            gfa.write(segment_line)
        gfa.write(self.get_all_gfa_link_lines())

    def get_depth_string(self, segment):
        '''
        Given a particular segment, this function returns a string with the segment's copy depths
        (if it has any).
        '''
        if segment.number not in self.copy_depths:
            return ''
        return ','.join(['%.4f' % x for x in self.copy_depths[segment.number]])

    def get_copy_number_colour(self, segment):
        '''
        Given a particular segment, this function returns a colour string based on the copy number.
        '''
        if segment.number not in self.copy_depths:
            return 'black'
        copy_number = len(self.copy_depths[segment.number])
        if copy_number == 1:
            return 'forestgreen'
        if copy_number == 2:
            return 'gold'
        if copy_number == 3:
            return 'firebrick'
        else: # 4+
            return 'red'

    def assign_copy_depth(self):
        while True:
            single_copy_count = self.assign_single_copy_segments()
            print(single_copy_count) # TEMP
            if not single_copy_count:
                break
  
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


    def assign_single_copy_segments(self):
        '''
        This function assigns a single copy to segments which satisfy all of these criteria:
          1) It has a read depth within the error margin of the median depth.  Only use segments
             without copy depth information for the median value to which we're comparing.
          2) It has a length of at least N95 for the whole graph.
          3) It has a maximum of one link connected to each end.
        It returns the number of segments for which it assigned a copy depth.
        '''
        n_95_size = self.get_n_segment_length(95)
        segments_without_copies = [x for x in self.segments.values() if x.number not in self.copy_depths]
        if not segments_without_copies:
            return 0
        median_depth = self.get_median_read_depth(segments_without_copies)
        assignment_count = 0
        for segment in segments_without_copies:
            long_enough = segment.get_length() >= n_95_size
            good_depth = within_error_margin(segment.depth, median_depth, self.error_margin)
            few_enough_links = self.at_most_one_link_per_end(segment)
            if long_enough and good_depth and few_enough_links:
                self.copy_depths[segment.number] = [segment.depth]
                assignment_count += 1
        return assignment_count
        
    def at_most_one_link_per_end(self, segment):
        '''
        Returns True if the given segment has no more than one link on either end.
        '''
        if segment in self.forward_links and len(self.forward_links[segment]) > 1:
            return False
        if segment in self.reverse_links and len(self.reverse_links[segment]) > 1:
            return False
        return True

def within_error_margin(val_1, val_2, error_margin):
    '''
    Returns whether val_1 is within the error margin of val_2.
    I.e. val_2 * (1 - em) <= val_1 <= val_2 * (1 + em)
    E.g. if val_2 is 100 and the error margin is 0.3, then val_1 must be in the range of 70 to 130
         (inclusive) for this function to return true.
    '''
    return val_1 >= val_2 * (1 - error_margin) and val_1 <= val_2 * (1 + error_margin)




