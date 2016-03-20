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
            segment_line = segment_line[:-1]
            ADD LABEL AND COLOUR TAGS HERE
            segment_line += '\n'
            gfa.write(segment_line)
        gfa.write(self.get_all_gfa_link_lines())

    def assign_single_copy_segments(self):
        '''
        This function assigns a single copy to segments which satisfy all of these criteria:
          1) It has a read depth within the error margin of the median depth.  Only use segments
             without copy depth information for the median value to which we're comparing.
          2) It has a length of at least N90 for the whole graph.
          3) It has a maximum of one link connected to each end.
        It returns the number of segments for which it assigned a copy depth.
        '''
        n_90_size = self.get_n_segment_length(90)
        segments_without_copies = [x for x in self.segments.values() if x not in self.copy_depths]
        if not segments_without_copies:
            return 0
        median_depth = self.get_median_read_depth(segments_without_copies)
        assignment_count = 0
        for segment in segments_without_copies:
            if segment.get_length() >= n_90_size and \
               within_error_margin(segment.depth, median_depth, self.error_margin) and \
               self.at_most_one_link_per_end(segment):
                self.copy_depths[segment.number] = [segment.depth]
                assignment_count += 0
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
    return val_1 >= val_2 * (1 - error_margin) and val_1 <= val_2 * (1 - error_margin)




