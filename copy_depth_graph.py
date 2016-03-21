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
        return ', '.join(['%.3f' % x for x in self.copy_depths[segment.number]])

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
            return 'darkorange'
        else: # 4+
            return 'red'

    def determine_copy_depth(self):
        '''
        This function iteratively applies the various methods for assigning copy depth to segments
        until no more assignments can be made.  It may not succeed in assigning copy depths to all
        segments, as some segments will have strange/difficult connections or depth which prevent
        automatic copy depth determination.
        '''
        while self.assign_single_copy_depths():
            self.determine_copy_depth_part_2()

    def determine_copy_depth_part_2(self):
        '''
        This function will recursively call itself to follow these rules:
           1) Run merge_copy_depths until it stops successfully assigning copy depths.
           2) Run redistribute_copy_depths.  If it succeeded in assigning any copy depths, go back
              to step 1.
           3) Run simple_loop_copy_depths.  If it succeeded in assigning any copy depths, go back
              to step 1.
        When this function completes, it means that no more copy depths can be assigned using those
        three functions.
        '''
        while self.merge_copy_depths():
            pass
        if self.redistribute_copy_depths():
            self.determine_copy_depth_part_2()
        if self.simple_loop_copy_depths():
            self.determine_copy_depth_part_2()

    def assign_single_copy_depths(self):
        '''
        This function assigns a single copy to segments which satisfy all of these criteria:
          1) It has a read depth within the error margin of the median depth.  Only use segments
             without copy depth information for the median value to which we're comparing.
          2) It has a length of at least N95 for the whole graph.
          3) It has a maximum of one link connected to each end.
        It returns the number of segments for which it assigned a copy depth.
        '''
        n_95_size = self.get_n_segment_length(95)
        segments = self.get_segments_without_copies()
        if not segments:
            return 0
        median_depth = self.get_median_read_depth(segments)
        assignment_count = 0
        for segment in segments:
            long_enough = segment.get_length() >= n_95_size
            good_depth = within_error_margin(segment.depth, median_depth, self.error_margin)
            few_enough_links = self.at_most_one_link_per_end(segment)
            if long_enough and good_depth and few_enough_links:
                self.copy_depths[segment.number] = [segment.depth]
                assignment_count += 1
        print('assign_single_copy_segments:', assignment_count) # TEMP
        return assignment_count

    def merge_copy_depths(self):
        '''
        This function looks for segments where they have input on one end where:
          1) All input segments have copy depth assigned.
          2) All input segments exclusively input to this segment.
        In these cases, if the sum of the input copy depths is within the error margin of the
        segment's depth, then we assign copy depths to the segment, scaling the inputs so their sum
        exactly matches the segment's depth.
        If both ends can potentially be used to assign copy depths to a node, we use the side which
        more closely matches the node's depth.
        '''
        segments = self.get_segments_without_copies()
        if not segments:
            return 0
        assignment_count = 0
        for segment in segments:
            num = segment.number
            exclusive_inputs = self.get_exclusive_inputs(num)
            exclusive_outputs = self.get_exclusive_outputs(num)
            in_depth_possible = exclusive_inputs and self.all_have_copy_depths(exclusive_inputs)
            out_depth_possible = exclusive_outputs and self.all_have_copy_depths(exclusive_outputs)
            in_depth_acceptable = False
            out_depth_acceptable = False
            if in_depth_possible:
                in_depths, in_error = self.scale_copy_depths(num, exclusive_inputs)
                in_depth_acceptable = in_error <= self.error_margin
            if out_depth_possible:
                out_depths, out_error = self.scale_copy_depths(num, exclusive_outputs)
                out_depth_acceptable = out_error <= self.error_margin
            if in_depth_acceptable or out_depth_acceptable:
                assignment_count += 1
                if in_depth_acceptable and (not out_depth_acceptable or in_error < out_error):
                    self.copy_depths[num] = in_depths
                else:
                    self.copy_depths[num] = out_depths
        print('merge_copy_depths:', assignment_count) # TEMP
        return assignment_count

    def redistribute_copy_depths_easy(self):
        '''
        This function deals with the easier case of copy depth redistribution: where one segments
        with copy depth leads exclusively to multiple segments without copy depth.
        We will then try to redistribute the source segment's copy depths among the destination
        segments.  If it can be done within the allowed error margin, the destination segments will
        get their copy depths.
        '''
        assignment_count = 0

        # TO DO
        # TO DO
        # TO DO
        #
        # LOOP THROUGH ALL SEGMENTS, LOOKING FOR THE FIRST INSTANCE OF A SEGMENT WITH DEPTH WHICH EXCLUSIVELY LEADS (IN EITHER DIRECTION) TO MULTIPLE SEGMENTS, AT LEAST ONE OF WHICH LACKS DEPTH.
        # IF THE SOURCE COPY NUMBER IS LESS THAN THE DESTINATION SEGMENT COUNT, IT ISN'T POSSIBLE SO CONTINUE.
        # SHUFFLE THE SOURCE DEPTHS INTO BINS, ONE FOR EACH DESTINATION SEGMENT, IN ALL POSSIBLE COMBINATIONS.
        # ASSESS THE ERROR OF EACH COMBINATION.  FOR DESTINATION SEGMENTS THAT ALREADY HAVE COPY DEPTHS, WE CAN COMPARE THE SOURCE DEPTHS TO THE DESTINATION DEPTHS DIRECTLY.  FOR DESTINATION SEGMENTS WITHOUT COPY DEPTHS, WE COMPARE THE PROPOSED SUM TO THE SEGMENT'S READ DEPTH.
        # CHOOSE THE LOWEST ERROR COMBINATIOIN.  IF ITS ERROR IS WITHIN THE ALLOWED ERROR MARGIN, SCALE THE COPY DEPTHS AND ASSIGN THEM.
        # THEN LEAVE THIS FUNCTION (I.E. ONLY DO ONE REDISTRIBUTION)
        #    
        # TO DO
        # TO DO
        # TO DO

        print('redistribute_copy_depths:', assignment_count) # TEMP
        return assignment_count

    def redistribute_copy_depths_complex(self):
        '''
        This function deals with a more complex case of copy depth redistribution: where multiple
        segments with copy depth lead exclusively to multiple segments without copy depth.
        We will then try to redistribute the source segment's copy depths among the destination
        segments.  If it can be done within the allowed error margin, the destination segments will
        get their copy depths.
        '''
        assignment_count = 0

        # TO DO
        # TO DO
        # TO DO
        #
        # LOOP THROUGH ALL SEGMENTS, LOOKING FOR ANY INSTANCES OF A SEGMENT WITH COPY DEPTHS WHICH CONNECTS TO A SEGMENT WITHOUT COPY DEPTHS.
        # EXPAND THE SOURCE AND DESTINATION GROUPS UNTIL EITHER WE GET A SOURCE NODE LACKING COPY DEPTHS (FAILURE) OR THEY STOP GROWING.
        # NOW WE SHOULD HAVE A SOURCE GROUP, ALL OF WHICH HAVE COPY DEPTHS, AND A DESTINATION GROUP, AT LEAST SOME OF WHICH DO NOT HAVE COPY DEPTHS.
        # DO A SIMILAR COMBINATORIAL APPROACH AS WITH THE EASY CASE.
        #
        # TO DO
        # TO DO
        # TO DO

        print('redistribute_copy_depths:', assignment_count) # TEMP
        return assignment_count

    def simple_loop_copy_depths(self):
        '''
        This function assigns copy depths to simple loop structures.  It will only assign copy
        depths in cases where the loop occurs once - higher repetition loops will not be given copy
        depths due to the increasing uncertainty in repetition counts.
        '''
        assignment_count = 0
        # TO DO
        # TO DO
        # TO DO
        # TO DO
        # TO DO
        # TO DO
        print('simple_loop_copy_depths:', assignment_count) # TEMP
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

    def all_have_copy_depths(self, segment_numbers):
        '''
        Takes a list of segment numbers and returns whether every segment in the list has copy
        depths assigned.
        '''
        for num in segment_numbers:
            if num not in self.copy_depths:
                return False
        return True

    def scale_copy_depths(self, segment_number, source_segment_numbers):
        '''
        Using a list of segments which are the source of copy depth, this function scales them so
        that their sum matches the depth of the given segment.
        It returns:
          1) a list of depth numbers
          2) the error (i.e. the degree of scaling which had to occur)
        It assumes that all of the source segments definitely have copy depths.
        '''
        source_depths = []
        for num in source_segment_numbers:
            source_depths += self.copy_depths[num]
        source_depth_sum = sum(source_depths)
        target_depth = self.segments[segment_number].depth
        scaling_factor = target_depth / source_depth_sum
        scaled_depths = sorted([scaling_factor * x for x in source_depths], reverse=True)
        if target_depth > 0.0:
            error = abs(source_depth_sum - target_depth) / target_depth
        else:
            error = 1.0
        return scaled_depths, error

    def get_segments_without_copies(self):
        '''
        Returns a list of the graph segments lacking copy depth information.
        '''
        return [x for x in self.segments.values() if x.number not in self.copy_depths]


def within_error_margin(val_1, val_2, error_margin):
    '''
    Returns whether val_1 is within the error margin of val_2.
    I.e. val_2 * (1 - em) <= val_1 <= val_2 * (1 + em)
    E.g. if val_2 is 100 and the error margin is 0.3, then val_1 must be in the range of 70 to 130
         (inclusive) for this function to return true.
    '''
    return val_1 >= val_2 * (1 - error_margin) and val_1 <= val_2 * (1 + error_margin)




