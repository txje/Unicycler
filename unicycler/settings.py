"""
Various hard-coded settings for Unicycler.

Author: Ryan Wick
email: rrwick@gmail.com
"""

# Unicycler will only work with read alignments if they are long enough. This values specifies
# the threshold relative to the graph overlap. E.g. if this value is 2 and the graph used was a
# SPAdes 95-mer graph, it would have an overlap of 95 bp and so the minimum used alignment would
# be 190 bp.
MIN_ALIGNMENT_LENGTH_RELATIVE_TO_GRAPH_OVERLAP = 2

# If less than this fraction of a read was aligned in the first aligning pass, Unicycler will try
# again using much more sensitive alignment settings. This helps to align reads which come from
# particularly difficult repetitive regions.
MIN_READ_FRACTION_ALIGNED = 0.9

# This is how much overlap is allowed between two alignments in a single read, relative to the
# graph's overlap. For example, if the graph has an overlap of 95 and this value is 1.1,
# then alignments within a read can go up to 105 bp, but alignments with more overlap will be
# filtered out.
ALLOWED_ALIGNMENT_OVERLAP = 1.1

# Unicycler will not use the lowest quality alignments for making bridges. This setting specifies
# the threshold. E.g. if it is 5, then any alignment with a scaled score of less than the 5th
# percentile scaled score will be thrown out.
MIN_SCALED_SCORE_PERCENTILE = 5.0

# When Unicycler is searching for paths connecting two graph segments which matches a read
# consensus sequence, it will only consider paths which have a length similar to the consensus
# sequence. These settings define the acceptable range. E.g. if they are 0.7 and 1.3, Unicycler
# will consider graph paths that range from 70% to 130% of the consensus sequence length.
MIN_RELATIVE_PATH_LENGTH = 0.8
MAX_RELATIVE_PATH_LENGTH = 1.2

# These settings are used when Unicycler is exhaustively searching for paths connecting two graph
# segments. If the number of working paths or the number of final paths gets too high during the
# search (exceeds these thresholds), Unicycler will give up and instead try a progressive path
# search.
ALL_PATH_SEARCH_MAX_WORKING_PATHS = 10000
ALL_PATH_SEARCH_MAX_FINAL_PATHS = 250

# These settings are used when Unicycler is progressively searching for paths connecting two graph
# segments. When its number of working paths reaches PROGRESSIVE_PATH_SEARCH_MAX_WORKING_PATHS, it
# will cull them down to the PROGRESSIVE_PATH_SEARCH_KEEP_COUNT best paths.
PROGRESSIVE_PATH_SEARCH_MAX_WORKING_PATHS = 100
PROGRESSIVE_PATH_SEARCH_KEEP_COUNT = 10

# These settings are used for Unicycler's copy number determination - the process by which it
# tries to figure out the depth of constituent components of each segment.
#   * INITIAL_SINGLE_COPY_TOLERANCE controls how much excess depth is acceptable for the first
#     single-copy assignment pass.
#   * COPY_PROPAGATION_TOLERANCE controls how much discrepancy is allowed when propagating copy
#     depths from one segment to the next.
#   * MIN_SINGLE_COPY_LENGTH is how short of a segment can be called single-copy when adding
#     additional single copy segments.
#   * MIN_HALF_MEDIAN_FOR_DIPLOID is used when determining whether a graph should be considered
#     diploid or not for copy depths. At least this fraction of the bases must be near in depth to
#     half of the graph's median value in order for the single-copy depth to be shifted down to
#     0.5.
#   * MAX_COPY_DEPTH_DISTRIBUTION_ARRANGEMENTS caps the number of possible ways to redistribute a
#     segment's copy depths to its neighbours. If there are more possibilities than this,
#     Unicycler won't bother trying.
INITIAL_SINGLE_COPY_TOLERANCE = 0.1
COPY_PROPAGATION_TOLERANCE = 0.5
MIN_SINGLE_COPY_LENGTH = 1000
MIN_HALF_MEDIAN_FOR_DIPLOID = 0.1
MAX_COPY_DEPTH_DISTRIBUTION_ARRANGEMENTS = 10000

# When Unicycler is cleaning up the graph after bridging, it can delete graph paths and graph
# components which are mostly (but not entirely) used up in bridges. This value controls the
# threshold. Making it smaller will make Unicycler more willing to delete stuff. Making it larger
# will make Unicycler less willing to delete stuff.
CLEANING_USEDUPNESS_THRESHOLD = 0.5

# When making a consensus sequence, Unicycler will use up to this many read sequences. It is
# capped because a consensus with too many reads will take too long and likely not be much
# better. E.g. a 100 read consensus will likely give a similar sequence as a 25 read consensus,
# but it will take much longer.
MAX_READS_FOR_CONSENSUS = 25
