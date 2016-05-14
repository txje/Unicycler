
// This is the k-mer size used in the line-finding process. It is small, which leads to lots of
// background noise in the alignment rectangle, but it means that alignments will still be found
// for high-error reads.
#define KMER_SIZE 5

// Alignment lines that have an excessively small or large slope will be rejected.
#define MIN_ALLOWED_SLOPE 0.5
#define MAX_ALLOWED_SLOPE 2.0

// Reference sequences are trimmed down before conducting an actual alignment.
#define PAD_SIZE 1000

// These are the Seqan band widths used in a banded alignment.
#define STARTING_BAND_SIZE 10
#define MAX_BAND_SIZE 160


// The band size for scoring common k-mers is fixed for all sensitivity levels, as it is used
// before any line-finding.
#define COMMON_KMER_BAND_SIZE 16

// Points in an alignment line are cleaned up by throwing out the points with the most divergent
// slopes.
#define WORST_SLOPE_FRACTION_TO_DISCARD 0.05
#define WORST_SLOPE_STEPS 3

// Any CommonKmerSet where the max score falls below this value is ignored.
#define MINIMUM_MAX_SCORE 5.0

// Line finding settings come in different sensitivity levels.
#define LOW_SCORE_THRESHOLD_LEVEL_1 0.05
#define HIGH_SCORE_THRESHOLD_LEVEL_1 0.5
#define MERGE_DISTANCE_LEVEL_1 100.0
#define MIN_ALIGNMENT_LENGTH_LEVEL_1 80.0
#define MIN_POINT_COUNT_LEVEL_1 16

#define LOW_SCORE_THRESHOLD_LEVEL_2 0.025
#define HIGH_SCORE_THRESHOLD_LEVEL_2 0.25
#define MERGE_DISTANCE_LEVEL_2 100.0
#define MIN_ALIGNMENT_LENGTH_LEVEL_2 40.0
#define MIN_POINT_COUNT_LEVEL_2 16

#define LOW_SCORE_THRESHOLD_LEVEL_3 0.0125
#define HIGH_SCORE_THRESHOLD_LEVEL_3 0.125
#define MERGE_DISTANCE_LEVEL_3 100.0
#define MIN_ALIGNMENT_LENGTH_LEVEL_3 20.0
#define MIN_POINT_COUNT_LEVEL_3 8
