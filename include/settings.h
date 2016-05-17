
// This is the k-mer size used in the line-finding process. It is small, which leads to lots of
// background noise in the alignment rectangle, but it means that alignments will still be found
// for high-error reads.
#define KMER_SIZE 7

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
#define COMMON_KMER_BAND_SIZE 25
#define COMMON_KMER_BAND_THICKNESS 100

// Points in an alignment line are cleaned up by throwing out the points with the most divergent
// slopes.
#define WORST_SLOPE_FRACTION_TO_DISCARD 0.05
#define WORST_SLOPE_STEPS 3

// Line finding stops when the score drops below this value.
#define MIN_LINE_SCORE 2.0

// Line finding stops when the relative line error (line error over the aligned ref length) exceeds
// this value.
#define MAX_ALLOWED_LINE_ERROR 0.05

// We will try to find alignment lines until this number of bad lines has been found.
#define BAD_LINE_COUNT 10

// Alignment lines must reach these levels to be used.
#define MIN_ALIGNMENT_LENGTH 40.0
#define MIN_POINT_COUNT 10
