
// This is the k-mer size used in the line-finding process. It is small, which leads to lots of
// background noise in the alignment rectangle, but it means that alignments will still be found
// for high-error reads.
#define KMER_SIZE 7

// Alignment lines that have an excessively small or large slope will be rejected.
#define MIN_ALLOWED_SLOPE 0.75
#define MAX_ALLOWED_SLOPE 1.3333

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
// either a relative error (relative to the reference alignment length) or an absolute error (in
// base pairs).
#define MAX_ALLOWED_LINE_RELATIVE_ERROR 0.2
#define MAX_ALLOWED_LINE_ABSOLUTE_ERROR 100

// We will try to find alignment lines until too many bad lines have been found. If we have aligned
// the entire read, then the threshold is lower (because we're probably done). If we haven't, then
// we try harder.
#define BAD_LINE_COUNT_SINGLE_ALIGNMENT 4
#define BAD_LINE_COUNT_ENTIRE_READ 50
#define BAD_LINE_COUNT_PARTIAL_READ 100
#define BAD_LINE_COUNT_SINGLE_ALIGNMENT_EXTRA_SENSITIVE 20
#define BAD_LINE_COUNT_ENTIRE_READ_EXTRA_SENSITIVE 250
#define BAD_LINE_COUNT_PARTIAL_READ_EXTRA_SENSITIVE 500

// Alignment lines must reach these levels to be used.
#define MIN_ALIGNMENT_LENGTH 40.0
#define MIN_POINT_COUNT 10

// A merged line will be tried for alignment lines that are closer than this.
#define ALIGNMENT_LINE_MERGE_DISTANCE 100.0
