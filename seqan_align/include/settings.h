
// This is the k-mer size used in the line-finding process. It is small, which leads to lots of
// background noise in the alignment rectangle, but it means that alignments will still be found
// for high-error reads.
#define KMER_SIZE 5

// Settings related to line finding.
#define BAND_SIZE 16
#define LOW_SCORE_THRESHOLD 2.0
#define HIGH_SCORE_THRESHOLD 20.0
#define MERGE_DISTANCE 100.0
#define MIN_ALIGNMENT_LENGTH 20.0
#define MIN_POINT_COUNT 4

// Alignment lines that have an excessively small or large slope will be rejected.
#define MIN_ALLOWED_SLOPE 0.5
#define MAX_ALLOWED_SLOPE 2.0

// Reference sequences are trimmed down before conducting an actual alignment.
#define PAD_SIZE 1000

// These are the Seqan band widths used in a banded alignment.
#define STARTING_BAND_SIZE 10
#define MAX_BAND_SIZE 160
