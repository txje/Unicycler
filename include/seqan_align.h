

#include <seqan/sequence.h>
#include <string>
#include <vector>
#include "kmers.h"
#include "commonkmerset.h"
#include "alignmentline.h"
#include "semiglobalalignment.h"

using namespace seqan;

// Functions that are called by the Python script must have C linkage, not C++ linkage.
extern "C" {

    char * semiGlobalAlignment(char * readNameC, char * readSeqC, int verbosity,
                               double expectedSlope, KmerPositions * kmerPositions,
                               int matchScore, int mismatchScore, int gapOpenScore, int gapExtensionScore,
                               double lowScoreThreshold);

    // char * semiGlobalAlignment(char * readNameC, char * readSeqC, char * refNameC, char * refSeqC,
    //                            double expectedSlope, int verbosity, KmerPositions * kmerPositions,
    //                            int matchScore, int mismatchScore, int gapOpenScore,
    //                            int gapExtensionScore, int sensitivityLevel);

    char * startExtensionAlignment(char * read, char * ref,
                                   int matchScore, int mismatchScore, int gapOpenScore,
                                   int gapExtensionScore);

    char * endExtensionAlignment(char * read, char * ref,
                                 int matchScore, int mismatchScore, int gapOpenScore,
                                 int gapExtensionScore);

    void freeCString(char * p) {free(p);}

    char * getRandomSequenceAlignmentScores(int seqLength, int n,
                                            int matchScore, int mismatchScore, int gapOpenScore, int gapExtensionScore);
}

SemiGlobalAlignment * semiGlobalAlignmentOneLine(std::string & readName, std::string & refName,
                                                 std::string * readSeq, std::string * refSeq,
                                                 AlignmentLine * line, int verbosity, std::string & output,
                                                 Score<int, Simple> & scoringScheme);

SemiGlobalAlignment * semiGlobalAlignmentOneLineOneBand(std::string & readName, std::string & refName,
                                                        Dna5String & readSeq, int readLen,
                                                        Dna5String & refSeq, int refLen,
                                                        AlignmentLine * line, int bandSize,
                                                        int verbosity, std::string & output,
                                                        Score<int, Simple> & scoringScheme);

SemiGlobalAlignment * fullyGlobalAlignment(std::string s1, std::string s2,
                                           int matchScore, int mismatchScore, int gapOpenScore, int gapExtensionScore);

char * cppStringToCString(std::string cpp_string);

std::string getReverseComplement(std::string sequence);

double fractionOfReadAlignedOverThreshold(std::vector<SemiGlobalAlignment *> & alignments, double scoreThreshold);

double fractionOfReadAligned(std::vector<SemiGlobalAlignment *> & alignments);

std::vector<std::pair<int, int> > simplifyRanges(std::vector<std::pair<int, int> > & ranges);

std::string getRandomSequence(int seqLength, std::mt19937 & gen, std::uniform_int_distribution<double> & dist);

char getRandomBase(std::mt19937 & gen, std::uniform_int_distribution<double> & dist);

CommonKmerSet * getHighestScoringSet(std::vector<CommonKmerSet *> & commonKmerSets);


