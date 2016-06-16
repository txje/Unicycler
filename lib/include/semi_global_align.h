#ifndef SEMI_GLOBAL_ALIGN_H
#define SEMI_GLOBAL_ALIGN_H

#include <seqan/sequence.h>
#include <string>
#include <vector>
#include "kmers.h"
#include "commonkmerset.h"
#include "alignmentline.h"
#include "scoredalignment.h"
#include "random_alignments.h"

using namespace seqan;

// Functions that are called by the Python script must have C linkage, not C++ linkage.
extern "C" {

    char * semiGlobalAlignment(char * readNameC, char * readSeqC, int verbosity,
                               double expectedSlope, KmerPositions * refKmerPositions,
                               int matchScore, int mismatchScore, int gapOpenScore, int gapExtensionScore,
                               double lowScoreThreshold, bool returnBad, int kSize);

    char * startExtensionAlignment(char * read, char * ref,
                                   int matchScore, int mismatchScore, int gapOpenScore,
                                   int gapExtensionScore);

    char * endExtensionAlignment(char * read, char * ref,
                                 int matchScore, int mismatchScore, int gapOpenScore,
                                 int gapExtensionScore);

    void freeCString(char * p);
}

ScoredAlignment * semiGlobalAlignmentOneLine(std::string & readName, std::string & refName,
                                             std::string * readSeq, std::string * refSeq,
                                             AlignmentLine * line, int verbosity, std::string & output,
                                             Score<int, Simple> & scoringScheme);

ScoredAlignment * semiGlobalAlignmentOneLineOneBand(std::string & readName, std::string & refName,
                                                    Dna5String & readSeq, int readLen,
                                                    Dna5String & refSeq, int refLen,
                                                    AlignmentLine * line, int bandSize,
                                                    int verbosity, std::string & output,
                                                    Score<int, Simple> & scoringScheme);

char * cppStringToCString(std::string cpp_string);

std::string getReverseComplement(std::string sequence);

double fractionOfReadAligned(std::vector<ScoredAlignment *> & alignments);

std::vector<std::pair<int, int> > simplifyRanges(std::vector<std::pair<int, int> > & ranges);

CommonKmerSet * getHighestScoringSet(std::vector<CommonKmerSet *> & commonKmerSets);

#endif // SEMI_GLOBAL_ALIGN_H
