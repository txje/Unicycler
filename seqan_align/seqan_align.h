

#include <seqan/sequence.h>
#include <string>
#include "kmers.h"
#include "alignmentline.h"
#include "semiglobalalignment.h"

using namespace seqan;

// Functions that are called by the Python script must have C linkage, not C++ linkage.
extern "C" {

	char * semiGlobalAlignment(char * readNameC, char * readSeqC, char * refNameC, char * refSeqC,
	                           double expectedSlope, int verbosity, KmerPositions * kmerPositions,
	                           int matchScore, int mismatchScore, int gapOpenScore,
	                           int gapExtensionScore);

	char * startExtensionAlignment(char * read, char * ref, int readLen, int refLen, int verbosity,
	                               int matchScore, int mismatchScore, int gapOpenScore,
	                               int gapExtensionScore);

	char * endExtensionAlignment(char * read, char * ref, int readLen, int refLen, int verbosity,
	                             int matchScore, int mismatchScore, int gapOpenScore,
	                             int gapExtensionScore);

	void freeCString(char * p) {free(p);}

}

SemiGlobalAlignment * semiGlobalAlignmentOneLine(std::string & readSeq, std::string & refSeq,
                                       AlignmentLine * line, int verbosity, std::string & output,
                                       Score<int, Simple> & scoringScheme);

SemiGlobalAlignment * semiGlobalAlignmentOneLineOneBand(Dna5String & readSeq, int readLen,
                                              Dna5String & refSeq, int refLen,
                                              AlignmentLine * line, int bandSize,
                                              int verbosity, std::string & output,
                                              Score<int, Simple> & scoringScheme);

char * cppStringToCString(std::string cpp_string);
