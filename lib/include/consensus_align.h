
#ifndef CONSENSUS_ALIGN_H
#define CONSENSUS_ALIGN_H

#include <seqan/basic.h>
#include <seqan/score.h>
#include <seqan/consensus.h>

using namespace seqan;

// Functions that are called by the Python script must have C linkage, not C++ linkage.
extern "C" {

	char * multipleSequenceAlignment(char * sequences[], char * qualities[], int sequenceCount,
		                             int matchScore, int mismatchScore, int gapOpenScore, int gapExtensionScore);

}


char getMostCommonBase(std::vector<char> & bases, std::vector<char> & qualities);


#endif // CONSENSUS_ALIGN_H
