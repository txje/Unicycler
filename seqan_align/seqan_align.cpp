

// This is the k-mer size used in the line-finding process. It is small, which leads to lots of
// background noise in the alignment rectangle, but it means that alignments will still be found
// for high-error reads.
#define KMER_SIZE 5

// Alignment lines that have an excessively small or large slope will be rejected.
#define MIN_ALLOWED_SLOPE 0.5
#define MAX_ALLOWED_SLOPE 1.5

// Reference sequences are trimmed down before conducting an actual alignment.
#define PAD_SIZE 1000

// These are the settings relating to line finding.
#define BAND_SIZE 16
#define LOW_SCORE_THRESHOLD 2.0
#define HIGH_SCORE_THRESHOLD 20.0
#define MERGE_DISTANCE 100.0
#define MIN_ALIGNMENT_LENGTH 20.0
#define MIN_POINT_COUNT 4

// These are the Seqan band widths used in a banded alignment.
#define STARTING_BAND_SIZE 10
#define MAX_BAND_SIZE 160

typedef Align<Dna5String, ArrayGaps> TAlign;
typedef Row<TAlign>::Type TRow;
typedef Seed<Simple> TSeed;
typedef SeedSet<TSeed> TSeedSet;
typedef Row<TAlign>::Type TRow;
typedef Iterator<TRow>::Type TRowIterator;
typedef Score<int, Simple> ScoringScheme
typedef std::unordered_map<std::string, std::vector<int> > KmerPosMap;


#include <stdio.h>
#include <chrono>
#include <seqan/basic.h>
#include <seqan/align.h>
#include <seqan/sequence.h>
#include <seqan/stream.h>
#include <seqan/seeds.h>
#include <string>
#include <set>
#include <unordered_set>
#include <unordered_map>
#include <tuple>
#include <vector>
#include <map>
#include <algorithm>
#include <cmath>
#include <utility>
#include <iterator>
#include <mutex>
#include <limits>

#include "alignment.h"
#include "alignmentline.h"
#include "kmers.h"





using namespace seqan;
using namespace std::chrono;

extern "C" {

char * semiGlobalAlignment(char * readNameC, char * readSeqC, char * refNameC, char * refSeqC,
                           double expectedSlope, int verbosity, KmerPositions * kmerPositions,
                           int matchScore, int mismatchScore, int gapOpenScore,
                           int gapExtensionScore);

Alignment * semiGlobalAlignmentOneLine(Dna5String & readSeq, int readLen,
                                       Dna5String & refSeq, int refLen,
                                       AlignmentLine * line, int verbosity, std::string & output,
                                       ScoringScheme & scoringScheme);

Alignment * semiGlobalAlignmentOneLineOneBand(Dna5String & readSeq, int readLen,
                                              Dna5String & refSeq, int refLen,
                                              AlignmentLine * line, int bandSize,
                                              int verbosity, std::string & output,
                                              ScoringScheme & scoringScheme);

char * startExtensionAlignment(char * read, char * ref, int readLen, int refLen, int verbosity,
                               int matchScore, int mismatchScore, int gapOpenScore,
                               int gapExtensionScore);

char * endExtensionAlignment(char * read, char * ref, int readLen, int refLen, int verbosity,
                             int matchScore, int mismatchScore, int gapOpenScore,
                             int gapExtensionScore);

void free_c_string(char * p) {free(p);}

char * cppStringToCString(std::string cpp_string);

long long getTime() {return duration_cast< milliseconds >(system_clock::now().time_since_epoch()).count();}




// std::string vectorToString(std::vector<int> * v);
// double getMedian(std::vector<double> & v);
// void printKmerSize(int kmerSize, int locationCount, std::string & output);









// This is the big function called by Python code. It conducts a semi-global Seqan alignment
// the given read and reference and returns the console output and all found alignments in a
// string.
char * semiGlobalAlignment(char * readNameC, char * readSeqC, char * refNameC, char * refSeqC,
                           double expectedSlope, int verbosity, KmerPositions * kmerPositions) {
    // This string will collect all of the console output for the alignment.
    std::string output;

    // Change the read/ref names and sequences to C++ strings.
    std::string readName(readNameC);
    std::string refName(refNameC);
    std::string readSeq(readSeqC);
    std::string refSeq(refSeqC);
    int readLength = readSeq.length();
    int refLength = refSeq.length();

    // Find all alignment lines in the read-ref rectangle. These will be used as guides for the 
    // Seqan alignments.
    LineFindingResults * lineFindingResults = findAlignmentLines(readName, refName,
                                                                 readLength, refLength,
                                                                 expectedSlope, verbosity,
                                                                 KmerPositions, output);

    // Now conduct an alignment for each line.
    std::vector<Alignment *> alignments;
    Score<int, Simple> scoringScheme(matchScore, mismatchScore, gapExtensionScore, gapOpenScore);
    for (int i = 0; i < lineFindingResults->m_lines.size(); ++i) {
        AlignmentLine * line = lineFindingResults->m_lines[i];
        Alignment * alignment = semiGlobalAlignmentOneLine(readSeq, refSeq, line, verbosity,
                                                           output, scoringScheme);
        alignments.push_back(alignment);
    }
    delete lineFindingResults;

    // The returned string is semicolon-delimited. The last part is the console output and the
    // other parts are alignment description strings.
    std::string returnString;
    for (int i = 0; i < alignments.size(); ++i)
    {
        returnString += alignments[i]->getFullString();
        delete alignments[i];
        returnString += ";";
    }
    returnString += output;
    return cppStringToCString(returnString);
}




 // Runs an alignment using Seqan between one read and one reference along one line.
 // It starts with a smallish band size (fast) and works up to larger ones to see if they improve
 // the alignment.
Alignment * semiGlobalAlignmentOneLine(std::string & readSeq, std::string & refSeq,
                                       AlignmentLine * line, int verbosity, std::string & output,
                                       ScoringScheme & scoringScheme) {
    long long startTime = getTime();

    int trimmedRefLength = line->m_trimmedRefEnd - line->m_trimmedRefStart
    std::string trimmedRefSeq = refSeq.substr(trimmedRefStart, trimmedRefLength);

    Dna5String readSeqSeqan(readSeq);
    Dna5String refSeqSeqan(trimmedRefSeq);
    int readLength = readSeq.length();
    int refLength = refSeq.length();

    int bandSize = STARTING_BAND_SIZE;
    Alignment * alignment = semiGlobalAlignmentOneLineOneBand(readSeqSeqan, readLength, refSeqSeqan, trimmedRefLength,
                                                              line, bandSize, verbosity, output, scoringScheme);

    Alignment * bestAlignment = alignment;

    // Now we try larger bands to see if that improves the alignment score. We keep trying bigger
    // bands until the score stops improving or we reach the max band size.
    while (true) {
        bandSize *= 2;
        if (bandSize > MAX_BAND_SIZE)
            break;

        Alignment * newAlignment = semiGlobalAlignmentOneLineOneBand(readSeqSeqan, readLength, refSeqSeqan, trimmedRefLength,
                                                                     line, bandSize, verbosity, output, scoringScheme);
        if (newAlignment->m_scaledScore <= bestAlignment.scaledScore) {
            delete newAlignment;
            break;
        }
        else {
            delete bestAlignment;
            bestAlignment = newAlignment;
        }
    }

    bestAlignment.m_milliseconds = getTime() - startTime;
    return bestAlignment;
}






// This function, given a line, will search for semi-global alignments around that line. The
// bandSize parameter specifies how far of an area around the line is searched. It is generally
// much faster than exhaustiveSemiGlobalAlignment, though it may not find the optimal alignment.
// A lower bandSize is faster with a larger chance of missing the optimal alignment.
Alignment * semiGlobalAlignmentOneLineOneBand(Dna5String & readSeq, int readLen,
                                              Dna5String & refSeq, int refLen,
                                              AlignmentLine * line, int bandSize,
                                              int verbosity, std::string & output,
                                              ScoringScheme & scoringScheme) {


    // I encountered a Seqan crash when the band size exceeded the sequence length, so don't let
    // that happen.
    int shortestSeqLen = std::min(readLen, refLen);
    if (bandSize > shortestSeqLen)
        bandSize = shortestSeqLen;

    // The reference sequence here is the trimmed reference sequence, not the whole reference
    // sequence. But the seed chain was made using the same offset as the trimming, so everything
    // should line up nicely (no offset adjustment needed).

    long long startTime = getTime();

    Align<Dna5String, ArrayGaps> alignment;
    resize(rows(alignment), 2);
    assignSource(row(alignment, 0), sequenceH);
    assignSource(row(alignment, 1), sequenceV);
    AlignConfig<true, true, true, true> alignConfig;

    bandedChainAlignment(alignment, line->m_bridgedSeedChain, scoringScheme, alignConfig,
                         bandSize);

    Alignment * returnedAlignment = new Alignment(alignment, refOffset, startTime, false, false);

    if (verbosity > 2)
        output += "  Seqan alignment, bandwidth = " + std::to_string(bandSize) + ": " + alignment.getShortDisplayString() + "\n"
    if (verbosity > 3)
        output += "    " + alignment.m_cigar + "\n";

    return returnedAlignment;
}




// This function is used to conduct a short alignment for the sake of extending a GraphMap
// alignment.
char * startExtensionAlignment(char * read, char * ref, int readLen, int refLen, int verbosity,
                               int matchScore, int mismatchScore, int gapOpenScore,
                               int gapExtensionScore) {
    long long startTime = getTime();
    std::string output;

    Dna5String sequenceH = read;
    Dna5String sequenceV = ref;

    Align<Dna5String, ArrayGaps> alignment;
    resize(rows(alignment), 2);
    assignSource(row(alignment, 0), sequenceH);
    assignSource(row(alignment, 1), sequenceV);
    Score<int, Simple> scoringScheme(matchScore, mismatchScore, gapExtensionScore, gapOpenScore);

    // The only free gaps are at the start of ref (the reference sequence).
    AlignConfig<false, true, false, false> alignConfig;
    globalAlignment(alignment, scoringScheme, alignConfig);

    Alignment returnedAlignment(alignment, 0, startTime, false, true);
    return cppStringToCString(returnedAlignment->getFullString());
}



// This function is used to conduct a short alignment for the sake of extending a GraphMap
// alignment.
char * endExtensionAlignment(char * read, char * ref, int readLen, int refLen, int verbosity,
                             int matchScore, int mismatchScore, int gapOpenScore,
                             int gapExtensionScore) {
    long long startTime = getTime();
    std::string output;

    Dna5String sequenceH = read;
    Dna5String sequenceV = ref;

    Align<Dna5String, ArrayGaps> alignment;
    resize(rows(alignment), 2);
    assignSource(row(alignment, 0), sequenceH);
    assignSource(row(alignment, 1), sequenceV);
    Score<int, Simple> scoringScheme(matchScore, mismatchScore, gapExtensionScore, gapOpenScore);

    // The only free gaps are at the end of ref (the reference sequence).
    AlignConfig<false, false, true, false> alignConfig;
    globalAlignment(alignment, scoringScheme, alignConfig);

    Alignment returnedAlignment(alignment, 0, startTime, true, false);
    return cppStringToCString(returnedAlignment->getFullString());
}



char * cppStringToCString(std::string cpp_string) {
    char * c_string = (char*)malloc(sizeof(char) * (cpp_string.size() + 1));
    std::copy(cpp_string.begin(), cpp_string.end(), c_string);
    c_string[cpp_string.size()] = '\0';
    return c_string;
}



// std::string vectorToString(std::vector<int> * v) {
//     std::stringstream ss;
//     for(size_t i = 0; i < v->size(); ++i) {
//         if (i != 0)
//             ss << ",";
//         ss << (*v)[i];
//     }
//     return ss.str();
// }

// // Gets the median from an already-sorted vector.
// double getMedian(std::vector<double> & v) {
//     size_t size = v.size();
//     if (size % 2 == 0)
//         return (v[size / 2 - 1] + v[size / 2]) / 2;
//     else 
//         return v[size / 2];
// }


// void printKmerSize(int kmerSize, int locationCount, std::string & output) {
//     output += "  " + std::to_string(locationCount) + " " + std::to_string(kmerSize) + "-mers in common\n";
// }




}
