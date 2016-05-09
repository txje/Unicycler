#include "seqan_align.h"

#include <seqan/align.h>
#include <iostream>
#include <vector>
#include "settings.h"


// This is the big function called by Python code. It conducts a semi-global Seqan alignment
// the given read and reference and returns the console output and all found alignments in a
// string.
char * semiGlobalAlignment(char * readNameC, char * readSeqC, char * refNameC, char * refSeqC,
                           double expectedSlope, int verbosity, KmerPositions * kmerPositions,
                           int matchScore, int mismatchScore, int gapOpenScore,
                           int gapExtensionScore) {
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
                                                                 kmerPositions, output);
    // Now conduct an alignment for each line.
    std::vector<SemiGlobalAlignment *> alignments;
    if (lineFindingResults != 0) {
        Score<int, Simple> scoringScheme(matchScore, mismatchScore, gapExtensionScore, gapOpenScore);
        for (size_t i = 0; i < lineFindingResults->m_lines.size(); ++i) {
            AlignmentLine * line = lineFindingResults->m_lines[i];
            SemiGlobalAlignment * alignment = semiGlobalAlignmentOneLine(readSeq, refSeq, line, verbosity,
                                                                         output, scoringScheme);
            alignments.push_back(alignment);
        }
        delete lineFindingResults;
    }

    // The returned string is semicolon-delimited. The last part is the console output and the
    // other parts are alignment description strings.
    std::string returnString;
    for (size_t i = 0; i < alignments.size(); ++i) {
        returnString += alignments[i]->getFullString() + ";";
        delete alignments[i];
    }
    returnString += output;
    return cppStringToCString(returnString);
}




 // Runs an alignment using Seqan between one read and one reference along one line.
 // It starts with a smallish band size (fast) and works up to larger ones to see if they improve
 // the alignment.
SemiGlobalAlignment * semiGlobalAlignmentOneLine(std::string & readSeq, std::string & refSeq,
                                       AlignmentLine * line, int verbosity, std::string & output,
                                       Score<int, Simple> & scoringScheme) {
    long long startTime = getTime();

    int trimmedRefLength = line->m_trimmedRefEnd - line->m_trimmedRefStart;
    std::string trimmedRefSeq = refSeq.substr(line->m_trimmedRefStart, trimmedRefLength);

    Dna5String readSeqSeqan(readSeq);
    Dna5String refSeqSeqan(trimmedRefSeq);
    int readLength = readSeq.length();

    int bandSize = STARTING_BAND_SIZE;
    SemiGlobalAlignment * alignment = semiGlobalAlignmentOneLineOneBand(readSeqSeqan, readLength, refSeqSeqan, trimmedRefLength,
                                                              line, bandSize, verbosity, output, scoringScheme);

    SemiGlobalAlignment * bestAlignment = alignment;

    // Now we try larger bands to see if that improves the alignment score. We keep trying bigger
    // bands until the score stops improving or we reach the max band size.
    while (true) {
        bandSize *= 2;
        if (bandSize > MAX_BAND_SIZE)
            break;

        SemiGlobalAlignment * newAlignment = semiGlobalAlignmentOneLineOneBand(readSeqSeqan, readLength, refSeqSeqan, trimmedRefLength,
                                                                     line, bandSize, verbosity, output, scoringScheme);
        if (newAlignment->m_scaledScore <= bestAlignment->m_scaledScore) {
            delete newAlignment;
            break;
        }
        else {
            delete bestAlignment;
            bestAlignment = newAlignment;
        }
    }

    bestAlignment->m_milliseconds = getTime() - startTime;
    return bestAlignment;
}






// This function, given a line, will search for semi-global alignments around that line. The
// bandSize parameter specifies how far of an area around the line is searched.
SemiGlobalAlignment * semiGlobalAlignmentOneLineOneBand(Dna5String & readSeq, int readLen,
                                              Dna5String & refSeq, int refLen,
                                              AlignmentLine * line, int bandSize,
                                              int verbosity, std::string & output,
                                              Score<int, Simple> & scoringScheme) {
    long long startTime = getTime();

    // I encountered a Seqan crash when the band size exceeded the sequence length, so don't let
    // that happen.
    int shortestSeqLen = std::min(readLen, refLen);
    if (bandSize > shortestSeqLen)
        bandSize = shortestSeqLen;

    // The reference sequence here is the trimmed reference sequence, not the whole reference
    // sequence. But the seed chain was made using the same offset as the trimming, so everything
    // should line up nicely (no offset adjustment needed).

    Align<Dna5String, ArrayGaps> alignment;
    resize(rows(alignment), 2);
    assignSource(row(alignment, 0), readSeq);
    assignSource(row(alignment, 1), refSeq);
    AlignConfig<true, true, true, true> alignConfig;

    bandedChainAlignment(alignment, line->m_bridgedSeedChain, scoringScheme, alignConfig,
                         bandSize);

    SemiGlobalAlignment * sgAlignment = new SemiGlobalAlignment(alignment, line->m_trimmedRefStart, startTime, false, false, scoringScheme);

    if (verbosity > 2)
        output += "  Seqan alignment, bandwidth = " + std::to_string(bandSize) + ": " + sgAlignment->getShortDisplayString() + "\n";
    if (verbosity > 3)
        output += "    " + sgAlignment->m_cigar + "\n";

    return sgAlignment;
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

    SemiGlobalAlignment startAlignment(alignment, 0, startTime, false, true, scoringScheme);
    return cppStringToCString(startAlignment.getFullString());
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

    SemiGlobalAlignment endAlignment(alignment, 0, startTime, true, false, scoringScheme);
    return cppStringToCString(endAlignment.getFullString());
}



char * cppStringToCString(std::string cpp_string) {
    char * c_string = (char*)malloc(sizeof(char) * (cpp_string.size() + 1));
    std::copy(cpp_string.begin(), cpp_string.end(), c_string);
    c_string[cpp_string.size()] = '\0';
    return c_string;
}

