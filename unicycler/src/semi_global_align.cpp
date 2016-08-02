#include "semi_global_align.h"

#include <seqan/align.h>
#include <iostream>
#include <vector>
#include "settings.h"
#include <limits>
#include <algorithm>
#include <utility>




// This is the big function called by Python code. It conducts a semi-global Seqan alignment
// the given read against all references and returns the console output and all found alignments in a
// string.
char * semiGlobalAlignment(char * readNameC, char * readSeqC, int verbosity,
                           double expectedSlope, KmerPositions * refKmerPositions,
                           int matchScore, int mismatchScore, int gapOpenScore, int gapExtensionScore,
                           double lowScoreThreshold, bool returnBad, int kSize,
                           bool extraSensitive) {
    // Choose settings based on the sensitivity.
    int badLineCountSingleAlignment;
    int badLineCountEntireRead;
    int badLineCountPartialRead;
    if (extraSensitive) {
        badLineCountSingleAlignment = BAD_LINE_COUNT_SINGLE_ALIGNMENT_EXTRA_SENSITIVE;
        badLineCountEntireRead = BAD_LINE_COUNT_ENTIRE_READ_EXTRA_SENSITIVE;
        badLineCountPartialRead = BAD_LINE_COUNT_PARTIAL_READ_EXTRA_SENSITIVE;
    }
    else {
        badLineCountSingleAlignment = BAD_LINE_COUNT_SINGLE_ALIGNMENT;
        badLineCountEntireRead = BAD_LINE_COUNT_ENTIRE_READ;
        badLineCountPartialRead = BAD_LINE_COUNT_PARTIAL_READ;
    }

    // This string will collect all of the console output for the alignment.
    std::string output;

    // Change the read name and sequence to C++ strings.
    std::string readName(readNameC);
    std::string posReadName = readName + "+";
    std::string negReadName = readName + "-";
    std::string posReadSeq(readSeqC);
    std::string negReadSeq = getReverseComplement(posReadSeq);
    int readLength = posReadSeq.length();

    // Make a new KmerPositions object for the reads.
    KmerPositions readKmerPositions;
    readKmerPositions.addPositions(posReadName, posReadSeq, kSize);
    readKmerPositions.addPositions(negReadName, negReadSeq, kSize);

    // Create a CommonKmerSet for the read (both forward and reverse complement) and every reference.
    std::vector<std::string> referenceNames = refKmerPositions->getAllNames();
    std::vector<CommonKmerSet *> commonKmerSets;
    for (size_t i = 0; i < referenceNames.size(); ++i) {
        std::string refName = referenceNames[i];
        int refLength = refKmerPositions->getLength(refName);
        CommonKmerSet * forwardCommonKmerSet = new CommonKmerSet(posReadName, refName, readLength, refLength, expectedSlope, &readKmerPositions, refKmerPositions, kSize);
        commonKmerSets.push_back(forwardCommonKmerSet);
        CommonKmerSet * reverseCommonKmerSet = new CommonKmerSet(negReadName, refName, readLength, refLength, expectedSlope, &readKmerPositions, refKmerPositions, kSize);
        commonKmerSets.push_back(reverseCommonKmerSet);
    }

    if (verbosity > 2)
        output += "Seqan alignment attempts (using expected slope of " + std::to_string(expectedSlope) + ")\n";

    // We now extract alignment lines and perform alignments until we have had too many failures.
    Score<int, Simple> scoringScheme(matchScore, mismatchScore, gapExtensionScore, gapOpenScore);
    std::vector<AlignmentLine *> alignmentLines;
    std::vector<ScoredAlignment *> allAlignments;
    std::vector<ScoredAlignment *> goodAlignments;
    int badAlignmentCount = 0;
    bool oneAlignmentWholeRead = false;
    bool entireReadAligned = false;
    bool needMoreAlignments = true;

    while (needMoreAlignments) {
        // Extract an alignment line from around the highest scoring point.
        CommonKmerSet * highestScoringSet = getHighestScoringSet(commonKmerSets);
        if (highestScoringSet == 0)
            break;
        std::string readName = highestScoringSet->m_readName;
        std::string refName = highestScoringSet->m_refName;
        AlignmentLine * line = highestScoringSet->extractAlignmentLine();

        // Check to see if the alignment line is good.
        if (line == 0) {
            ++badAlignmentCount;
            if (verbosity > 2)
                output += "  line: " + refName + ", none, BAD\n";
        }
        else if (!line->buildSeedChain(MIN_POINT_COUNT, MIN_ALIGNMENT_LENGTH, kSize)) {
            ++badAlignmentCount;
            if (verbosity > 2)
                output += "  line: " + line->getDescriptiveString() + ", no seed chain, BAD\n";
            delete line;
        }
        else {
            // If the code got here, then we should have an alignment line with a seed chain ready to
            // go, so let's perform the alignment!
            if (verbosity > 2)
                output += "  line: " + line->getDescriptiveString() + ", GOOD\n";
            std::string * readSeq = readKmerPositions.getSequence(readName);
            std::string * refSeq = refKmerPositions->getSequence(refName);
            ScoredAlignment * alignment = semiGlobalAlignmentOneLine(readName, refName, readSeq, refSeq,
                                                                     line, verbosity, output, scoringScheme);
            alignmentLines.push_back(line);
            allAlignments.push_back(alignment);

            // Check to see if the alignment failed or if it scored too low.
            if (alignment == 0) {
                ++badAlignmentCount;
                if (verbosity > 2)
                    output += "    alignment: failed, BAD\n";
            }
            else if (alignment->m_scaledScore < lowScoreThreshold) {
                ++badAlignmentCount;
                if (verbosity > 2)
                    output += "    alignment: " + alignment->getShortDisplayString() + ", BAD\n";
            }

            // If the alignment is good, we add it to the results.
            else {
                goodAlignments.push_back(alignment);
                if (verbosity > 2)
                    output += "    alignment: " + alignment->getShortDisplayString() + ", GOOD\n";
                if (alignment->getReadAlignmentLength() == readLength)
                    oneAlignmentWholeRead = true;
            }

            // Finally, we want to check whether this alignment line is near any of the previous
            // alignment lines.
            AlignmentLine * bestNearbyLine = 0;
            ScoredAlignment * bestNearbyLineAlignment = 0;
            for (size_t i = 0; i < alignmentLines.size() - 1; ++i) {
                AlignmentLine * previousLine = alignmentLines[i];
                ScoredAlignment * previousLineAlignment = allAlignments[i];

                if (previousLine->isNear(line)) {

                    // If this is the first nearby line found, or...
                    if (bestNearbyLine == 0 ||

                        // this is the first nearby line found with an alignment, or...
                        (bestNearbyLineAlignment == 0 && previousLineAlignment != 0) ||

                        // this nearby line has a better score than our previous best...
                        (bestNearbyLineAlignment != 0 && previousLineAlignment != 0 &&
                         previousLineAlignment->m_scaledScore > bestNearbyLineAlignment->m_scaledScore) ) 
                    {
                        bestNearbyLine = previousLine;
                        bestNearbyLineAlignment = previousLineAlignment;
                    }
                }
            }

            // If the alignment line was indeed near another one, we try merging the two lines to
            // make a new line and aligning to that. This is good for cases of long alignments that
            // aren't well captured by a single line.
            if (bestNearbyLine != 0) {
                if (verbosity > 2) {
                    output += "  merging lines: 1) " + line->getDescriptiveString() + "\n";
                    output += "                 2) " + bestNearbyLine->getDescriptiveString() + "\n";
                }
                AlignmentLine * mergedLine = new AlignmentLine(line, bestNearbyLine);
                if (!mergedLine->buildSeedChain(MIN_POINT_COUNT, MIN_ALIGNMENT_LENGTH, kSize)) {
                    if (verbosity > 2)
                        output += "                 failed\n";
                    delete mergedLine;
                }
                else {
                    if (verbosity > 2)
                        output += "                 result: " + mergedLine->getDescriptiveString() + "\n";
                    ScoredAlignment * mergedLineAlignment = semiGlobalAlignmentOneLine(readName, refName, readSeq, refSeq,
                                                                                       mergedLine, verbosity, output, scoringScheme);
                    alignmentLines.push_back(mergedLine);
                    allAlignments.push_back(mergedLineAlignment);
                    if (mergedLineAlignment == 0) {
                        if (verbosity > 2)
                            output += "                 merged line alignment failed\n";
                    }
                    else if (mergedLineAlignment->m_scaledScore < lowScoreThreshold) {
                        if (verbosity > 2)
                            output += "                 merged line alignment: " + mergedLineAlignment->getShortDisplayString() + ", BAD\n";
                    }
                    else {
                        goodAlignments.push_back(mergedLineAlignment);
                        if (verbosity > 2)
                            output += "                 merged line alignment: " + mergedLineAlignment->getShortDisplayString() + ", GOOD\n";
                    }
                }
            }
        }

        if (oneAlignmentWholeRead)
            needMoreAlignments = (badAlignmentCount < badLineCountSingleAlignment);
        else {
            entireReadAligned = (fractionOfReadAligned(goodAlignments) == 1.0);
            if (entireReadAligned)
                needMoreAlignments = (badAlignmentCount < badLineCountEntireRead);
            else
                needMoreAlignments = (badAlignmentCount < badLineCountPartialRead);
        }
    }

    // Either all alignments or only good alignments are returned, depending on a parameter.
    std::vector<ScoredAlignment *> * returnedAlignments;
    if (returnBad)
        returnedAlignments = &allAlignments;
    else
        returnedAlignments = &goodAlignments;

    // The returned string is semicolon-delimited. The last part is the console output and the
    // other parts are alignment description strings.
    std::string returnString;
    for (size_t i = 0; i < returnedAlignments->size(); ++i) {
        ScoredAlignment * alignment = (*returnedAlignments)[i];
        if (alignment != 0)
            returnString += alignment->getFullString() + ";";
    }
    returnString += output;

    // Clean up.
    for (size_t i = 0; i < commonKmerSets.size(); ++i)
        delete commonKmerSets[i];
    for (size_t i = 0; i < allAlignments.size(); ++i) {
        if (allAlignments[i] != 0)
            delete allAlignments[i];
        delete alignmentLines[i];
    }

    return cppStringToCString(returnString);
}


 // Runs an alignment using Seqan between one read and one reference along one line.
 // It starts with a smallish band size (fast) and works up to larger ones to see if they improve
 // the alignment.
ScoredAlignment * semiGlobalAlignmentOneLine(std::string & readName, std::string & refName,
                                             std::string * readSeq, std::string * refSeq,
                                             AlignmentLine * line, int verbosity, std::string & output,
                                             Score<int, Simple> & scoringScheme) {
    long long startTime = getTime();

    int trimmedRefLength = line->m_trimmedRefEnd - line->m_trimmedRefStart;
    std::string trimmedRefSeq = refSeq->substr(line->m_trimmedRefStart, trimmedRefLength);

    Dna5String readSeqSeqan(*readSeq);
    Dna5String refSeqSeqan(trimmedRefSeq);
    int readLength = readSeq->length();

    int bandSize = STARTING_BAND_SIZE;
    ScoredAlignment * bestAlignment = 0;
    double bestAlignmentScore = std::numeric_limits<double>::min();

    // We perform the alignment with increasing band sizes until the score stops improving or we
    // reach the max band size.
    while (true) {
        ScoredAlignment * alignment = semiGlobalAlignmentOneLineOneBand(readName, refName,
                                                                        readSeqSeqan, readLength,
                                                                        refSeqSeqan, trimmedRefLength,
                                                                        line, bandSize, verbosity,
                                                                        output, scoringScheme);
        if (alignment != 0) {
            double alignmentScore = alignment->m_scaledScore;
            bool badScore = (alignmentScore <= bestAlignmentScore);
            bool tooShort = (alignment->getReadAlignmentLength() < MIN_ALIGNMENT_LENGTH) ||
                            (alignment->getRefAlignmentLength() < MIN_ALIGNMENT_LENGTH);
            if (badScore || tooShort) {
                delete alignment;
                break;
            }
            else {
                if (bestAlignment != 0)
                    delete bestAlignment;
                bestAlignment = alignment;
                bestAlignmentScore = alignmentScore;
            }
        }
        bandSize *= 2;
        if (bandSize > MAX_BAND_SIZE)
            break;
    }

    if (bestAlignment != 0)
        bestAlignment->m_milliseconds = getTime() - startTime;

    return bestAlignment;
}






// This function, given a line, will search for semi-global alignments around that line. The
// bandSize parameter specifies how far of an area around the line is searched.
ScoredAlignment * semiGlobalAlignmentOneLineOneBand(std::string & readName, std::string & refName,
                                                        Dna5String & readSeq, int readLen,
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

    ScoredAlignment * sgAlignment;
    try {
        bandedChainAlignment(alignment, line->m_bridgedSeedChain, scoringScheme, alignConfig,
                             bandSize);
        sgAlignment = new ScoredAlignment(alignment, readName, refName, readLen, refLen,
                                          line->m_trimmedRefStart, startTime, bandSize, false,
                                          false, false, scoringScheme);
    }
    catch (...) {
        if (verbosity > 2)
            output += "  Alignment failed, bandwidth = " + std::to_string(bandSize) + "\n";
        sgAlignment = 0;
    }

    return sgAlignment;
}




// This function is used to conduct a short alignment for the sake of extending a GraphMap
// alignment.
char * startExtensionAlignment(char * read, char * ref,
                               int matchScore, int mismatchScore, int gapOpenScore,
                               int gapExtensionScore) {
    long long startTime = getTime();
    std::string output;

    Dna5String sequenceH = read;
    Dna5String sequenceV = ref;
    std::string readName = "";
    std::string refName = "";

    Align<Dna5String, ArrayGaps> alignment;
    resize(rows(alignment), 2);
    assignSource(row(alignment, 0), sequenceH);
    assignSource(row(alignment, 1), sequenceV);
    Score<int, Simple> scoringScheme(matchScore, mismatchScore, gapExtensionScore, gapOpenScore);

    // The only free gaps are at the start of ref (the reference sequence).
    AlignConfig<false, true, false, false> alignConfig;
    globalAlignment(alignment, scoringScheme, alignConfig);

    ScoredAlignment startAlignment(alignment, readName, refName, length(read), length(ref),
                                   0, startTime, 0, false, true, true, scoringScheme);
    return cppStringToCString(startAlignment.getFullString());
}



// This function is used to conduct a short alignment for the sake of extending a GraphMap
// alignment.
char * endExtensionAlignment(char * read, char * ref,
                             int matchScore, int mismatchScore, int gapOpenScore,
                             int gapExtensionScore) {
    long long startTime = getTime();
    std::string output;

    Dna5String sequenceH = read;
    Dna5String sequenceV = ref;
    std::string readName = "";
    std::string refName = "";

    Align<Dna5String, ArrayGaps> alignment;
    resize(rows(alignment), 2);
    assignSource(row(alignment, 0), sequenceH);
    assignSource(row(alignment, 1), sequenceV);
    Score<int, Simple> scoringScheme(matchScore, mismatchScore, gapExtensionScore, gapOpenScore);

    // The only free gaps are at the end of ref (the reference sequence).
    AlignConfig<false, false, true, false> alignConfig;
    globalAlignment(alignment, scoringScheme, alignConfig);

    ScoredAlignment endAlignment(alignment, readName, refName, length(read), length(ref),
                                 0, startTime, 0, true, false, false, scoringScheme);
    return cppStringToCString(endAlignment.getFullString());
}



void freeCString(char * p) {
    free(p);
}


char * cppStringToCString(std::string cpp_string) {
    char * c_string = (char*)malloc(sizeof(char) * (cpp_string.size() + 1));
    std::copy(cpp_string.begin(), cpp_string.end(), c_string);
    c_string[cpp_string.size()] = '\0';
    return c_string;
}




std::string getReverseComplement(std::string sequence) {
    std::string reverseComplement;
    reverseComplement.reserve(sequence.length());
    for (int i = sequence.length() - 1; i >= 0; --i) {
        char letter = sequence[i];
        switch (letter) {
        case 'A': reverseComplement.push_back('T'); break;
        case 'T': reverseComplement.push_back('A'); break;
        case 'G': reverseComplement.push_back('C'); break;
        case 'C': reverseComplement.push_back('G'); break;
        case 'R': reverseComplement.push_back('Y'); break;
        case 'Y': reverseComplement.push_back('R'); break;
        case 'S': reverseComplement.push_back('S'); break;
        case 'W': reverseComplement.push_back('W'); break;
        case 'K': reverseComplement.push_back('M'); break;
        case 'M': reverseComplement.push_back('K'); break;
        case 'B': reverseComplement.push_back('V'); break;
        case 'D': reverseComplement.push_back('H'); break;
        case 'H': reverseComplement.push_back('D'); break;
        case 'V': reverseComplement.push_back('B'); break;
        case 'N': reverseComplement.push_back('N'); break;
        case '.': reverseComplement.push_back('.'); break;
        case '-': reverseComplement.push_back('-'); break;
        case '?': reverseComplement.push_back('?'); break;
        case '*': reverseComplement.push_back('*'); break;
        }
    }
    return reverseComplement;
}


double fractionOfReadAligned(std::vector<ScoredAlignment *> & alignments) {
    if (alignments.size() == 0)
        return true;
    std::vector<std::pair<int, int> > ranges;
    for (size_t i = 0; i < alignments.size(); ++i) {
        ScoredAlignment * alignment = alignments[i];
        int start, end;
        if (alignment->isRevComp()) {
            start = alignment->m_readLength - alignment->m_readEndPos;
            end = alignment->m_readLength - alignment->m_readStartPos;
        }
        else {
            start = alignment->m_readStartPos;
            end = alignment->m_readEndPos;
        }
        ranges.push_back(std::pair<int, int>(start, end));
    }
    std::vector<std::pair<int, int> > simplifiedRanges = simplifyRanges(ranges);
    int alignedLength = 0;
    for (size_t i = 0; i < simplifiedRanges.size(); ++i)
        alignedLength += simplifiedRanges[i].second - simplifiedRanges[i].first;
    return double(alignedLength) / alignments[0]->m_readLength;
}


std::vector<std::pair<int, int> > simplifyRanges(std::vector<std::pair<int, int> > & ranges) {
    std::sort(ranges.begin(),ranges.end());
    std::vector<std::pair<int, int> > simplifiedRanges;
    std::vector<std::pair<int, int> >::iterator it = ranges.begin();
    std::pair<int,int> current = *(it)++;
    while (it != ranges.end()){
       if (current.second >= it->first){
           current.second = std::max(current.second, it->second); 
       } else {
           simplifiedRanges.push_back(current);
           current = *(it);
       }
       it++;
    }
    simplifiedRanges.push_back(current);
    return simplifiedRanges;
}


// Given a vector of CommonKmerSet pointers, this function returns the one with the highest score.
CommonKmerSet * getHighestScoringSet(std::vector<CommonKmerSet *> & commonKmerSets) {
    double allSetsMaxScore = std::numeric_limits<double>::min();
    CommonKmerSet * bestSet = 0;
    for (size_t i = 0; i < commonKmerSets.size(); ++i) {
        double setMaxScore = commonKmerSets[i]->m_maxScore;
        if (setMaxScore > allSetsMaxScore) {
            allSetsMaxScore = setMaxScore;
            bestSet = commonKmerSets[i];
        }
    }

    // The best set has to at least reach the low threshold.
    if (allSetsMaxScore < MIN_LINE_SCORE)
        return 0;

    return bestSet;
}



