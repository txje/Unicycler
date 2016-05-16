#include "seqan_align.h"

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
                           double expectedSlope, KmerPositions * kmerPositions,
                           int matchScore, int mismatchScore, int gapOpenScore, int gapExtensionScore,
                           double lowScoreThreshold) {
    // This string will collect all of the console output for the alignment.
    std::string output;

    // Change the read name and sequence to C++ strings.
    std::string readName(readNameC);
    std::string posReadName = readName + "+";
    std::string negReadName = readName + "-";
    std::string posReadSeq(readSeqC);
    std::string negReadSeq = getReverseComplement(posReadSeq);
    int readLength = posReadSeq.length();


    // std::cout << "READ: " << readName << std::endl << std::flush; // TEMP


    // At this point, the kmerPositions object should have only the reference sequences.
    std::vector<std::string> referenceNames = kmerPositions->getAllNames();

    // Add both the forward and reverse read sequences to the KmerPositions object.
    kmerPositions->addPositions(posReadName, posReadSeq);
    kmerPositions->addPositions(negReadName, negReadSeq);

    // Create a CommonKmerSet for the read (both forward and reverse complement) and every reference.
    std::vector<CommonKmerSet *> commonKmerSets;
    float maxScoreAllSets = 0.0f;
    CommonKmerSet * commonKmerSet;
    for (size_t i = 0; i < referenceNames.size(); ++i) {
        std::string refName = referenceNames[i];
        int refLength = kmerPositions->getLength(refName);

        // Forward read sequence
        commonKmerSet = new CommonKmerSet(posReadName, refName, readLength, refLength, expectedSlope, kmerPositions);
        if (commonKmerSet->m_maxScore < MINIMUM_MAX_SCORE)
            delete commonKmerSet;
        else {
            commonKmerSets.push_back(commonKmerSet);
            maxScoreAllSets = std::max(maxScoreAllSets, commonKmerSet->m_maxScore);
        }

        // Reverse read sequence
        commonKmerSet = new CommonKmerSet(negReadName, refName, readLength, refLength, expectedSlope, kmerPositions);
        if (commonKmerSet->m_maxScore < MINIMUM_MAX_SCORE)
            delete commonKmerSet;
        else {
            commonKmerSets.push_back(commonKmerSet);
            maxScoreAllSets = std::max(maxScoreAllSets, commonKmerSet->m_maxScore);
        }
    }

    // std::cout << "  MAX SCORE: " << maxScoreAllSets << std::endl << std::flush; // TEMP

    // Sort the common k-mer sets by their max score so high-scoring sets are used first.
    std::sort(commonKmerSets.begin(), commonKmerSets.end(), [](const CommonKmerSet * a, const CommonKmerSet * b) {
        return a->m_maxScore > b->m_maxScore;   
    });

    // Now for the alignments! We first try at sensitivity level 1.
    std::vector<SemiGlobalAlignment *> alignments = semiGlobalAlignmentOneLevel(commonKmerSets, kmerPositions,
                                                                                verbosity, output, expectedSlope,
                                                                                matchScore, mismatchScore,
                                                                                gapOpenScore, gapExtensionScore,
                                                                                1, maxScoreAllSets);

    // If none of the read aligned well, then we give up (assume that the read is rubbish and
    // doesn't) align anywhere. If all of the read aligned well, then we've succeeded and are
    // finished. But if some intermediate fraction of the read aligned well, it seems like the read
    // isn't rubbish but we failed to find all of its alignments. In this case we try again using
    // more sensitive settings.
    double fractionAlignedWell = fractionOfReadAlignedOverThreshold(alignments, lowScoreThreshold);
    if (fractionAlignedWell > 0.0 && fractionAlignedWell < 1.0) {
        std::vector<SemiGlobalAlignment *> l2Alignments = semiGlobalAlignmentOneLevel(commonKmerSets, kmerPositions,
                                                                                      verbosity, output, expectedSlope,
                                                                                      matchScore, mismatchScore,
                                                                                      gapOpenScore, gapExtensionScore,
                                                                                      2, maxScoreAllSets);
        alignments.insert(alignments.end(), l2Alignments.begin(), l2Alignments.end());
    }

    // Clean up.
    kmerPositions->deletePositions(posReadName);
    kmerPositions->deletePositions(negReadName);
    for (size_t i = 0; i < commonKmerSets.size(); ++i)
        delete commonKmerSets[i];

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




std::vector<SemiGlobalAlignment *> semiGlobalAlignmentOneLevel(std::vector<CommonKmerSet *> & commonKmerSets,
                                                               KmerPositions * kmerPositions,
                                                               int verbosity, std::string & output, float expectedSlope,
                                                               int matchScore, int mismatchScore,
                                                               int gapOpenScore, int gapExtensionScore,
                                                               int sensitivityLevel, float maxScoreAllSets) {
    if (verbosity > 2)
        output += "Seqan alignments at sensitivity level " + std::to_string(sensitivityLevel) + ":\n";

    // Set the algorithm settings using the sentitivity level.
    double lowScoreThreshold, highScoreThreshold, minAlignmentLength;
    int minPointCount;
    if (sensitivityLevel == 1) {
        lowScoreThreshold = LOW_SCORE_THRESHOLD_LEVEL_1;
        highScoreThreshold = HIGH_SCORE_THRESHOLD_LEVEL_1;
        minAlignmentLength = MIN_ALIGNMENT_LENGTH_LEVEL_1;
        minPointCount = MIN_POINT_COUNT_LEVEL_1;
    }
    else { // sensitivityLevel == 2
        lowScoreThreshold = LOW_SCORE_THRESHOLD_LEVEL_2;
        highScoreThreshold = HIGH_SCORE_THRESHOLD_LEVEL_2;
        minAlignmentLength = MIN_ALIGNMENT_LENGTH_LEVEL_2;
        minPointCount = MIN_POINT_COUNT_LEVEL_2;
    }

    // The low and high score thresholds are initially expressed as a fraction of the max score.
    // Now turn them into absolute values.
    lowScoreThreshold *= maxScoreAllSets;
    highScoreThreshold *= maxScoreAllSets;

    // std::cout << "  LEVEL: " << sensitivityLevel << std::endl << std::flush; // TEMP
    // std::cout << "  lowScoreThreshold: " << lowScoreThreshold << std::endl << std::flush; // TEMP
    // std::cout << "  highScoreThreshold: " << highScoreThreshold << std::endl << std::flush; // TEMP

    Score<int, Simple> scoringScheme(matchScore, mismatchScore, gapExtensionScore, gapOpenScore);

    // Go through the common k-mer sets and perform line-finding and then aligning.
    std::vector<SemiGlobalAlignment *> alignments;
    for (size_t i = 0; i < commonKmerSets.size(); ++i) {
        CommonKmerSet * commonKmerSet = commonKmerSets[i];

        // If a common k-mer set's max score is below the high threshold, then we know there won't
        // be any alignment lines, so don't bother continuing.
        if (commonKmerSet->m_maxScore < highScoreThreshold)
            continue;

        std::string readName = commonKmerSet->m_readName;
        std::string refName = commonKmerSet->m_refName;
        std::string * readSeq = kmerPositions->getSequence(readName);
        std::string * refSeq = kmerPositions->getSequence(refName);
        int readLength = readSeq->length();
        int refLength = refSeq->length();

        // std::cout << "  REF: " << refName << std::endl << std::flush; // TEMP

        std::vector<AlignmentLine *> alignmentLines = findAlignmentLines(commonKmerSet, readLength, refLength, expectedSlope,
                                                                         verbosity, output, 
                                                                         lowScoreThreshold, highScoreThreshold);
        if (alignmentLines.size() == 0)
            continue;

        for (size_t j = 0; j < alignmentLines.size(); ++j) {
            AlignmentLine * line = alignmentLines[j];

            // std::cout << "    LINE " << j << ": points = " << line->m_linePoints.size() << std::endl << std::flush; // TEMP

            bool seedChainSuccess = line->buildSeedChain(minPointCount, minAlignmentLength);
            if (seedChainSuccess) {

                // std::cout << "      slope = " << line->m_slope << ", intercept = " << line->m_intercept << std::endl << std::flush; // TEMP

                SemiGlobalAlignment * alignment = semiGlobalAlignmentOneLine(readName, refName, readSeq, refSeq,
                                                                             line, verbosity, output, scoringScheme);

                // std::cout << "      alignment: " << alignment->getShortDisplayString() << std::endl << std::flush; // TEMP

                if (alignment != 0)
                    alignments.push_back(alignment);
            }
        }

        // Clean up.
        for (size_t j = 0; j < alignmentLines.size(); ++j)
            delete alignmentLines[j];
    }

    return alignments;
}



 // Runs an alignment using Seqan between one read and one reference along one line.
 // It starts with a smallish band size (fast) and works up to larger ones to see if they improve
 // the alignment.
SemiGlobalAlignment * semiGlobalAlignmentOneLine(std::string & readName, std::string & refName,
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
    SemiGlobalAlignment * bestAlignment = 0;
    double bestAlignmentScore = std::numeric_limits<double>::min();

    // We perform the alignment with increasing band sizes until the score stops improving or we
    // reach the max band size.
    while (true) {
        SemiGlobalAlignment * alignment = semiGlobalAlignmentOneLineOneBand(readName, refName,
                                                                            readSeqSeqan, readLength,
                                                                            refSeqSeqan, trimmedRefLength,
                                                                            line, bandSize, verbosity,
                                                                            output, scoringScheme);
        if (alignment != 0) {
            double alignmentScore = alignment->m_scaledScore;
            if (alignmentScore <= bestAlignmentScore) {
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

    if (bestAlignment != 0) {
        if (verbosity > 2)
            output += "  " + bestAlignment->getShortDisplayString() + ", band size = " + std::to_string(bandSize) + "\n";
        if (verbosity > 3)
            output += "    " + bestAlignment->m_cigar + "\n";
        bestAlignment->m_milliseconds = getTime() - startTime;
    }
    return bestAlignment;
}






// This function, given a line, will search for semi-global alignments around that line. The
// bandSize parameter specifies how far of an area around the line is searched.
SemiGlobalAlignment * semiGlobalAlignmentOneLineOneBand(std::string & readName, std::string & refName,
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

    SemiGlobalAlignment * sgAlignment;
    try {
        bandedChainAlignment(alignment, line->m_bridgedSeedChain, scoringScheme, alignConfig,
                             bandSize);
        sgAlignment = new SemiGlobalAlignment(alignment, readName, refName, readLen, refLen,
                                              line->m_trimmedRefStart, startTime, false, false, scoringScheme);
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

    SemiGlobalAlignment startAlignment(alignment, readName, refName, length(read), length(ref),
                                       0, startTime, false, true, scoringScheme);
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

    SemiGlobalAlignment endAlignment(alignment, readName, refName, length(read), length(ref),
                                     0, startTime, true, false, scoringScheme);
    return cppStringToCString(endAlignment.getFullString());
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


    // def get_fraction_aligned(self):
    //     '''
    //     This function returns the fraction of the read which is covered by any of the read's
    //     alignments.
    //     '''
    //     read_ranges = [x.read_start_end_positive_strand() \
    //                    for x in self.alignments]
    //     read_ranges = simplify_ranges(read_ranges)
    //     aligned_length = sum([x[1] - x[0] for x in read_ranges])
    //     return aligned_length / len(self.sequence)


// Returns true if some parts of the read fail to meet the score threshold
double fractionOfReadAlignedOverThreshold(std::vector<SemiGlobalAlignment *> & alignments, double scoreThreshold) {
    std::vector<SemiGlobalAlignment *> goodScoreAlignments;
    for (size_t i = 0; i < alignments.size(); ++i) {
        if (alignments[i]->m_scaledScore >= scoreThreshold)
            goodScoreAlignments.push_back(alignments[i]);
    }
    return fractionOfReadAligned(goodScoreAlignments);
}

double fractionOfReadAligned(std::vector<SemiGlobalAlignment *> & alignments) {
    if (alignments.size() == 0)
        return true;
    std::vector<std::pair<int, int> > ranges;
    for (size_t i = 0; i < alignments.size(); ++i) {
        SemiGlobalAlignment * alignment = alignments[i];
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

