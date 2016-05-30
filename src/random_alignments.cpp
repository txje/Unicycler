#include "random_alignments.h"

#include <seqan/align.h>
#include <iostream>
#include <vector>
#include <limits>
#include <algorithm>
#include <utility>
#include <random>
#include "seqan_align.h"



// This function runs a bunch of alignments between random sequences to get a mean and std dev of
// the scaled scores. It return them in a C string (for Python).
char * getRandomSequenceAlignmentScores(int seqLength, int n,
                                        int matchScore, int mismatchScore, int gapOpenScore, int gapExtensionScore) {
    std::vector<double> scores;

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<int> dist(0, 3);

    for (int i = 0; i < n; ++i) {
        std::string s1 = getRandomSequence(seqLength, gen, dist);
        std::string s2 = getRandomSequence(seqLength, gen, dist);
        SemiGlobalAlignment * alignment = fullyGlobalAlignment(s1, s2, matchScore, mismatchScore, gapOpenScore, gapExtensionScore);

        if (alignment != 0) {
            scores.push_back(alignment->m_scaledScore);
            delete alignment;
        }
    }

    double mean, stdev;
    getMeanAndStDev(scores, mean, stdev);
    return cppStringToCString(std::to_string(mean) + "," + std::to_string(stdev));
}

// This function returns lots of information about random global alignments.
char * getRandomSequenceAlignmentErrorRates(int seqLength, int n,
                                            int matchScore, int mismatchScore, int gapOpenScore, int gapExtensionScore) {
    std::vector<double> matchesOverAlignmentLength;
    std::vector<double> errorsOverAlignmentLength;
    std::vector<double> matchesOverSeqLength;
    std::vector<double> errorsOverSeqLength;
    std::vector<double> errorsOverSeqLengthCountAllInsertionsAsOne;
    std::vector<double> scores;
    std::vector<double> alignmentLengths;

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<int> dist(0, 3);

    for (int i = 0; i < n; ++i) {
        int totalMatches = 0;
        int totalErrors = 0;
        int totalErrorsCountAllInsertionsAsOne = 0;

        std::string s1 = getRandomSequence(seqLength, gen, dist);
        std::string s2 = getRandomSequence(seqLength, gen, dist);

        Dna5String sequenceH(s1);
        Dna5String sequenceV(s2);

        Align<Dna5String, ArrayGaps> alignment;
        resize(rows(alignment), 2);
        assignSource(row(alignment, 0), sequenceH);
        assignSource(row(alignment, 1), sequenceV);

        Score<int, Simple> scoringScheme(matchScore, mismatchScore, gapExtensionScore, gapOpenScore);

        AlignConfig<false, false, false, false> alignConfig;
        int score = globalAlignment(alignment, scoringScheme, alignConfig);

        // Extract the alignment sequences into C++ strings for constant time random access.
        std::ostringstream stream1;
        stream1 << row(alignment, 0);
        std::string s1Alignment =  stream1.str();
        std::ostringstream stream2;
        stream2 << row(alignment, 1);
        std::string s2Alignment =  stream2.str();
        int alignmentLength = s1Alignment.size();

        int insertionLength = 0;
        for (size_t i = 0; i < s1Alignment.size(); ++i)
        {
            if (s1Alignment[i] != '-' && s2Alignment[i] == '-')
                ++insertionLength;
            else
                insertionLength = 0;

            if (s1Alignment[i] == '-' || s2Alignment[i] == '-' || s1Alignment[i] != s2Alignment[i]) {
                ++totalErrors;
                if (insertionLength < 2)
                    ++totalErrorsCountAllInsertionsAsOne;
            }
            else // match
                ++totalMatches;
        }
        matchesOverAlignmentLength.push_back(double(totalMatches) / double(alignmentLength));
        errorsOverAlignmentLength.push_back(double(totalErrors) / double(alignmentLength));
        matchesOverSeqLength.push_back(double(totalMatches) / double(seqLength));
        errorsOverSeqLength.push_back(double(totalErrors) / double(seqLength));
        errorsOverSeqLengthCountAllInsertionsAsOne.push_back(double(totalErrorsCountAllInsertionsAsOne) / double(seqLength));
        scores.push_back(double(score));
        alignmentLengths.push_back(double(s1Alignment.size()));
    }

    std::string returnString;
    double mean, stdev;

    returnString += "identity (using alignment length) mean\t";
    returnString += "identity (using alignment length) std dev\t";
    returnString += "errors (using alignment length) mean\t";
    returnString += "errors (using alignment length) std dev\t";
    returnString += "identity (using sequence length) mean\t";
    returnString += "identity (using sequence length) std dev\t";
    returnString += "errors (using sequence length) mean\t";
    returnString += "errors (using sequence length) std dev\t";
    returnString += "errors (using sequence length, count all insertions as one) mean\t";
    returnString += "errors (using sequence length, count all insertions as one) std dev\t";
    returnString += "score mean\t";
    returnString += "score std dev\t";
    returnString += "alignment length mean\t";
    returnString += "alignment length std dev\n";

    getMeanAndStDev(matchesOverAlignmentLength, mean, stdev);
    returnString += std::to_string(mean) + "\t";
    returnString += std::to_string(stdev) + "\t";
    getMeanAndStDev(errorsOverAlignmentLength, mean, stdev);
    returnString += std::to_string(mean) + "\t";
    returnString += std::to_string(stdev) + "\t";
    getMeanAndStDev(matchesOverSeqLength, mean, stdev);
    returnString += std::to_string(mean) + "\t";
    returnString += std::to_string(stdev) + "\t";
    getMeanAndStDev(errorsOverSeqLength, mean, stdev);
    returnString += std::to_string(mean) + "\t";
    returnString += std::to_string(stdev) + "\t";
    getMeanAndStDev(errorsOverSeqLengthCountAllInsertionsAsOne, mean, stdev);
    returnString += std::to_string(mean) + "\t";
    returnString += std::to_string(stdev) + "\t";
    getMeanAndStDev(scores, mean, stdev);
    returnString += std::to_string(mean) + "\t";
    returnString += std::to_string(stdev) + "\t";
    getMeanAndStDev(alignmentLengths, mean, stdev);
    returnString += std::to_string(mean) + "\t";
    returnString += std::to_string(stdev) + "\n";

    return cppStringToCString(returnString);
}

std::string getRandomSequence(int seqLength, std::mt19937 & gen, std::uniform_int_distribution<int> & dist) {
    std::string seq;
    seq.reserve(seqLength);
    for (int i = 0 ; i < seqLength; ++i)
        seq += getRandomBase(gen, dist);
    return seq;
}

char getRandomBase(std::mt19937 & gen, std::uniform_int_distribution<int> & dist) {
    int baseNum = dist(gen);
    if (baseNum == 0)
        return 'A';
    else if (baseNum == 1)
        return 'C';
    else if (baseNum == 2)
        return 'G';
    else // baseNum == 3
        return 'T';
}

// This function runs a global alignment (slow) between two sequences.
SemiGlobalAlignment * fullyGlobalAlignment(std::string s1, std::string s2,
                                           int matchScore, int mismatchScore, int gapOpenScore, int gapExtensionScore) {
    long long startTime = getTime();

    Dna5String sequenceH(s1);
    Dna5String sequenceV(s2);

    Align<Dna5String, ArrayGaps> alignment;
    resize(rows(alignment), 2);
    assignSource(row(alignment, 0), sequenceH);
    assignSource(row(alignment, 1), sequenceV);
    Score<int, Simple> scoringScheme(matchScore, mismatchScore, gapExtensionScore, gapOpenScore);

    AlignConfig<false, false, false, false> alignConfig;
    globalAlignment(alignment, scoringScheme, alignConfig);

    std::string s1Name = "s1";
    std::string s2Name = "s2";

    return new SemiGlobalAlignment(alignment, s1Name, s2Name, s1.length(), s2.length(),
                                   0, startTime, 0, true, true, scoringScheme);
}



void getMeanAndStDev(std::vector<double> & v, double & mean, double & stdDev) {
    mean = 0.0;
    stdDev = 0.0;
    int count = v.size();
    if (count < 1)
        return;
    double devSum = 0.0;
    for (int i = 0; i < count; ++i)
        mean += v[i];
    mean /= count;
    for (int i = 0; i < count; ++i) {
        double dev = v[i] - mean;
        devSum += dev * dev;
    }
    stdDev = sqrt(devSum / v.size());
}

