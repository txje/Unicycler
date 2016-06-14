#include "consensus_align.h"

#include <iostream>
#include <map>
#include "seqan_align.h"
#include <algorithm>
#include <cmath>


// template <typename TScoreValue, typename TSequenceValue, typename TSpec>
// void showScoringMatrix(Score<TScoreValue, ScoreMatrix<TSequenceValue, TSpec> > const & scoringScheme)
// {
//     // Print top row.
//     for (unsigned i = 0; i < ValueSize<TSequenceValue>::VALUE; ++i)
//         std::cout << "\t" << TSequenceValue(i);
//     std::cout << std::endl;
//     // Print each row.
//     for (unsigned i = 0; i < ValueSize<TSequenceValue>::VALUE; ++i)
//     {
//         std::cout << TSequenceValue(i);
//         for (unsigned j = 0; j < ValueSize<TSequenceValue>::VALUE; ++j)
//         {
//             std::cout << "\t" << score(scoringScheme, TSequenceValue(i), TSequenceValue(j));
//         }
//         std::cout << std::endl;
//     }
// }

char * multipleSequenceAlignment(char * sequences[], char * qualities[], int sequenceCount,
                                 int matchScore, int mismatchScore, int gapOpenScore, int gapExtensionScore) {

    // We want to use a scoring scheme which is almost like the simple scoring scheme we use for
    // read alignment, but with the addition of free Ns in the middle.
    Score<int, ScoreMatrix<Dna5, Default> > scoringScheme(gapExtensionScore, gapOpenScore);
    for (int i = 0; i < 5; ++i) {
        for (int j = 0; j < 5; ++j) {
            int score;
            if (i == 4 || j == 4) // either base is N
                score = 0;
            else if (i == j)
                score = matchScore;
            else
                score = mismatchScore;
            setScore(scoringScheme, Dna5(i), Dna5(j), score);
        }
    }

    // Perform the multiple sequence alignment.
    Align<Dna5String> align;
    resize(rows(align), sequenceCount);
    for (int i = 0; i < sequenceCount; ++i)
        assignSource(row(align, i), sequences[i]);
    globalMsaAlignment(align, scoringScheme);

    // std::cout << "\n" << align << "\n"; // TEMP

    // Extract the alignment sequences into C++ strings.
    std::vector<std::string> gappedSequences;
    for (int i = 0; i < sequenceCount; ++i) {
        std::ostringstream stream;
        stream << row(align, i);
        gappedSequences.push_back(stream.str());
    }
    int alignmentLength = gappedSequences[0].length();

    // Add gaps to the quality scores so they match up with the bases.
    std::vector<std::string> gappedQualities;
    for (int i = 0; i < sequenceCount; ++i) {
        std::string ungappedQuality(qualities[i]);
        std::string gappedQuality;
        gappedQuality.resize(gappedSequences[i].length(), ' ');
        int pos = 0;
        for (int j = 0; j < alignmentLength; ++j) {
            if (gappedSequences[i][j] != '-')
                gappedQuality[j] = ungappedQuality[pos++];
        }
        gappedQualities.push_back(gappedQuality);
    }

    // Build a consensus sequence. Sequences are ignored before their first non-N base was seen
    // (for end-only sequences) and after their last non-N base was seen (for start-only
    // sequences).
    std::vector<bool> encounteredNonNBase;
    encounteredNonNBase.resize(sequenceCount, false);
    std::vector<bool> lastBaseWasN;
    lastBaseWasN.resize(sequenceCount, false);

    std::string consensus;
    // std::string gappedConsensus; // TEMP
    for (int i = 0; i < alignmentLength; ++i) {
        std::vector<char> bases;
        std::vector<char> qualities;
        bases.reserve(sequenceCount);
        qualities.reserve(sequenceCount);

        // std::cout << "Position:  " << i << "\n"; // TEMP

        for (int j = 0; j < sequenceCount; ++j) {
            char base = toupper(gappedSequences[j][i]);
            char quality = gappedQualities[j][i];

            if (base == 'N')
                lastBaseWasN[j] = true;
            else if (base != '-') // is A, C, G or T
                lastBaseWasN[j] = false;
                encounteredNonNBase[j] = true;

            if (base != 'N' and encounteredNonNBase[j] and !lastBaseWasN[j]) {
                bases.push_back(base);
                qualities.push_back(quality);
            }
        }
        if (bases.size() > 0) {
            char mostCommonBase = getMostCommonBase(bases, qualities);
            // std::cout << "Call:      " << mostCommonBase << "\n\n"; // TEMP
            if (mostCommonBase != '-')
                consensus.push_back(mostCommonBase);
            // gappedConsensus.push_back(mostCommonBase); // TEMP
        }
    }

    // std::cout << "\n" << gappedConsensus << "\n"; // TEMP

    return cppStringToCString(consensus);
}

char getMostCommonBase(std::vector<char> & bases, std::vector<char> & qualities) {
    std::string baseValues = "ACGT-";

    // std::cout << "Bases:     "; // TEMP
    // for (size_t i = 0; i < bases.size(); ++i) // TEMP
    //     std::cout << bases[i]; // TEMP
    // std::cout << "\n"; // TEMP

    // std::cout << "Qualities: "; // TEMP
    // for (size_t i = 0; i < qualities.size(); ++i) // TEMP
    //     std::cout << qualities[i]; // TEMP
    // std::cout << "\n"; // TEMP

    // Tally the count for each base.
    std::map<char, int> baseCounts;
    for (int i = 0; i < 5; ++i)
        baseCounts[baseValues[i]] = 0;
    for (size_t i = 0; i < bases.size(); ++i)
        baseCounts[bases[i]]++;

    // std::cout << "Counts:    "; // TEMP
    // std::cout << "A: " << baseCounts['A'] << "  "; // TEMP
    // std::cout << "C: " << baseCounts['C'] << "  "; // TEMP
    // std::cout << "G: " << baseCounts['G'] << "  "; // TEMP
    // std::cout << "T: " << baseCounts['T'] << "  "; // TEMP
    // std::cout << "-: " << baseCounts['-'] << "\n"; // TEMP

    // If only one base (or a gap) is the most common, return that.
    int largestCount = 0;
    for (int i = 0; i < 5; ++i)
        largestCount = std::max(largestCount, baseCounts[baseValues[i]]);
    std::vector<char> mostCommonBases;
    for (int i = 0; i < 5; ++i) {
        if (baseCounts[baseValues[i]] == largestCount)
            mostCommonBases.push_back(baseValues[i]);
    }
    if (mostCommonBases.size() == 1)
        return mostCommonBases[0];

    // If the code got here, then there's a tie. Check to see if one is a gap, as that should
    // always lose to an actual base.
    if (mostCommonBases.size() == 2 and mostCommonBases[1] == '-')
        return mostCommonBases[0];

    // If the code got here, then there's a tie between two or more bases, so we'll need to use the
    // qualities to decide. The base with the largest Phred sum wins. Since we're just summing the
    // Phred scores, there's no need to worry about +33 vs +64 offset.
    std::map<char, int> phredSums;
    for (int i = 0; i < 4; ++i)
        phredSums[baseValues[i]] = 0;
    for (size_t i = 0; i < bases.size(); ++i)
        phredSums[bases[i]] += int(qualities[i]);
    int largestPhredSum = 0;
    for (size_t i = 0; i < mostCommonBases.size(); ++i)
        largestPhredSum = std::max(largestPhredSum, phredSums[mostCommonBases[i]]);

    // While looping through the bases in their original order, return the first one which has the
    // largest Phred sum and was one of the most common bases. This ensure that a tie (when there
    // are two equally common bases with the same Phred sum) it is broken by the original order in
    // the bases vector, which is arbitrary but at least consistent and not biased towards a
    // particular base.
    for (size_t i = 0; i < bases.size(); ++i) {
        char base = bases[i];
        if (base != '-' && phredSums[base] == largestPhredSum &&
            std::find(mostCommonBases.begin(), mostCommonBases.end(), base) != mostCommonBases.end())
            return base;
    }

    // The code should never get here, as the most common base with the biggest Phred sum should
    // found above.
    return '-';
}

