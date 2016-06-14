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

char * multipleSequenceAlignment(char * sequences[], char * qualities[], int sequenceCount, int pieceSize,
                                 int matchScore, int mismatchScore, int gapOpenScore, int gapExtensionScore) {
    // Convert the input C-string arrays to C++ string vectors.
    std::vector<std::string> ungappedSequences;
    ungappedSequences.reserve(sequenceCount);
    for (int i = 0; i < sequenceCount; ++i)
        ungappedSequences.push_back(std::string(sequences[i]));
    std::vector<std::string> ungappedQualities;
    ungappedQualities.reserve(sequenceCount);
    for (int i = 0; i < sequenceCount; ++i)
        ungappedQualities.push_back(std::string(qualities[i]));

    // These vectors will hold the aligned sequences and qualities.
    std::vector<std::string> gappedSequences;
    gappedSequences.reserve(sequenceCount);
    std::vector<std::string> gappedQualities;
    gappedQualities.reserve(sequenceCount);

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

    int totalSeqLength = 0;
    for (int i = 0; i < sequenceCount; ++i)
        totalSeqLength += ungappedSequences[i].length();
    int meanSeqLength = int(0.5 + double(totalSeqLength) / sequenceCount);

    // If the sequences are long, we will do the alignment in pieces. This will make it run much
    // faster, and since we expect the alignment to follow the diagonal pretty well, this shouldn't
    // affect the alignment's quality.
    if (meanSeqLength > pieceSize) {
        gappedSequences.resize(sequenceCount, "");
        std::vector<std::string> ungappedSequencePieces;
        ungappedSequencePieces.reserve(sequenceCount);
        std::vector<int> ungappedSequencePieceIndices;
        ungappedSequencePieceIndices.reserve(sequenceCount);
        std::vector<std::string> gappedSequencePieces;
        gappedSequencePieces.reserve(sequenceCount);
        std::vector<int> startPositions;
        startPositions.resize(sequenceCount, 0);

        bool lastPiece = false;
        while (!lastPiece) {
            // std::cout << "\n\n\n\n"; // TEMP

            ungappedSequencePieces.clear();
            ungappedSequencePieceIndices.clear();
            lastPiece = true;
            for (int i = 0; i < sequenceCount; ++i) {
                std::string ungappedSequencePiece = ungappedSequences[i].substr(startPositions[i], pieceSize);
                if (!ungappedSequencePiece.empty()) {
                    ungappedSequencePieces.push_back(ungappedSequencePiece);
                    ungappedSequencePieceIndices.push_back(i);
                    // std::cout << ungappedSequencePieces.back() << "\n"; // TEMP
                }
                if (int(ungappedSequencePiece.length()) == pieceSize)
                    lastPiece = false;
            }

            Align<Dna5String> align;
            resize(rows(align), ungappedSequencePieces.size());
            for (size_t i = 0; i < ungappedSequencePieces.size(); ++i)
                assignSource(row(align, i), ungappedSequencePieces[i]);
            globalMsaAlignment(align, scoringScheme);

            // std::cout << "\n" << align << "\n"; // TEMP

            gappedSequencePieces.clear();
            gappedSequencePieces.resize(sequenceCount, "");
            int longestGappedSequencePiece = 0;
            for (size_t i = 0; i < ungappedSequencePieces.size(); ++i) {
                std::ostringstream stream;
                stream << row(align, i);
                std::string gappedSequencePiece = stream.str();
                gappedSequencePieces[ungappedSequencePieceIndices[i]] = gappedSequencePiece;
                longestGappedSequencePiece = std::max(longestGappedSequencePiece, int(gappedSequencePiece.length()));
            }
            for (int i = 0; i < sequenceCount; ++i) {
                if (gappedSequencePieces[i].empty())
                    gappedSequencePieces[i] = std::string(longestGappedSequencePiece, '-');
            }

            // If this was the last piece, then we keep the whole alignment.
            if (lastPiece) {
                for (int i = 0; i < sequenceCount; ++i)
                    gappedSequences[i] += gappedSequencePieces[i];
            }

            // If the wasn't the last piece, we keep the first half of the alignment and advance the
            // start positions as appropriate.
            else {
                int pieceAlignmentLength = gappedSequencePieces[0].length();
                int halfPieceAlignmentLength = pieceAlignmentLength / 2;
                for (int i = 0; i < sequenceCount; ++i) {
                    std::string keptGappedSequence = gappedSequencePieces[i].substr(0, halfPieceAlignmentLength);
                    gappedSequences[i] += keptGappedSequence;
                    startPositions[i] += (halfPieceAlignmentLength - std::count(keptGappedSequence.begin(), keptGappedSequence.end(), '-'));
                }
            }
        }
    }

    // If the sequences are short enough, we can just do the alignment all at once.
    else {
        Align<Dna5String> align;
        resize(rows(align), sequenceCount);
        for (int i = 0; i < sequenceCount; ++i)
            assignSource(row(align, i), ungappedSequences[i]);
        globalMsaAlignment(align, scoringScheme);

        // std::cout << "\n" << align << "\n"; // TEMP

        for (int i = 0; i < sequenceCount; ++i) {
            std::ostringstream stream;
            stream << row(align, i);
            gappedSequences.push_back(stream.str());
        }
    }

    // Add gaps to the quality scores so they match up with the bases.
    int alignmentLength = gappedSequences[0].length();
    for (int i = 0; i < sequenceCount; ++i) {
        std::cout << gappedSequences[i] << "\n"; // TEMP
        std::string gappedQuality;
        gappedQuality.resize(gappedSequences[i].length(), ' ');
        int pos = 0;
        for (int j = 0; j < alignmentLength; ++j) {
            if (gappedSequences[i][j] != '-')
                gappedQuality[j] = ungappedQualities[i][pos++];
        }
        gappedQualities.push_back(gappedQuality);
    }

    // For each gapped sequence, get the position of the first and last non-N base. This is useful
    // for sequences which don't span the whole alignment but are just at the start or end.
    std::vector<int> firstNonN;
    firstNonN.reserve(sequenceCount);
    std::vector<int> lastNonN;
    lastNonN.reserve(sequenceCount);
    for (int i = 0; i < sequenceCount; ++i) {
        size_t firstNPos = gappedSequences[i].find_first_of('N');
        if (firstNPos == std::string::npos) {
            firstNonN.push_back(0);
            lastNonN.push_back(gappedSequences[i].length() - 1);
        }
        else {
            size_t firstNonNPos = gappedSequences[i].find_first_of("ACGTacgt");
            size_t lastNonNPos = gappedSequences[i].find_last_of("ACGTacgt");
            if (firstNonNPos != std::string::npos) {
                firstNonN.push_back(int(firstNonNPos));
                lastNonN.push_back(int(lastNonNPos));
            }
            else {
                firstNonN.push_back(0);
                lastNonN.push_back(gappedSequences[i].length() - 1);
            }
        }
        // std::cout << firstNonN.back() << "  " << lastNonN.back() << "\n"; // TEMP
    }

    // Build a consensus sequence. Sequences are ignored before their first non-N base was seen
    // (for end-only sequences) and after their last non-N base was seen (for start-only
    // sequences).
    std::string consensus;
    std::string gappedConsensus;
    for (int i = 0; i < alignmentLength; ++i) {
        std::vector<char> bases;
        std::vector<char> qualities;
        bases.reserve(sequenceCount);
        qualities.reserve(sequenceCount);

        // std::cout << "Position:  " << i << "\n"; // TEMP

        for (int j = 0; j < sequenceCount; ++j) {
            if (i < firstNonN[j] || i > lastNonN[j])
                continue;
            char base = toupper(gappedSequences[j][i]);
            char quality = gappedQualities[j][i];
            if (base != 'N') {
                bases.push_back(base);
                qualities.push_back(quality);
            }
        }
        if (bases.size() > 0) {
            char mostCommonBase = getMostCommonBase(bases, qualities);
            // std::cout << "Call:      " << mostCommonBase << "\n\n"; // TEMP
            if (mostCommonBase != '-')
                consensus.push_back(mostCommonBase);
            gappedConsensus.push_back(mostCommonBase);
        }
    }
    std::cout << "\n" << gappedConsensus << "\n"; // TEMP

    // Score each sequence against the consensus.
    std::vector<double> scaledScores;
    scaledScores.reserve(sequenceCount);
    for (int i = 0; i < sequenceCount; ++i) {
        int rawScore = scoreAlignment(gappedConsensus, gappedSequences[i], firstNonN[i], lastNonN[i],
                                      matchScore, mismatchScore, gapOpenScore, gapExtensionScore);
        int seqLength = ungappedSequences[i].length();
        int seqLengthNoN = seqLength;
        for (int j = 0; j < seqLength; ++j) {
            if (ungappedSequences[i][j] == 'N')
                --seqLengthNoN;
        }
        int perfectScore = matchScore * seqLengthNoN;
        int worstScore = mismatchScore * seqLengthNoN;
        double scaledScore;
        if (perfectScore > 0)
            scaledScore = 100.0 * double(rawScore - worstScore) / double(perfectScore - worstScore);
        else
            scaledScore = 0.0;
        scaledScores.push_back(scaledScore);
        // std::cout << rawScore << "   " << scaledScore << "\n"; // TEMP
    }

    std::string returnString = consensus;
    returnString += ';' + std::to_string(scaledScores[0]);
    for (int i = 1; i < sequenceCount; ++i)
        returnString += ',' + std::to_string(scaledScores[i]);
    return cppStringToCString(returnString);
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


// Takes two sequences (with gaps) and produces a raw alignment score.
int scoreAlignment(std::string & seq1, std::string & seq2, int startPos, int endPos,
                   int matchScore, int mismatchScore, int gapOpenScore, int gapExtensionScore) {

    // std::cout << "\n"; // TEMP

    bool insertionInProgress = false;
    bool deletionInProgress = false;
    int score = 0;
    for (int i = startPos; i <= endPos; ++i) {
        char base1 = seq1[i];
        char base2 = seq2[i];
        if (base1 == '-' && base2 == '-') {
            // std::cout << " "; // TEMP
            continue;
        }
        if (base1 == 'N' || base2 == 'N') {
            // std::cout << " "; // TEMP
            continue;
        }
        if (base1 == '-') {
            if (insertionInProgress)
                score += gapExtensionScore;
            else
                score += gapOpenScore;
            // std::cout << "I"; // TEMP
            insertionInProgress = true;
            deletionInProgress = false;
        }
        else if (base2 == '-') {
            if (deletionInProgress)
                score += gapExtensionScore;
            else
                score += gapOpenScore;
            // std::cout << "D"; // TEMP
            insertionInProgress = false;
            deletionInProgress = true;
        }
        else { // match or mismatch
            if (base1 == base2) {
                score += matchScore;
                // std::cout << "M"; // TEMP
            }
            else { // base1 != base2
                score += mismatchScore;
                // std::cout << "X"; // TEMP
            }
            insertionInProgress = false;
            deletionInProgress = false;
        }
    }

    // std::cout << "\n"; // TEMP

    return score;
}

