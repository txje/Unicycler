#include <stdio.h>
#include <chrono>
#include <seqan/basic.h>
#include <seqan/align.h>
#include <seqan/sequence.h>
#include <seqan/stream.h>
#include <seqan/seeds.h>
#include <string>
#include <set>
#include <tuple>
#include <vector>
#include <map>
#include <algorithm>
#include <cmath>

using namespace seqan;
using namespace std::chrono;

extern "C" {

typedef Dna5String TSequence;
typedef Align<TSequence, ArrayGaps> TAlign;
typedef Row<TAlign>::Type TRow;
typedef Seed<Simple> TSeed;
typedef SeedSet<TSeed> TSeedSet;
typedef Row<TAlign>::Type TRow;
typedef Iterator<TRow>::Type TRowIterator;

enum CigarType {MATCH, INSERTION, DELETION, CLIP, NOTHING};

class CommonKmer
{
public:
    CommonKmer(std::string sequence, int hPosition, int vPosition, double angle);
    static double getRotationAngle(double slope) {return -atan(slope);}

    std::string m_sequence;
    int m_hPosition;
    int m_vPosition;
    double m_rotatedHPosition;
    double m_rotatedVPosition;
    double m_bandLength; // Length of band at expected slope
    int m_bandCount; // Number of neighbour points in band
    double m_score; // Score calculated from bandCount, bandLength and overall kmer density
};


// These are the functions that will be called by the Python script.
char * bandedSemiGlobalAlignment(char * s1, char * s2, int s1Len, int s2Len, double slope,
                                 double intercept, int kSize, int bandSize, int verbosity,
                                 int matchScore, int mismatchScore, int gapOpenScore,
                                 int gapExtensionScore, char * kmerLocations);
char * findAlignmentLines(char * s1, char * s2, int s1Len, int s2Len, double expectedSlope,
                          int verbosity);
// char * exhaustiveSemiGlobalAlignment(char * s1, char * s2, int s1Len, int s2Len);
void free_c_string(char * p);

// These functions are internal to this C++ code.
void addBridgingSeed(TSeedSet & seedSet, int hStart, int vStart, int hEnd, int vEnd);
std::vector<CommonKmer> getCommonKmers(std::string * s1, std::string * s2, double expectedSlope,
                                       int verbosity, std::string & output);
std::map<std::string, std::vector<int> > getCommonLocations(std::string * s1, std::string * s2, int kSize);
char * turnAlignmentIntoDescriptiveString(Align<Dna5String, ArrayGaps> * alignment, int score,
                                          int matchScore, long long startTime, std::string output);
CigarType getCigarType(char b1, char b2, bool alignmentStarted);
std::string getCigarPart(CigarType type, int length);
char * cppStringToCString(std::string cpp_string);
std::string vectorToString(std::vector<int> * v);
void fixSeedChainToLine(String<TSeed> * seedChain, double bandSize);
void getSeedChainSlopeAndIntercept(String<TSeed> * seedChain, int firstI, int lastI,
                                   double * slope, double * intercept);
void getSlopeAndIntercept(int hStart, int hEnd, int vStart, int vEnd,
                                   double * slope, double * intercept);
double getMedian(std::vector<double> * v);
void printKmerSize(int kmerSize, int locationCount, std::string & output);
double getLineLength(double x, double y, double slope, double xSize, double ySize);
void linearRegression(std::vector<CommonKmer> & pts, double * slope, double * intercept);
void parseKmerLocationsFromString(std::string & str, std::vector<int> & v1, std::vector<int> & v2);
std::string getKmerTable(std::vector<CommonKmer> & commonKmers);
std::string getSeedChainTable(String<TSeed> & seedChain);
void addSeedMerge(TSeedSet & seedSet, TSeed & seed);


// // This function does the full semi-global alignment using the entirety of both sequences. It will
// // be slow but will always find the ideal alignment.
// char * exhaustiveSemiGlobalAlignment(char * s1, char * s2, int s1Len, int s2Len)
// {
//     std::string output = "";
//     long long startTime = duration_cast< milliseconds >(system_clock::now().time_since_epoch()).count();

//     TSequence sequenceH = s1;
//     TSequence sequenceV = s2;

//     Align<Dna5String, ArrayGaps> alignment;
//     resize(rows(alignment), 2);
//     assignSource(row(alignment, 0), sequenceH);
//     assignSource(row(alignment, 1), sequenceV);
//     Score<int, Simple> scoringScheme(matchScore, mismatchScore, gapExtensionScore, gapOpenScore);
//     AlignConfig<true, true, true, true> alignConfig;
//     int score = globalAlignment(alignment, scoringScheme, alignConfig);

//     return turnAlignmentIntoDescriptiveString(&alignment, score, startTime, output);
// }




// This function, given a line, will search for semi-global alignments around that line. The
// bandSize parameter specifies how far of an area around the line is searched. It is generally
// much faster than exhaustiveSemiGlobalAlignment, though it may not find the optimal alignment.
// A lower bandSize is faster with a larger chance of missing the optimal alignment.
char * bandedSemiGlobalAlignment(char * s1, char * s2, int s1Len, int s2Len, double slope,
                                 double intercept, int kSize, int bandSize, int verbosity,
                                 int matchScore, int mismatchScore, int gapOpenScore,
                                 int gapExtensionScore, char * kmerLocations)
{
    std::string output = "";
    long long startTime = duration_cast< milliseconds >(system_clock::now().time_since_epoch()).count();

    // Extreme slope values will not work.
    if (slope > 1.5)
        return cppStringToCString(output + ";Failed: slope too large");
    if (slope < 0.5)
        return cppStringToCString(output + ";Failed: slope too small");

    TSequence sequenceH = s1;
    TSequence sequenceV = s2;

    std::string kmerLocationsStr(kmerLocations);
    std::vector<int> s1KmerLocations;
    std::vector<int> s2KmerLocations;
    parseKmerLocationsFromString(kmerLocationsStr, s1KmerLocations, s2KmerLocations);

    // Build a Seqan seed set using our common k-mers.
    TSeedSet seedSet;
    for (int i = 0; i < s1KmerLocations.size(); ++i)
    {
        TSeed seed(s1KmerLocations[i], s2KmerLocations[i], kSize);
        addSeedMerge(seedSet, seed);
    }

    // We now get a Seqan global chain of the seeds.
    String<TSeed> seedChain;
    chainSeedsGlobally(seedChain, seedSet, SparseChaining());
    int seedsInChain = length(seedChain);
    if (seedsInChain == 0)
        return cppStringToCString(output + ";Failed: no global seed chain");

    if (verbosity > 3)
    {
        output += "  Globally chained seeds before bridging\n";
        output += getSeedChainTable(seedChain);
    }

    // Now we create a new seed chain will all of the gaps bridged. This will help keep alignment
    // in a narrow band, even when the seeds are spaced apart.
    TSeedSet bridgedSeedSet;

    // Create a seed bridge for the start of the chain by following the slope backwards from the
    // first seed.
    TSeed firstSeed = seedChain[0];
    int firstSeedHStart = beginPositionH(firstSeed);
    int firstSeedVStart = beginPositionV(firstSeed);
    if (firstSeedHStart > 0 && firstSeedVStart > 0)
    {
        double vPosAtStartOfH = firstSeedVStart - (slope * firstSeedHStart);
        if (vPosAtStartOfH >= 0.0)
            addBridgingSeed(bridgedSeedSet, 0, std::round(vPosAtStartOfH), firstSeedHStart, firstSeedVStart);
        else
            addBridgingSeed(bridgedSeedSet, std::round(-vPosAtStartOfH / slope), 0, firstSeedHStart, firstSeedVStart);
    }
    addSeedMerge(bridgedSeedSet, firstSeed);
    
    // Fill in any gaps in the middle of the seed chain.
    for (int i = 1; i < seedsInChain; ++i)
    {
        TSeed seed1 = seedChain[i-1];
        TSeed seed2 = seedChain[i];
        int seed1HEnd = endPositionH(seedChain[i-1]);
        int seed1VEnd = endPositionV(seedChain[i-1]);
        int seed2HStart = beginPositionH(seedChain[i]);
        int seed2VStart = beginPositionV(seedChain[i]);
        addBridgingSeed(bridgedSeedSet, seed1HEnd, seed1VEnd, seed2HStart, seed2VStart);
        addSeedMerge(bridgedSeedSet, seed2);
    }

    // Create a seed bridge for the end of the chain by following the slope forwards from the
    // last seed.
    TSeed lastSeed = seedChain[seedsInChain - 1];
    int lastSeedHEnd = endPositionH(lastSeed);
    int lastSeedVEnd = endPositionV(lastSeed);
    if (lastSeedHEnd < s1Len && lastSeedVEnd < s2Len)
    {
        double vPosAtStartOfH = lastSeedVEnd - (slope * lastSeedHEnd);
        double vPosAtEndOfH = (slope * s1Len) + vPosAtStartOfH;
        if (vPosAtEndOfH <= s2Len)
            addBridgingSeed(bridgedSeedSet, lastSeedHEnd, lastSeedVEnd, s1Len, std::round(vPosAtEndOfH));
        else
            addBridgingSeed(bridgedSeedSet, lastSeedHEnd, lastSeedVEnd, std::round((s2Len - vPosAtStartOfH) / slope), s2Len);
    }

    String<TSeed> bridgedSeedChain;
    chainSeedsGlobally(bridgedSeedChain, bridgedSeedSet, SparseChaining());

    if (verbosity > 3)
    {
        output += "  Seed chain after bridging\n";
        output += getSeedChainTable(bridgedSeedChain);
    }

    Align<Dna5String, ArrayGaps> alignment;
    resize(rows(alignment), 2);
    assignSource(row(alignment, 0), sequenceH);
    assignSource(row(alignment, 1), sequenceV);
    Score<int, Simple> scoringScheme(matchScore, mismatchScore, gapExtensionScore, gapOpenScore);
    AlignConfig<true, true, true, true> alignConfig;
    int score = bandedChainAlignment(alignment, bridgedSeedChain, scoringScheme, alignConfig, bandSize);

    return turnAlignmentIntoDescriptiveString(&alignment, score, matchScore, startTime, output);
}

// This function takes the seed chain which should end at the point hStart, vStart. New seeds will be
// added to the chain to reach the point hEnd, vEnd.
void addBridgingSeed(TSeedSet & seedSet, int hStart, int vStart, int hEnd, int vEnd)
{
    int hDiff = hEnd - hStart;
    int vDiff = vEnd - vStart;

    // If the start and end are against each other, then there's nothing to bridge.
    if (hDiff == 0 || vDiff == 0)
        return;

    TSeed seed(hStart, vStart, hEnd, vEnd);
    addSeedMerge(seedSet, seed);
}

// This function searches for lines in the 2D ref-read space that represent likely semi-global
// alignments.
char * findAlignmentLines(char * s1, char * s2, int s1Len, int s2Len, double expectedSlope,
                          int verbosity)
{
    std::string s1Str(s1);
    std::string s2Str(s2);
    std::string output;

    std::vector<CommonKmer> commonKmers = getCommonKmers(&s1Str, &s2Str, expectedSlope,
                                                         verbosity, output);
    if (commonKmers.size() < 2)
        return cppStringToCString(output + ";Failed: too few common kmers");

    int kSize = commonKmers[0].m_sequence.length();

    double commonKmerDensity = double(commonKmers.size()) / (double(s1Str.length()) * double(s2Str.length()));
    if (verbosity > 3)
        output += "  Common k-mer density: " + std::to_string(commonKmerDensity) + "\n";

    // Sort by rotated vertical position so lines should be roughly horizontal.
    std::sort(commonKmers.begin(), commonKmers.end(), [](const CommonKmer & a, const CommonKmer & b) {
        return a.m_rotatedVPosition < b.m_rotatedVPosition;   
    });

    // Score each point based on the number of other points in its band. The score is scaled by
    // the length of the band so short bands aren't penalised.
    double bandSize = 200.0; // TO DO: MAKE THIS A PARAMETER?
    int windowStart = 0;
    int windowEnd = 0;
    double maxScore = 0.0;
    for (int i = 0; i < commonKmers.size(); ++i)
    {
        double y = commonKmers[i].m_rotatedVPosition;
        while (y - commonKmers[windowStart].m_rotatedVPosition > bandSize)
            ++windowStart;
        while (windowEnd < commonKmers.size() &&
               commonKmers[windowEnd].m_rotatedVPosition - y < bandSize)
            ++windowEnd;
        double bandLength = getLineLength(commonKmers[i].m_hPosition, commonKmers[i].m_vPosition,
                                          expectedSlope, s1Len, s2Len);
        int bandCount = windowEnd - windowStart;
        double bandKmerDensity = bandCount / (bandLength * 2.0 * bandSize);
        double score = bandKmerDensity / commonKmerDensity;
        if (score > maxScore)
            maxScore = score;

        commonKmers[i].m_bandLength = bandLength;
        commonKmers[i].m_bandCount = bandCount;
        commonKmers[i].m_score = score;
    }

    if (verbosity > 3)
    {
        output += "  Common k-mer positions:\n";
        output += getKmerTable(commonKmers);
        output += "  Max score: " + std::to_string(maxScore) + "\n";
    }

    // Now group all of the line points. A line group begins when the score exceeds a threshold and
    // it ends when the score drops below the threshold.
    std::vector<std::vector<CommonKmer> > lineGroups;
    double scoreThreshold = 10.0; // TO DO: MAKE THIS A PARAMETER?
    bool lineInProgress = false;
    for (int i = 0; i < commonKmers.size(); ++i)
    {
        if (commonKmers[i].m_score > scoreThreshold)
        {
            if (!lineInProgress)
            {
                lineGroups.push_back(std::vector<CommonKmer>());
                lineInProgress = true;
            }
            lineGroups.back().push_back(commonKmers[i]);
        }
        else // score is below threshold
            lineInProgress = false;
    }

    // For each line group, throw out any point which is too divergent from its neighbours. To do
    // this, we determine the slope between each point and its nearest neighbours and how divergent
    // it is from the expected slope. We discard the points with the most divergent slopes.

    // Parameters for this filtering step. TO DO: MAKE THESE ADJUSTABLE?
    double fractionToDiscard = 0.05;
    int steps = 3;

    for (int i = 0; i < lineGroups.size(); ++i)
    {
        std::vector<CommonKmer> * lineGroup = &(lineGroups[i]);

        // Sort by rotated horizontal position so we proceed along the line from start to end. 
        std::sort(lineGroup->begin(), lineGroup->end(), [](const CommonKmer & a, const CommonKmer & b) {
            return a.m_rotatedHPosition < b.m_rotatedHPosition;   
        });

        // Store the max slope for each point in the line group.
        std::vector<double> maxSlopes;
        maxSlopes.reserve(lineGroup->size());
        int groupSize = lineGroup->size();
        for (int j = 0; j < groupSize; ++j)
        {
            CommonKmer * thisPoint = &((*lineGroup)[j]);
            double maxSlope = 0.0;
            for (int k = -steps; k <= steps; ++k)
            {
                int neighbourIndex = j + k;
                if (neighbourIndex >= 0 && neighbourIndex < groupSize)
                {
                    CommonKmer * neighbour = &((*lineGroup)[neighbourIndex]);
                    double slope = (thisPoint->m_rotatedVPosition - neighbour->m_rotatedVPosition) /
                                   (thisPoint->m_rotatedHPosition - neighbour->m_rotatedHPosition);
                    maxSlope = std::max(maxSlope, fabs(slope));
                }
            }
            maxSlopes.push_back(maxSlope);
        }

        // Determine a slope cutoff that will exclude the correct fraction of points.
        std::vector<double> sortedMaxSlopes = maxSlopes;
        std::sort(sortedMaxSlopes.begin(), sortedMaxSlopes.end());
        double slopeCutoff = sortedMaxSlopes[int(groupSize * (1.0 - fractionToDiscard))];
  
        // Create a new line group, excluding with excessively high slopes.
        std::vector<CommonKmer> fixedLineGroup;
        for (int j = 0; j < groupSize; ++j)
        {
            double maxSlope = maxSlopes[j];
            if (maxSlope < slopeCutoff)
                fixedLineGroup.push_back((*lineGroup)[j]);
        }
        lineGroups[i] = fixedLineGroup;
    }

    // Remove any line groups with fewer than two points.
    lineGroups.erase(std::remove_if(lineGroups.begin(), lineGroups.end(), 
                                    [](std::vector<CommonKmer> i) {return i.size() < 2;}),
                     lineGroups.end());

    if (lineGroups.size() == 0)
        return cppStringToCString(output + ";Failed: no lines found");
    
    if (verbosity > 2)
        output += "  Lines found:\n";
    
    std::string linesString;
    for (int i = 0; i < lineGroups.size(); ++i)
    {
        if (i > 0)
            linesString += ";";

        // Perform a simple least-squares linear regression on each line group.
        double slope, intercept;
        linearRegression(lineGroups[i], &slope, &intercept);
        linesString += std::to_string(slope) + "," + std::to_string(intercept) + "," + std::to_string(kSize);

        // Add the k-mer locations to the returned string.
        for (int j = 0; j < lineGroups[i].size(); ++j)
            linesString += "," + std::to_string(lineGroups[i][j].m_hPosition) + "," + std::to_string(lineGroups[i][j].m_vPosition);

        if (verbosity > 2)
            output += "    slope = " + std::to_string(slope) + ", intercept = " + std::to_string(intercept) + "\n";

        if (verbosity > 3)
            output += getKmerTable(lineGroups[i]);
    }
    return cppStringToCString(output + ";" + linesString);
}





// This function returns a list of the k-mers common to the two sequences.
std::vector<CommonKmer> getCommonKmers(std::string * s1, std::string * s2, double expectedSlope,
                                       int verbosity, std::string & output)
{
    std::vector<CommonKmer> commonKmers;
    double rotationAngle = CommonKmer::getRotationAngle(expectedSlope);

    // We will dynamically choose a k-mer size that gives a useful density of common locations.
    int targetKCount = double(s1->length()) * double(s2->length()) * 0.00001; // MIGHT NEED TO TUNE THIS CONSTANT
    if (targetKCount < 40)
        targetKCount = 40;
    int minimumKCount = targetKCount / 4;
    if (verbosity > 2)
        output += "  Target k-mer range: " + std::to_string(minimumKCount) + " to " + std::to_string(targetKCount) + "\n";
    int kSize = 10; // Starting k-mer size

    // The order used is based not on s1 or s2 but shorter sequence vs longer sequence.
    std::string * shorter = s1;
    std::string * longer = s2;
    bool seq1IsShorter = true;
    if (s1->length() > s2->length())
    {
        std::swap(shorter, longer);
        seq1IsShorter = false;
    }

    std::map<std::string, std::vector<int> > commonLocations;
    commonLocations = getCommonLocations(shorter, longer, kSize);
    if (verbosity > 2)
        printKmerSize(kSize, commonLocations.size(), output);

    // If the starting k-mer gave too many locations, we increase it until we're in the correct
    // range.
    if (commonLocations.size() > targetKCount)
    {
        while (true)
        {
            ++kSize;
            commonLocations = getCommonLocations(shorter, longer, kSize);
            if (verbosity > 2) printKmerSize(kSize, commonLocations.size(), output);

            // If we're reached the target range, that's good...
            if (commonLocations.size() <= targetKCount)
            {
                // But if we went under the minimum, then we need to back down one k-mer step, even
                // if it takes us over our target.
                if (commonLocations.size() < minimumKCount)
                {
                    --kSize;
                    commonLocations = getCommonLocations(shorter, longer, kSize);
                    if (verbosity > 2) printKmerSize(kSize, commonLocations.size(), output);
                }
                break;
            }
        }
    }

    // If the starting k-mer gave too few locations, we decrease it until we're in the correct
    // range. But if at any step our smaller k-mer size reduces our common k-mer count, then 
    if (commonLocations.size() < minimumKCount)
    {
        int lastCount = commonLocations.size();
        while (true)
        {
            --kSize;
            if (kSize < 2)
                break;
            commonLocations = getCommonLocations(shorter, longer, kSize);
            if (verbosity > 2)
                printKmerSize(kSize, commonLocations.size(), output);
            if (commonLocations.size() >= minimumKCount)
                break;

            // If making the k-mer size smaller reduced the number of common k-mers, then we've
            // gone too far. Back up and break out of the loop.
            if (commonLocations.size() < lastCount)
            {
                ++kSize;
                commonLocations = getCommonLocations(shorter, longer, kSize);
                if (verbosity > 2)
                    printKmerSize(kSize, commonLocations.size(), output);
                break;
            }
            lastCount = commonLocations.size();
        }
    }

    // Build the vector of CommonKmer objects.
    int kCount = shorter->size() - kSize;
    for (int i = 0; i < kCount; ++i)
    {
        std::string kmer = shorter->substr(i, kSize);
        if (commonLocations.find(kmer) != commonLocations.end()) // If kmer is in commonLocations
        {
            int shorterPos = i;
            for (int j = 0; j < commonLocations[kmer].size(); ++j)
            {
                int longerPos = commonLocations[kmer][j];
                if (seq1IsShorter)
                    commonKmers.push_back(CommonKmer(kmer, shorterPos, longerPos, rotationAngle));
                else
                    commonKmers.push_back(CommonKmer(kmer, longerPos, shorterPos, rotationAngle));
            }
        }
    }

    return commonKmers;
}


// This function creates a map where the key is a kmer sequence and the value is a vector of 
// sequence positions for the longer sequence. If the positions are for sequence 1, then it sets
// *seq1Positions to true. If the positions are for sequence 2, then it sets *seq1Positions to
// false.
std::map<std::string, std::vector<int> > getCommonLocations(std::string * shorter, std::string * longer,
                                                            int kSize)
{
    // Build a set of all k-mers in the shorter sequence.
    std::set<std::string> shorterSeqKmers;
    int kCount = shorter->size() - kSize;
    for (int i = 0; i < kCount; ++i)
        shorterSeqKmers.insert(shorter->substr(i, kSize));

    // Loop through the longer sequence, and for all kmers also present in the shorter sequence,
    // add to the map.
    std::map<std::string, std::vector<int> > commonLocations;
    kCount = longer->size() - kSize;
    for (int i = 0; i < kCount; ++i)
    {
        std::string kmer = longer->substr(i, kSize);
        if (shorterSeqKmers.find(kmer) != shorterSeqKmers.end()) // If kmer is in shorterSeqKmers
        {
            if (commonLocations.find(kmer) == commonLocations.end()) // If kmer is not in commonLocations
                commonLocations[kmer] = std::vector<int>();
            commonLocations[kmer].push_back(i);
        }
    }
    return commonLocations;
}





char * turnAlignmentIntoDescriptiveString(Align<Dna5String, ArrayGaps> * alignment, int score,
                                          int matchScore, long long startTime, std::string output)
{
    // Extract the alignment sequences into C++ strings, as the TRow type doesn't seem to have
    // constant time random access.
    std::ostringstream stream1;
    stream1 << row(*alignment, 0);;
    std::string s1Alignment =  stream1.str();
    std::ostringstream stream2;
    stream2 << row(*alignment, 1);;
    std::string s2Alignment =  stream2.str();

    int alignmentLength = std::max(s1Alignment.size(), s2Alignment.size());
    if (alignmentLength == 0)
        return cppStringToCString(output + ";Failed: alignment length zero");

    // Build a CIGAR string of the alignment.
    std::string cigarString;
    CigarType currentCigarType;
    int currentCigarLength = 0;
    int matchCount = 0;
    int mismatchCount = 0;
    int insertionCount = 0;
    int deletionCount = 0;
    std::vector<int> s2MismatchPositions;
    std::vector<int> s2InsertionPositions;
    std::vector<int> s2DeletionPositions;
    int s1Start = -1, s2Start = -1;
    int s1Bases = 0, s2Bases = 0;
    bool alignmentStarted = false;
    bool s1Started = false, s2Started = false;
    for (int i = 0; i < alignmentLength; ++i)
    {
        char base1 = s1Alignment[i];
        char base2 = s2Alignment[i];

        // We consider the alignment to have started when we've encountered a base in both
        // sequences (though not necessarily at the same time).
        if (base1 != '-')
            s1Started = true;
        if (base2 != '-')
            s2Started = true;
        if (s1Started && s2Started && !alignmentStarted)
        {
            s1Start = s1Bases;
            s2Start = s2Bases;
            alignmentStarted = true;
        }

        CigarType cigarType = getCigarType(base1, base2, alignmentStarted);
        if (i == 0)
            currentCigarType = cigarType;

        // Tally up counts and positions for matches, mismatches, insertions and deletions.
        if (cigarType == MATCH)
        {
            if (base1 == base2)
                ++matchCount;
            else
            {
                ++mismatchCount;
                s2MismatchPositions.push_back(s2Bases);
            }
        }
        else if (cigarType == DELETION)
        {
            ++deletionCount;
            s2DeletionPositions.push_back(s2Bases);
        }
        else if (cigarType == INSERTION)
        {
            ++insertionCount;
            s2InsertionPositions.push_back(s2Bases);
        }

        if (cigarType == currentCigarType)
            ++currentCigarLength;
        else
        {
            cigarString.append(getCigarPart(currentCigarType, currentCigarLength));
            currentCigarType = cigarType;
            currentCigarLength = 1;
        }

        if (base1 != '-')
            ++s1Bases;
        if (base2 != '-')
            ++s2Bases;
    }

    int s1End = s1Bases;
    int s2End = s2Bases;
    if (currentCigarType == INSERTION)
    {
        currentCigarType = CLIP;
        insertionCount -= currentCigarLength;
        s2InsertionPositions.resize(insertionCount);
        s1End -= currentCigarLength;
    }
    else if (currentCigarType == DELETION)
    {
        currentCigarType = NOTHING;
        deletionCount -= currentCigarLength;
        s2DeletionPositions.resize(deletionCount);
        s2End -= currentCigarLength;
    }    
    cigarString.append(getCigarPart(currentCigarType, currentCigarLength));

    int editDistance = mismatchCount + insertionCount + deletionCount;
    int alignedLength = matchCount + mismatchCount + insertionCount + deletionCount;
    double percentIdentity = 100.0 * matchCount / alignedLength;

    // Scale the score and adjust to the alignment length. 100 = perfect score.
    double perfectScore = double(matchScore) * double(alignedLength);
    double scaledScore = 100.0 * double(score) / perfectScore;

    long long endTime = duration_cast< milliseconds >(system_clock::now().time_since_epoch()).count(); 
    int milliseconds = endTime - startTime;

    std::string finalString = output + ";" + 
                              cigarString + ";" +
                              std::to_string(s1Start) + ";" + 
                              std::to_string(s1End) + ";" + 
                              std::to_string(s2Start) + ";" + 
                              std::to_string(s2End) + ";" + 
                              std::to_string(alignedLength) + ";" + 
                              std::to_string(matchCount) + ";" + 
                              std::to_string(mismatchCount) + ";" + 
                              vectorToString(&s2MismatchPositions) + ";" + 
                              std::to_string(insertionCount) + ";" + 
                              vectorToString(&s2InsertionPositions) + ";" + 
                              std::to_string(deletionCount) + ";" + 
                              vectorToString(&s2DeletionPositions) + ";" + 
                              std::to_string(editDistance) + ";" + 
                              std::to_string(percentIdentity) + ";" + 
                              std::to_string(score) + ";" + 
                              std::to_string(scaledScore) + ";" + 
                              std::to_string(milliseconds);

    return cppStringToCString(finalString);
}

// Frees dynamically allocated memory for a c string. Called by Python after the string has been
// received.
void free_c_string(char * p)
{
    free(p);
}

char * cppStringToCString(std::string cpp_string)
{
    char * c_string = (char*)malloc(sizeof(char) * (cpp_string.size() + 1));
    std::copy(cpp_string.begin(), cpp_string.end(), c_string);
    c_string[cpp_string.size()] = '\0';
    return c_string;
}

std::string vectorToString(std::vector<int> * v)
{
    std::stringstream ss;
    for(size_t i = 0; i < v->size(); ++i)
    {
        if (i != 0)
            ss << ",";
        ss << (*v)[i];
    }
    return ss.str();
}

CigarType getCigarType(char b1, char b2, bool alignmentStarted)
{
    if (b1 == '-')
    {
        if (alignmentStarted)
            return DELETION;
        else
            return NOTHING;
    }
    else if (b2 == '-')
    {
        if (alignmentStarted)
            return INSERTION;
        else
            return CLIP;
    }
    else
        return MATCH;
}

std::string getCigarPart(CigarType type, int length)
{
    std::string cigarPart = std::to_string(length);
    if (type == DELETION)
        cigarPart.append("D");
    else if (type == INSERTION)
        cigarPart.append("I");
    else if (type == CLIP)
        cigarPart.append("S");
    else if (type == MATCH)
        cigarPart.append("M");
    else //type == NOTHING
        return "";
    return cigarPart;
}


// This function gets the slope and intercept of the seed chain line as defined by the seeds at the
// given indices. Slope is defined as the change in V position over the change in H position.
// Intercept is defined as the V position at the H position of 0.
void getSeedChainSlopeAndIntercept(String<TSeed> * seedChain, int i1, int i2,
                                   double * slope, double * intercept)
{
    int hStart = beginPositionH((*seedChain)[i1]);
    int hEnd = endPositionH((*seedChain)[i2]);
    int vStart = beginPositionV((*seedChain)[i1]);
    int vEnd = endPositionV((*seedChain)[i2]);

    getSlopeAndIntercept(hStart, hEnd, vStart, vEnd, slope, intercept);
}


void getSlopeAndIntercept(int hStart, int hEnd, int vStart, int vEnd,
                          double * slope, double * intercept)
{
    *slope = -1.0;
    *intercept = -1.0;

    int hDiff = hEnd - hStart;
    int vDiff = vEnd - vStart;

    if (hDiff > 0)
    {
        *slope = double(vDiff) / hDiff;
        *intercept = vStart - (*slope * hStart);
    }
}

double getMedian(std::vector<double> * v)
{
    size_t size = v->size();
    std::sort(v->begin(), v->end());
    if (size % 2 == 0)
        return ((*v)[size / 2 - 1] + (*v)[size / 2]) / 2;
    else 
        return (*v)[size / 2];
}


// Creates the CommonKmer object using the position in the two sequences.
// It then rotates the point relative to the origin by the given angle (in radians).
CommonKmer::CommonKmer(std::string sequence, int hPosition, int vPosition, double angle) :
    m_sequence(sequence),
    m_hPosition(hPosition),
    m_vPosition(vPosition),
    m_score(0.0)
{
    double s = sin(angle);
    double c = cos(angle);
    m_rotatedHPosition = (m_hPosition * c) - (m_vPosition * s);
    m_rotatedVPosition = (m_hPosition * s) + (m_vPosition * c);
}

void printKmerSize(int kmerSize, int locationCount, std::string & output)
{
    output += "  " + std::to_string(locationCount) + " " + std::to_string(kmerSize) + "-mers in common\n";
}

// Given a point, a slope and rectangle bounds, this function returns the length of the line
// segment which goes through that point with that slope and is in those bounds.
double getLineLength(double x, double y, double slope, double xSize, double ySize)
{
    double xStart, yStart, xEnd, yEnd;

    double yIntercept = y - (slope * x);
    if (yIntercept >= 0.0)
    {
        xStart = 0.0;
        yStart = yIntercept;
    }
    else
    {
        xStart = -yIntercept / slope;
        yStart = 0.0;
    }

    double yAtXEnd = (slope * xSize) + yIntercept;
    if (yAtXEnd <= ySize)
    {
        xEnd = xSize;
        yEnd = yAtXEnd;
    }
    else
    {
        xEnd = (ySize - yIntercept) / slope;
        yEnd = ySize;
    }

    double xLength = xEnd - xStart;
    double yLength = yEnd - yStart;
    return sqrt((xLength * xLength) + (yLength * yLength));
}


// Adapted from:
// http://stackoverflow.com/questions/11449617/how-to-fit-the-2d-scatter-data-with-a-line-with-c
void linearRegression(std::vector<CommonKmer> & pts, double * slope, double * intercept)
{
    int n = pts.size();
    double sumH = 0.0, sumV = 0.0, sumHV = 0.0, sumHH = 0.0;
    for (int i = 0; i < n; ++i)
    {
        double hPos = pts[i].m_hPosition;
        double vPos = pts[i].m_vPosition;
        sumH += hPos;
        sumV += vPos;
        sumHV += hPos * vPos;
        sumHH += hPos * hPos;
    }
    double hMean = sumH / n;
    double vMean = sumV / n;
    double denominator = sumHH - (sumH * hMean);
    *slope = (sumHV - sumH * vMean) / denominator;
    *intercept = vMean - (*slope * hMean);
}


void parseKmerLocationsFromString(std::string & str, std::vector<int> & v1, std::vector<int> & v2)
{
    std::stringstream ss(str);
    while(ss.good())
    {
        std::string s1Location;
        std::string s2Location;
        std::getline(ss, s1Location, ',');
        std::getline(ss, s2Location, ',');
        v1.push_back(std::stoi(s1Location));
        v2.push_back(std::stoi(s2Location));
    }
}

std::string getKmerTable(std::vector<CommonKmer> & commonKmers)
{
    std::string table = "\tSeq 1 pos\tSeq 2 pos\tRotated seq 1 pos\tRotated seq 2 pos\tBand length\tBand count\tScore\n";
    for (int i = 0; i < commonKmers.size(); ++i)
    {
        table += "\t" + std::to_string(commonKmers[i].m_hPosition) +
                 "\t" + std::to_string(commonKmers[i].m_vPosition) + 
                 "\t" + std::to_string(commonKmers[i].m_rotatedHPosition) +
                 "\t" + std::to_string(commonKmers[i].m_rotatedVPosition) +
                 "\t" + std::to_string(commonKmers[i].m_bandLength) +
                 "\t" + std::to_string(commonKmers[i].m_bandCount) +
                 "\t" + std::to_string(commonKmers[i].m_score) + "\n";
    }
    return table;
}

std::string getSeedChainTable(String<TSeed> & seedChain)
{
    std::string table = "\tH start\tH end\tV start\tV end\n";
    int seedsInChain = length(seedChain);
    for (int i = 0; i < seedsInChain; ++i)
    {
        table += "\t" + std::to_string(beginPositionH(seedChain[i])) +
                 "\t" + std::to_string(endPositionH(seedChain[i])) +
                 "\t" + std::to_string(beginPositionV(seedChain[i])) +
                 "\t" + std::to_string(endPositionV(seedChain[i])) + "\n";
    }
    return table;
}

// This function adds a seed to a seed set. First it tries to merge the seed, and if that doesn't
// work it adds it as a separate seed.
void addSeedMerge(TSeedSet & seedSet, TSeed & seed)
{
    if (!addSeed(seedSet, seed, 1, Merge()))
        addSeed(seedSet, seed, Single());
}

}
