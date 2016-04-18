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
    bool operator<(const CommonKmer& other);

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
char * semiGlobalAlignmentAroundLine(char * s1, char * s2, int s1Len, int s2Len, double slope,
                                     double intercept, int bandSize, int debugOutput);
char * findAlignmentLines(char * s1, char * s2, int s1Len, int s2Len, double expectedSlope,
                         int debugOutput);
char * exhaustiveSemiGlobalAlignment(char * s1, char * s2, int s1Len, int s2Len);
void free_c_string(char * p);

// These functions are internal to this C++ code.
std::vector<CommonKmer> getCommonKmers(std::string * s1, std::string * s2, double expectedSlope,
                                       int debugOutput);
std::map<std::string, std::vector<int> > getCommonLocations(std::string * s1, std::string * s2, int kSize);
char * turnAlignmentIntoDescriptiveString(Align<Dna5String, ArrayGaps> * alignment,
                                          long long startTime);
CigarType getCigarType(char b1, char b2, bool alignmentStarted);
std::string getCigarPart(CigarType type, int length);
char * cppStringToCString(std::string cpp_string);
std::string vectorToString(std::vector<int> * v);
void fixSeedChainToLine(String<TSeed> * seedChain, double bandSize);
void getSeedChainSlopeAndIntercept(String<TSeed> * seedChain, int firstI, int lastI,
                                   double * slope, double * intercept);
void getSlopeAndIntercept(int hStart, int hEnd, int vStart, int vEnd,
                                   double * slope, double * intercept);
char * getSlopeAndInterceptString(double slope, double intercept);
double getMedian(std::vector<double> * v);
void printKmerSize(int kmerSize, int locationCount);
double getLineLength(double x, double y, double slope, double xSize, double ySize);




// This function does the full semi-global alignment using the entirety of both sequences. It will
// be slow but will always find the ideal alignment.
char * exhaustiveSemiGlobalAlignment(char * s1, char * s2, int s1Len, int s2Len)
{
    long long startTime = duration_cast< milliseconds >(system_clock::now().time_since_epoch()).count();

    TSequence sequenceH = s1;
    TSequence sequenceV = s2;

    Align<Dna5String, ArrayGaps> alignment;
    resize(rows(alignment), 2);
    assignSource(row(alignment, 0), sequenceH);
    assignSource(row(alignment, 1), sequenceV);
    Score<int, Simple> scoringScheme(1, -1, -1);
    AlignConfig<true, true, true, true> alignConfig;
    globalAlignment(alignment, scoringScheme, alignConfig);

    return turnAlignmentIntoDescriptiveString(&alignment, startTime);
}




// This function, given a line, will search for semi-global alignments around that line. The
// bandSize parameter specifies how far of an area around the line is searched. It is generally
// much faster than exhaustiveSemiGlobalAlignment, though it may not find the optimal alignment.
// A lower bandSize is faster with a larger chance of missing the optimal alignment.
char * semiGlobalAlignmentAroundLine(char * s1, char * s2, int s1Len, int s2Len, double slope,
                                     double intercept, int bandSize, int debugOutput)
{
    long long startTime = duration_cast< milliseconds >(system_clock::now().time_since_epoch()).count();

    // Extreme slope values will not work, so check for them now.
    if (slope > 1.5)
        return strdup("Failed: slope too large");
    if (slope < 0.5)
        return strdup("Failed: slope too small");

    TSequence sequenceH = s1;
    TSequence sequenceV = s2;

    // Convert the slope and intercept into start/end coordinates for the two sequences.
    int hStart, hEnd, vStart, vEnd;
    if (intercept >= 0.0)
    {
        hStart = 0;
        vStart = std::round(intercept);
    }
    else
    {
        hStart = std::round(-intercept / slope);
        vStart = 0;
    }
    double vPosAtEndOfH = (slope * s1Len) + intercept;
    if (vPosAtEndOfH <= s2Len)
    {
        hEnd = s1Len;
        vEnd = std::round(vPosAtEndOfH);
    }
    else
    {
        hEnd = std::round((s2Len - intercept) / slope);
        vEnd = s2Len;
    }
    int hDiff = hEnd - hStart;
    int vDiff = vEnd - vStart;

    if (debugOutput > 0)
    {
        std::cout << std::endl;
        std::cout << "ALIGNMENT RANGE FROM SLOPE" << std::endl;
        std::cout << "--------------------------" << std::endl;
        std::cout << "hStart: " << hStart << ", " << "hEnd: " << hEnd << ", " << "vStart: " << vStart << ", " << "vEnd: " << vEnd << std::endl << std::endl;;
    }

    // Now that we have an estimate for the start/end of both sequences, we build a single seed
    // chain around that line.
    String<TSeed> seedChain;

    // If the two sequences are the same length, then a single seed will do.
    if (hDiff == vDiff)
        appendValue(seedChain, TSeed(hStart, vStart, hEnd, vEnd));

    // If one of the distances is longer, then we need multiple seeds to step across at a slope
    // not equal to 1.
    else
    {
        int additionalSteps = std::abs(vDiff - hDiff);
        bool additionalHSteps = hDiff > vDiff;

        int seedCount = additionalSteps + 1;
        double pieceSize = double(std::min(hDiff, vDiff)) / seedCount;        
        double targetPosition = 0;
        int distanceSoFar = 0;
        int hPos = hStart;
        int vPos = vStart;
        for (int i = 0; i < seedCount; ++i)
        {
            targetPosition += pieceSize;
            int thisPieceSize = std::round(targetPosition - distanceSoFar);
            distanceSoFar += thisPieceSize;
            appendValue(seedChain, TSeed(hPos, vPos, hPos + thisPieceSize, vPos + thisPieceSize));
            hPos += thisPieceSize;
            vPos += thisPieceSize;
            if (additionalHSteps)
                hPos += 1;
            else
                vPos += 1;
        }
    }

    if (debugOutput > 1)
    {
        int seedsInChain = length(seedChain);
        std::cout << std::endl;
        std::cout << "SEED POSITIONS" << std::endl;
        std::cout << "--------------" << std::endl;
        std::cout << "H start\tH end\tV start\tV end" << std::endl;
        for (int i = 0; i < seedsInChain; ++i)
        {
            std::cout << beginPositionH(seedChain[i]) << "\t" << endPositionH(seedChain[i]) << "\t";
            std::cout << beginPositionV(seedChain[i]) << "\t" << endPositionV(seedChain[i]) << std::endl;
        }
        std::cout << std::endl;
    }

    Align<Dna5String, ArrayGaps> alignment;
    resize(rows(alignment), 2);
    assignSource(row(alignment, 0), sequenceH);
    assignSource(row(alignment, 1), sequenceV);
    Score<int, Simple> scoringScheme(1, -1, -1);
    AlignConfig<true, true, true, true> alignConfig;
    bandedChainAlignment(alignment, seedChain, scoringScheme, alignConfig, bandSize);

    return turnAlignmentIntoDescriptiveString(&alignment, startTime);
}




// This function searches for lines in the 2D ref-read space that represent likely semi-global
// alignments.
char * findAlignmentLines(char * s1, char * s2, int s1Len, int s2Len, double expectedSlope,
                          int debugOutput)
{
    std::string s1Str(s1);
    std::string s2Str(s2);

    std::vector<CommonKmer> commonKmers = getCommonKmers(&s1Str, &s2Str, expectedSlope,
                                                         debugOutput);
    if (commonKmers.size() < 2)
        return strdup("Failed: too few common kmers");

    double commonKmerDensity = double(commonKmers.size()) / (double(s1Str.length()) * double(s2Str.length()));
    if (debugOutput > 1)
        std::cout << std::endl << "COMMON K-MER DENSITY: " << commonKmerDensity << std::endl << std::endl;

    std::sort(commonKmers.begin(), commonKmers.end());

    // Score each point based on the number of other points in its band. The score is scaled by
    // the length of the band so short bands aren't penalised.
    double bandSize = 20.0; // TO DO: MAKE THIS A PARAMETER?
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

    if (debugOutput > 1)
    {
        std::cout << std::endl;
        std::cout << "COMMON K-MER POSITIONS" << std::endl;
        std::cout << "----------------------" << std::endl;
        std::cout << "Seq 1 pos\tSeq 2 pos\tRotated seq 1 pos\tRotated seq 2 pos\tBand length\tBand count\tScore" << std::endl;
        for (int i = 0; i < commonKmers.size(); ++i)
            std::cout << commonKmers[i].m_hPosition << "\t" << commonKmers[i].m_vPosition << "\t"
                      << commonKmers[i].m_rotatedHPosition << "\t" << commonKmers[i].m_rotatedVPosition
                      << "\t" << commonKmers[i].m_bandLength << "\t" << commonKmers[i].m_bandCount
                      << "\t" << commonKmers[i].m_score << std::endl;
        std::cout << std::endl;
        std::cout << "MAX SCORE: " << maxScore << std::endl << std::endl;
    }

    // Now group all of the line points. A line group begins when the score exceeds a threshold and
    // it ends when the score drops below the threshold.
    std::vector<std::vector<CommonKmer> > lineGroups;
    double scoreThreshold = 100.0; // TO DO: MAKE THIS A PARAMETER?
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
        {
            //If the previous group had a single point, remove it as lines need at least two
            // points.
            if (lineInProgress && lineGroups.back().size() == 1)
                lineGroups.pop_back();
            lineInProgress = false;
        }
    }

    if (debugOutput > 0)
        std::cout << std::endl << "LINES FOUND: " << lineGroups.size() << std::endl;

    

    // TO DO
    // TO DO
    // TO DO
    // TO DO
    // TO DO
    // TO DO
    // TO DO
    // TO DO
    // TO DO
    // TO DO
    // TO DO
    // TO DO
    // TO DO
    // TO DO
    // TO DO




    return getSlopeAndInterceptString(1.23, 4.56); //TEMP
}





// This function returns a list of the k-mers common to the two sequences.
std::vector<CommonKmer> getCommonKmers(std::string * s1, std::string * s2, double expectedSlope,
                                       int debugOutput)
{
    std::vector<CommonKmer> commonKmers;
    double rotationAngle = CommonKmer::getRotationAngle(expectedSlope);

    // We will dynamically choose a k-mer size that gives a useful density of common locations.
    int targetKCount = double(s1->length()) * double(s2->length()) * 0.000005; // MIGHT NEED TO TUNE THIS CONSTANT
    if (targetKCount < 10)
    {
        if (debugOutput > 0)
            std::cout << std::endl << "COMMON KMER FAILURE: SEQUENCES TOO SHORT " << std::endl;
        return commonKmers;
    }
    int minimumKCount = targetKCount / 4;
    if (debugOutput > 0)
        std::cout << std::endl << "TARGET K-MER RANGE: " << minimumKCount << " to " <<
                     targetKCount << std::endl;
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
    if (debugOutput > 0) printKmerSize(kSize, commonLocations.size());

    // If the starting k-mer gave too many locations, we increase it until we're in the correct
    // range.
    if (commonLocations.size() > targetKCount)
    {
        while (true)
        {
            ++kSize;
            commonLocations = getCommonLocations(shorter, longer, kSize);
            if (debugOutput > 0) printKmerSize(kSize, commonLocations.size());

            // If we're reached the target range, that's good...
            if (commonLocations.size() <= targetKCount)
            {
                // But if we went under the minimum, then we need to back down one k-mer step, even
                // if it takes us over our target.
                if (commonLocations.size() < minimumKCount)
                {
                    --kSize;
                    commonLocations = getCommonLocations(shorter, longer, kSize);
                    if (debugOutput > 0) printKmerSize(kSize, commonLocations.size());
                }
                break;
            }
        }
    }

    // If the starting k-mer gave too few locations, we decrease it until we're in the correct
    // range.
    if (commonLocations.size() < minimumKCount)
    {
        while (true)
        {
            --kSize;
            if (kSize < 2)
                break;
            commonLocations = getCommonLocations(shorter, longer, kSize);
            if (debugOutput > 0)
                printKmerSize(kSize, commonLocations.size());
            if (commonLocations.size() >= minimumKCount)
                break;
        }
    }
    if (debugOutput > 0) std::cout << std::endl;

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





char * turnAlignmentIntoDescriptiveString(Align<Dna5String, ArrayGaps> * alignment,
                                          long long startTime)
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
        return strdup("Failed: alignment length zero");

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
    for (int i = 0; i < alignmentLength; ++i)
    {
        char base1 = s1Alignment[i];
        char base2 = s2Alignment[i];

        if (base1 != '-' && base2 != '-' && !alignmentStarted)
        {
            s1Start = s1Bases;
            s2Start = s2Bases;
            alignmentStarted = true;
        }

        CigarType cigarType = getCigarType(base1, base2, alignmentStarted);
        if (i == 0)
            currentCigarType = cigarType;

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

    long long endTime = duration_cast< milliseconds >(system_clock::now().time_since_epoch()).count(); 
    int milliseconds = endTime - startTime;

    std::string finalString = cigarString + "," +
                              std::to_string(s1Start) + "," + 
                              std::to_string(s1End) + "," + 
                              std::to_string(s2Start) + "," + 
                              std::to_string(s2End) + "," + 
                              std::to_string(alignedLength) + "," + 
                              std::to_string(matchCount) + "," + 
                              std::to_string(mismatchCount) + "," + 
                              vectorToString(&s2MismatchPositions) + "," + 
                              std::to_string(insertionCount) + "," + 
                              vectorToString(&s2InsertionPositions) + "," + 
                              std::to_string(deletionCount) + "," + 
                              vectorToString(&s2DeletionPositions) + "," + 
                              std::to_string(editDistance) + "," + 
                              std::to_string(percentIdentity) + "," + 
                              std::to_string(milliseconds);

    return cppStringToCString(finalString);
}

char * getSlopeAndInterceptString(double slope, double intercept)
{
    std::string slopeAndInterceptString = std::to_string(slope) + "," + std::to_string(intercept);
    return cppStringToCString(slopeAndInterceptString);

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
            ss << ";";
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

bool CommonKmer::operator<(const CommonKmer& other)
{
    return this->m_rotatedVPosition < other.m_rotatedVPosition;
}


void printKmerSize(int kmerSize, int locationCount)
{
    std::cout << "K-MER: " << kmerSize << ", " << locationCount << " SITES IN COMMON" << std::endl;
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

}
