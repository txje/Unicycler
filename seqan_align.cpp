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
typedef std::tuple<std::string, int, int> Kmer;
typedef std::map<std::string, std::tuple<int, int>> KmerDict;
typedef std::tuple<int, int, int, int> CommonLocation;

enum CigarType {MATCH, INSERTION, DELETION, CLIP, NOTHING};


char * turnAlignmentIntoDescriptiveString(Align<Dna5String, ArrayGaps> * alignment,
                                          long long startTime);
std::vector<Kmer> getSeqKmers(std::string seq, int strLen, int kSize);
std::vector<CommonLocation> getCommonLocations(std::vector<Kmer> s1Kmers,
                                               std::vector<Kmer> s2Kmers);
CigarType getCigarType(char b1, char b2, bool alignmentStarted);
std::string getCigarPart(CigarType type, int length);
char * cppStringToCString(std::string cpp_string);
std::string vectorToString(std::vector<int> * v);
double getSeedChainMedianSlope(String<TSeed> * seedChain);
double getMedian(std::vector<double> * v);



//
// kSize = the kmer size used to find alignment seeds.
// bandSize = the margin around seeds used for alignment. Larger values are more likely to find the
//         best alignment, at a performance cost.
// allowedLengthDiscrepancy = how much the sequences are allowed to vary in length as judged by the
//         seed chain. E.g. 0.1 means that a ratio between 0.9 and 1.11111 is acceptable. Anything
//         outside that ratio will not be aligned. This value must be in the range [0, 1)
char * semiGlobalAlign(char * s1, char * s2, int s1Len, int s2Len, int kSize, int bandSize,
                       double allowedLengthDiscrepancy)
{
    long long startTime = duration_cast< milliseconds >(system_clock::now().time_since_epoch()).count();

    std::string s1Str(s1);
    std::string s2Str(s2);
    std::vector<Kmer> s1Kmers = getSeqKmers(s1Str, s1Len, kSize);
    std::vector<Kmer> s2Kmers = getSeqKmers(s2Str, s2Len, kSize);
    std::vector<CommonLocation> commonLocations = getCommonLocations(s1Kmers, s2Kmers);

    TSeedSet seedSet;
    for (int i = 0; i < commonLocations.size(); ++i)
    {
        CommonLocation l = commonLocations[i];
        TSeed seed(std::get<0>(l), std::get<2>(l), std::get<1>(l), std::get<3>(l));
        if (!addSeed(seedSet, seed, 1, Merge()))
            addSeed(seedSet, seed, Single());
    }

    TSequence sequenceH = s1;
    TSequence sequenceV = s2;

    // // Debugging code: output seed positions before chaining
    // std::cout << std::endl << "Seed positions before global chaining" << std::endl;
    // for (Iterator<TSeedSet>::Type it = begin(seedSet, Standard()); it != end(seedSet, Standard()); ++it)
    // {
    //     std::cout << beginPositionH(*it) << "\t" << beginPositionV(*it) << std::endl;
    //     std::cout << endPositionH(*it) << "\t" << endPositionV(*it) << std::endl;
    // }

    String<TSeed> seedChain;
    chainSeedsGlobally(seedChain, seedSet, SparseChaining());

    // Quit before doing the alignment if the seed chain doesn't look good.
    if (length(seedChain) == 0)
        return strdup("Failed: no global seed chain");

    int seedsInChain = length(seedChain);
    int seq1Span = endPositionH(seedChain[seedsInChain-1]) - beginPositionH(seedChain[0]);
    int seq2Span = endPositionV(seedChain[seedsInChain-1]) - beginPositionV(seedChain[0]);
    if (seq2Span == 0)
        return strdup("Failed: sequence 2 span zero");

    // // Debugging code: output seed positions after chaining
    // std::cout << std::endl << "Seed positions after global chaining" << std::endl;
    // for (int i = 0; i < seedsInChain; ++i)
    // {
    //     std::cout << beginPositionH(seedChain[i]) << "\t" << beginPositionV(seedChain[i]) << std::endl; //TEMP
    //     std::cout << endPositionH(seedChain[i]) << "\t" << endPositionV(seedChain[i]) << std::endl; //TEMP
    // }

    double ratio = getSeedChainMedianSlope(&seedChain);
    double minRatio = 1.0 - allowedLengthDiscrepancy;
    double maxRatio = 1.0 / minRatio;
    if (ratio < minRatio || ratio > maxRatio)
        return strdup("Failed: bad sequence length ratio");

    Align<Dna5String, ArrayGaps> alignment;
    resize(rows(alignment), 2);
    assignSource(row(alignment, 0), sequenceH);
    assignSource(row(alignment, 1), sequenceV);
    Score<int, Simple> scoringScheme(1, -1, -1);
    AlignConfig<true, true, true, true> alignConfig;
    bandedChainAlignment(alignment, seedChain, scoringScheme, alignConfig, bandSize);

    return turnAlignmentIntoDescriptiveString(&alignment, startTime);
}


char * semiGlobalAlignExact(char * s1, char * s2, int s1Len, int s2Len)
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

// Returns a list of all Kmers in a sequence.
std::vector<Kmer> getSeqKmers(std::string seq, int strLen, int kSize)
{
    std::vector<Kmer> kmers;
    int kCount = strLen - kSize;
    kmers.reserve(kCount);
    for (int i = 0; i < kCount; ++i)
    {
        Kmer kmer(seq.substr(i, kSize), i, i + kSize);
        kmers.push_back(kmer);
    }
    return kmers;
}

// Returns a list of all Kmers common to both lists.
std::vector<CommonLocation> getCommonLocations(std::vector<Kmer> s1Kmers, std::vector<Kmer> s2Kmers)
{
    std::vector<CommonLocation> commonLocations;

    // Store all s1 kmers in a map of seq -> positions.
    KmerDict s1KmerPositions;
    for (int i = 0; i < s1Kmers.size(); ++i)
    {
        std::string s1KmerSeq = std::get<0>(s1Kmers[i]);
        std::tuple<int, int> s1Positions(std::get<1>(s1Kmers[i]), std::get<2>(s1Kmers[i]));
        s1KmerPositions[s1KmerSeq] = s1Positions;
    }

    // For all s2 kmers, see if they are in the s1 map. If so, they are common.
    for (int i = 0; i < s2Kmers.size(); ++i)
    {
        std::string s2KmerSeq = std::get<0>(s2Kmers[i]);
        if (s1KmerPositions.count(s2KmerSeq))
        {
            std::tuple<int, int> s1Position = s1KmerPositions[s2KmerSeq];
            int s1Start = std::get<0>(s1Position);
            int s1End = std::get<1>(s1Position);
            int s2Start = std::get<1>(s2Kmers[i]);
            int s2End = std::get<2>(s2Kmers[i]);
            CommonLocation commonLocation(s1Start, s1End, s2Start, s2End);
            commonLocations.push_back(commonLocation);
        }
    }

    return commonLocations;
}

// Looks at all of the seeds in the chain and gets the slopes (defined as the slope from the start
// of one seed to the start of the next). Returns the median.
double getSeedChainMedianSlope(String<TSeed> * seedChain)
{
    int seedsInChain = length(*seedChain);
    if (seedsInChain == 0)
        return 0.0;

    int lastHStart = beginPositionH((*seedChain)[0]);
    int lastVStart = beginPositionV((*seedChain)[0]);

    std::vector<double> slopes;
    for (int i = 1; i < seedsInChain; ++i)
    {
        int hStart = beginPositionH((*seedChain)[i]);
        int vStart = beginPositionV((*seedChain)[i]);
        int hDiff = hStart - lastHStart;
        int vDiff = vStart - lastVStart;
        if (vDiff > 0)
        	slopes.push_back(double(hDiff) / double(vDiff));
    }

    if (slopes.size() == 0)
        return 0.0;
    return getMedian(&slopes);
}

double getMedian(std::vector<double> * v)
{
    size_t size = v->size();
    sort(v->begin(), v->end());
    if (size % 2 == 0)
        return ((*v)[size / 2 - 1] + (*v)[size / 2]) / 2.0;
    else
        return (*v)[size / 2];
}

}
