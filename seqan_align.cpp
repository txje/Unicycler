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


std::vector<Kmer> getSeqKmers(std::string seq, int strLen, int kSize);
std::vector<CommonLocation> getCommonLocations(std::vector<Kmer> s1Kmers, std::vector<Kmer> s2Kmers);
CigarType getCigarType(char b1, char b2, bool alignmentStarted);
std::string getCigarPart(CigarType type, int length);
char * cpp_string_to_c_string(std::string cpp_string);


char * semiGlobalAlign(char * s1, char * s2, int s1Len, int s2Len, int kSize, int bandSize)
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

	String<TSeed> seedChain;
	chainSeedsGlobally(seedChain, seedSet, SparseChaining());
	if (length(seedChain) == 0)
		return strdup("");

	Align<Dna5String, ArrayGaps> alignment;
    resize(rows(alignment), 2);
    assignSource(row(alignment, 0), sequenceH);
    assignSource(row(alignment, 1), sequenceV);
    Score<int, Simple> scoringScheme(1, -1, -1);
    AlignConfig<true, true, true, true> alignConfig;
    bandedChainAlignment(alignment, seedChain, scoringScheme, alignConfig, bandSize);

    // Extract the alignment sequences into C++ strings, as the TRow type doesn't seem to have
    // constant time random access.
	std::ostringstream stream1;
    stream1 << row(alignment, 0);;
    std::string s1Alignment =  stream1.str();
	std::ostringstream stream2;
    stream2 << row(alignment, 1);;
    std::string s2Alignment =  stream2.str();

	int alignmentLength = std::max(s1Alignment.size(), s2Alignment.size());
	if (alignmentLength == 0)
		return strdup("");

	// Build a CIGAR string of the alignment.
	std::string cigarString;
	CigarType currentCigarType;
	int currentCigarLength = 0;
	int matchCount = 0;
	int mismatchCount = 0;
	int insertionCount = 0;
	int deletionCount = 0;
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
		if (base1 != '-')
			++s1Bases;
		if (base2 != '-')
			++s2Bases;

		CigarType cigarType = getCigarType(base1, base2, alignmentStarted);
		if (cigarType == MATCH)
		{
			if (base1 == base2)
				++matchCount;
			else
				++mismatchCount;
		}
		if (cigarType == DELETION)
			++deletionCount;
		if (cigarType == INSERTION)
			++insertionCount;

		if (i == 0)
			currentCigarType = cigarType;
		if (cigarType == currentCigarType)
			++currentCigarLength;
		else
		{
			cigarString.append(getCigarPart(currentCigarType, currentCigarLength));
			currentCigarType = cigarType;
			currentCigarLength = 1;
		}
	}

	int s1End = s1Bases;
	int s2End = s2Bases;
	if (currentCigarType == INSERTION)
	{
		currentCigarType = CLIP;
		insertionCount -= currentCigarLength;
		s1End -= currentCigarLength;
	}
	else if (currentCigarType == DELETION)
	{
		currentCigarType = NOTHING;
		deletionCount -= currentCigarLength;
		s2End -= currentCigarLength;
	}	
	cigarString.append(getCigarPart(currentCigarType, currentCigarLength));

    long long endTime = duration_cast< milliseconds >(system_clock::now().time_since_epoch()).count();
	int editDistance = mismatchCount + insertionCount + deletionCount;
	int alignedLength = matchCount + mismatchCount + insertionCount + deletionCount;
	double percentIdentity = 100.0 * matchCount / alignedLength;
	int milliseconds = endTime - startTime;

	std::string finalString = cigarString + "," +
	                          std::to_string(s1Start) + "," + std::to_string(s1End) + "," + 
	                          std::to_string(s2Start) + "," + std::to_string(s2End) + "," + 
	                          std::to_string(alignedLength) + "," + 
	                          std::to_string(matchCount) + "," + std::to_string(mismatchCount) + "," + 
	                          std::to_string(insertionCount) + "," + std::to_string(deletionCount) + "," + 
	                          std::to_string(editDistance) + "," + std::to_string(percentIdentity) + "," + 
	                          std::to_string(milliseconds);



	return cpp_string_to_c_string(finalString);
}

// Frees dynamically allocated memory for a c string. Called by Python after the string has been
// received.
void free_c_string(char * p)
{
    free(p);
}

char * cpp_string_to_c_string(std::string cpp_string)
{
	char * c_string = (char*)malloc(sizeof(char) * (cpp_string.size() + 1));
	std::copy(cpp_string.begin(), cpp_string.end(), c_string);
	c_string[cpp_string.size()] = '\0';
	return c_string;
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
	std::string cigarPart;
	if (type == DELETION)
		cigarPart = "D";
	else if (type == INSERTION)
		cigarPart = "I";
	else if (type == CLIP)
		cigarPart = "S";
	else if (type == MATCH)
		cigarPart = "M";
	else //type == NOTHING
		return "";
	cigarPart.append(std::to_string(length));
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

}
