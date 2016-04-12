//#include <chrono>
#include <seqan/basic.h>
#include <seqan/align.h>
#include <seqan/sequence.h>
#include <seqan/stream.h>  // for printint strings

using namespace seqan;
using namespace std::chrono;


//WILL RETURN A PYTHON STRING LATER
std::string semiGlobalAlign(char * s1, char * s2)
{
    typedef Dna5String TSequence;
    typedef Align<Dna5String, ArrayGaps> TAlign;
    typedef Row<TAlign>::Type TRow;
//    typedef Iterator<TRow>::Type TRowIterator;

    TSequence seq1 = s1;
    TSequence seq2 = s2;

    TAlign align;
    resize(rows(align), 2);
    assignSource(row(align, 0), seq1);
    assignSource(row(align, 1), seq2);

    Score<int, Simple> scoringScheme(1, -1, -1);
    AlignConfig<true, true, true, true> alignConfig;

//    long long ms_before = duration_cast< milliseconds >(system_clock::now().time_since_epoch()).count();

    int score = globalAlignment(align, scoringScheme, alignConfig);

//    long long ms_after = duration_cast< milliseconds >(system_clock::now().time_since_epoch()).count();
//    long long elapsedTime = ms_after - ms_before;

    int alignmentLength = _max(length(row(align, 0)), length(row(align, 1)));

//    std::cout << "seq1 size: " << length(seq1) << std::endl;
//    std::cout << "seq2 size: " << length(seq2) << std::endl;
//    std::cout << "Score: " << score << std::endl;
//    std::cout << "Milliseconds: " << elapsedTime << std::endl;
//    std::cout << align << std::endl;

    //Get the Seqan results in a normal string.
    std::ostringstream stream1;
    std::ostringstream stream2;
    std::cout << align;
    stream1 << row(align, 0);
    stream2 << row(align, 1);
    std::string alignedSeq1 =  stream1.str();
    std::string alignedSeq2 =  stream2.str();

    int seq1FirstBase = -1, seq1LastBase = 0;
    int seq2FirstBase = -1, seq2LastBase = 0;
    for (int i = 0; i < alignmentLength; ++i)
    {
        char b1 = alignedSeq1[i];
        char b2 = alignedSeq2[i];

        if (b1 != '-')
        {
            if (seq1FirstBase == -1)
                seq1FirstBase = i;
            seq1LastBase = i;
        }
        if (b2 != '-')
        {
            if (seq2FirstBase == -1)
                seq2FirstBase = i;
            seq2LastBase = i;
        }
    }

    int alignmentStart = std::max(seq1FirstBase, seq2FirstBase);
    int alignmentEnd = std::min(seq1LastBase, seq2LastBase);
//    int alignedPartLength = alignmentEnd - alignmentStart + 1;
//    std::string seq1AlignedPart = alignedSeq1.substr(alignmentStart, alignedPartLength);
//    std::string seq2AlignedPart = alignedSeq2.substr(alignmentStart, alignedPartLength);

    int matches = 0;
    int mismatches = 0;
    int insertions = 0; //seq2 base not in seq1
    int deletions = 0; //seq1 base not in seq2

    for (int i = alignmentStart; i <= alignmentEnd; ++i)
    {
        char b1 = alignedSeq1[i];
        char b2 = alignedSeq2[i];

        if (b1 == '-')
            ++insertions;
        else if (b2 == '-')
            ++deletions;
        else if (b1 == b2)
            ++matches;
        else
            ++mismatches;
    }
//    std::cout << std::endl;

//    std::cout << seq1FirstBase << " " << seq1LastBase << std::endl;
//    std::cout << seq2FirstBase << " " << seq2LastBase << std::endl;
//    std::cout << seq1AlignedPart << std::endl;
//    std::cout << seq2AlignedPart << std::endl;

//    std::cout << "matches: " << matches << std::endl;
//    std::cout << "mismatches: " << mismatches << std::endl;
//    std::cout << "insertions: " << insertions << std::endl;
//    std::cout << "deletions: " << deletions << std::endl;


    std::string returnedString = "";
    returnedString += std::to_string(seq1FirstBase) + '\t';
    returnedString += std::to_string(seq1LastBase + 1) + '\t';
    returnedString += std::to_string(seq2FirstBase) + '\t';
    returnedString += std::to_string(seq2LastBase + 1) + '\t';
    returnedString += std::to_string(matches) + '\t';
    returnedString += std::to_string(mismatches) + '\t';
    returnedString += std::to_string(insertions) + '\t';
    returnedString += std::to_string(deletions) + '\t';

    std::cout << returnedString;

    return returnedString;
}


int main()
{
    semiGlobalAlign("AACATGGCTAAATTTCATGCCAGAATCAGTC",
                    "GCGCTGAACATGGCTACATTCATCCAGAA");

    return 0;
}


