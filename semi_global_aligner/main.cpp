#include <chrono>
#include <seqan/basic.h>
#include <seqan/align.h>
#include <seqan/sequence.h>
#include <seqan/stream.h>  // for printint strings

using namespace seqan;
using namespace std::chrono;


int main()
{
    typedef Dna5String TSequence;
    typedef Align<Dna5String, ArrayGaps> TAlign;
    typedef Row<TAlign>::Type TRow;
    typedef Iterator<TRow>::Type TRowIterator;

    TSequence seq1 = "CGGAGAAATGCCGCGGCTCACGCAGTCTGACTGCCTGAGTAAGCGCATTTCCAGTTGTCACAAATCATAATGCGGGCCTCCGCTAATGTAAGATGGGG";
    TSequence seq2 = "GCTCACGCAGTCTGACTGCCTGAGTAAGCGCATTTCCAGTTGTCACAAATCATAATGCGGGCCTCCGCTAATGTAAGATGGGGGCCGATCTGGAGCTGCGGTGGGT";

    TAlign align;
    resize(rows(align), 2);
    assignSource(row(align, 0), seq1);
    assignSource(row(align, 1), seq2);

    Score<int, Simple> scoringScheme(1, -1, -1);
    AlignConfig<true, true, true, true> alignConfig;

    long long ms_before = duration_cast< milliseconds >(system_clock::now().time_since_epoch()).count();

    int score = globalAlignment(align, scoringScheme, alignConfig);

    long long ms_after = duration_cast< milliseconds >(system_clock::now().time_since_epoch()).count();
    long long elapsedTime = ms_after - ms_before;

    int alignmentLength = _max(length(row(align, 0)), length(row(align, 1)));

    std::cout << "seq1 size: " << length(seq1) << std::endl;
    std::cout << "seq2 size: " << length(seq2) << std::endl;
    std::cout << "Score: " << score << std::endl;
    std::cout << "Milliseconds: " << elapsedTime << std::endl;
    std::cout << align << std::endl;

    int seq1FirstBase = 0, seq1LastBase = 0, seq2FirstBase = 0, seq2LastBase = 0;
    for (unsigned i = 0; i < length(rows(align)); ++i)
    {
        TRowIterator it = iter(row(align, i), 0);
        TRowIterator itEnd = iter(row(align, i), alignmentLength);
        unsigned pos = 0;
        int firstBase = 0;
        int lastBase = 0;
        bool seenBase = false;
        while (it != itEnd)
        {
            if (!isGap(it))
            {
                if (!seenBase)
                {
                    firstBase = pos;
                    seenBase = true;
                }
                lastBase = pos;
            }
            ++it;
            ++pos;
        }

        if (i == 0)
        {
            seq1FirstBase = firstBase;
            seq1LastBase = lastBase;
        }
        else
        {
            seq2FirstBase = firstBase;
            seq2LastBase = lastBase;
        }
    }

    int seq1StartGap = seq1FirstBase;
    int seq2StartGap = seq2FirstBase;
    int seq1EndGap = alignmentLength - seq1LastBase - 1;
    int seq2EndGap = alignmentLength - seq2LastBase - 1;

    std::cout << "seq1StartGap: " << seq1StartGap << std::endl;
    std::cout << "seq1EndGap:   " << seq1EndGap << std::endl;
    std::cout << "seq2StartGap: " << seq2StartGap << std::endl;
    std::cout << "seq2EndGap:   " << seq2EndGap << std::endl;

    return 0;
}
