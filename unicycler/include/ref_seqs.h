
#ifndef REF_SEQS_H
#define REF_SEQS_H

#include <string>
#include <unordered_map>

typedef std::unordered_map<std::string, std::string> SeqMap;

// Functions that are called by the Python script must have C linkage, not C++ linkage.
extern "C" {
    SeqMap * newRefSeqs();
    void addRefSeq(SeqMap * seqMap, char * nameC, char * sequenceC);
    void deleteRefSeqs(SeqMap * seqMap);
}

#endif // REF_SEQS_H
