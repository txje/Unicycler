
#include "ref_seqs.h"

SeqMap * newRefSeqs() {
    return new SeqMap();
}

void addRefSeq(SeqMap * seqMap, char * nameC, char * sequenceC) {
    seqMap->emplace(nameC, sequenceC);
}

void deleteRefSeqs(SeqMap * seqMap) {
    delete seqMap;
}
