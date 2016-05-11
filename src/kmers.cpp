
#include "kmers.h"




// Creates the CommonKmer object using the position in the two sequences.
// It then rotates the point relative to the origin by the given angle (in radians).
CommonKmer::CommonKmer(std::string sequence, int hPosition, int vPosition, double angle) :
    m_sequence(sequence),
    m_hPosition(hPosition),
    m_vPosition(vPosition),
    m_bandArea(0.0),
    m_score(1.0) {
    double s = sin(angle);
    double c = cos(angle);
    m_rotatedHPosition = (m_hPosition * c) - (m_vPosition * s);
    m_rotatedVPosition = (m_hPosition * s) + (m_vPosition * c);
}


// This is the destructor for KmerPositions. It cleans up all the KmerPosMaps which were allocated
// on the heap.
KmerPositions::~KmerPositions() {
    for (std::unordered_map<std::string, KmerPosMap *>::iterator i = m_kmerPositions.begin(); i != m_kmerPositions.end(); ++i)
        delete i->second;
}

// This function adds a sequence to the KmerPositions object. It creates a new KmerPosMap on the
// heap (will be deleted in destructor), fills it up and adds it to m_kmerPositions.
void KmerPositions::addPositions(char * nameC, char * sequenceC) {
    std::string name(nameC);
    std::string sequence(sequenceC);

    KmerPosMap * posMap = new KmerPosMap();
    int kCount = sequence.size() - KMER_SIZE + 1;
    for (int i = 0; i < kCount; ++i) {
        std::string kmer = sequence.substr(i, KMER_SIZE);
        if (posMap->find(kmer) == posMap->end())
            (*posMap)[kmer] = std::vector<int>();
        (*posMap)[kmer].push_back(i);
    }

    m_kmerPositions[name] = posMap;
}

void KmerPositions::deletePositions(char * nameC) {
    std::string name(nameC);
    KmerPosMap * kmerPosMap = getKmerPositions(name);
    if (kmerPosMap != 0) {
        m_kmerPositions.erase(name);
        delete kmerPosMap;
    }
}

// This function retrieves a KmerPosMap from the object using the name as a key. If the name isn't
// in the map, it returns 0.
KmerPosMap * KmerPositions::getKmerPositions(std::string & name) {
    if (m_kmerPositions.find(name) == m_kmerPositions.end())
        return 0;
    else
        return m_kmerPositions[name];
}



// This function returns a list of the k-mers common to the two sequences.
std::vector<CommonKmer> getCommonKmers(std::string & readName, std::string & refName,
                                       double expectedSlope, KmerPositions * kmerPositions) {
    std::vector<CommonKmer> commonKmers;
    double rotationAngle = CommonKmer::getRotationAngle(expectedSlope);

    KmerPosMap * readKmerPositions = kmerPositions->getKmerPositions(readName);
    KmerPosMap * refKmerPositions = kmerPositions->getKmerPositions(refName);

    KmerPosMap * smaller = readKmerPositions;
    KmerPosMap * larger = refKmerPositions;
    bool refKmersSmaller = false;
    if (smaller->size() > larger->size()) {
        std::swap(smaller, larger);
        refKmersSmaller = true;
    }

    for (KmerPosMap::iterator i = smaller->begin(); i != smaller->end(); ++i) {
        std::string kmer = i->first;
        KmerPosMap::iterator j = larger->find(kmer);
        if (j != larger->end()) {
            // If the code got here, then a common k-mer was found!
            std::vector<int> * readPositions = &(i->second);
            std::vector<int> * refPositions = &(j->second);
            if (refKmersSmaller)
                std::swap(readPositions, refPositions);

            for (size_t k = 0; k < readPositions->size(); ++k) {
                for (size_t l = 0; l < refPositions->size(); ++l)
                    commonKmers.push_back(CommonKmer(kmer, (*readPositions)[k], (*refPositions)[l], rotationAngle));
            }
        }
    }

    return commonKmers;
}



KmerPositions * newKmerPositions() {
    return new KmerPositions();
}

void addKmerPositions(KmerPositions * kmerPositions, char * name, char * sequence) {
    kmerPositions->addPositions(name, sequence);
}

void deleteKmerPositions(KmerPositions * kmerPositions, char * name) {
    kmerPositions->deletePositions(name);
}

void deleteAllKmerPositions(KmerPositions * kmerPositions) {
    delete kmerPositions;
}


