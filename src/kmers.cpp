
#include "kmers.h"




// Creates the CommonKmer object using the position in the two sequences.
// It then rotates the point relative to the origin by the given angle (in radians).
CommonKmer::CommonKmer(int hPosition, int vPosition, float angle) :
    m_hPosition(hPosition),
    m_vPosition(vPosition),
    m_score(1.0) {
    float s = sin(angle);
    float c = cos(angle);
    m_rotatedHPosition = (m_hPosition * c) - (m_vPosition * s);
    m_rotatedVPosition = (m_hPosition * s) + (m_vPosition * c);
}


// This is the destructor for KmerPositions. It cleans up all the KmerPosMaps which were allocated
// on the heap.
KmerPositions::~KmerPositions() {
    for (std::unordered_map<std::string, KmerPosMap *>::iterator i = m_kmerPositions.begin(); i != m_kmerPositions.end(); ++i)
        delete i->second;
}


// Returns a vector all of k-mer position names (should be exact the same as the the sequence names).
std::vector<std::string> KmerPositions::getAllNames() {
    std::vector<std::string> returnVector;
    for (std::unordered_map<std::string, KmerPosMap *>::iterator i = m_kmerPositions.begin(); i != m_kmerPositions.end(); ++i)
        returnVector.push_back(i->first);
    return returnVector;
}

// Returns the length of the sequence with the given name.
int KmerPositions::getLength(std::string & name) {
    if (m_sequences.find(name) == m_sequences.end())
        return 0;
    else
        return m_sequences[name].length();
}

// This function adds a sequence to the KmerPositions object. It creates a new KmerPosMap on the
// heap (will be deleted in destructor), fills it up and adds it to m_kmerPositions.
void KmerPositions::addPositions(std::string & name, std::string & sequence) {
    m_sequences[name] = sequence;

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

void KmerPositions::deletePositions(std::string & name) {
    if (m_sequences.find(name) != m_sequences.end())
        m_sequences.erase(name);
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


std::string * KmerPositions::getSequence(std::string & name) {
    if (m_kmerPositions.find(name) == m_kmerPositions.end())
        return 0;
    else
        return &(m_sequences[name]);

}



CommonKmerSet::CommonKmerSet(std::string & readName, std::string & refName,
                             int readLength, int refLength, int bandSize,
                             float expectedSlope, KmerPositions * kmerPositions) :
    m_readName(readName), m_refName(refName), m_expectedSlope(expectedSlope), m_maxScore(0)
{
    float rotationAngle = CommonKmer::getRotationAngle(expectedSlope);

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
                    m_commonKmers.push_back(CommonKmer((*readPositions)[k], (*refPositions)[l], rotationAngle));
            }
        }
    }

    // We scale the scores relative to the expected k-mer density.
    float expectedDensity = 1.0 / pow(4.0f, KMER_SIZE);

    // Sort by rotated vertical position so lines should be roughly horizontal.
    std::sort(m_commonKmers.begin(), m_commonKmers.end(), [](const CommonKmer & a, const CommonKmer & b) {
        return a.m_rotatedVPosition < b.m_rotatedVPosition;   
    });

    // Score each point based on the number of other points in its band.
    int halfBandSize = bandSize / 2;

    // There are four corners of the alignment rectangle which we also need to rotate.
    CommonKmer c1(0, 0, rotationAngle);
    CommonKmer c2(0, refLength, rotationAngle);
    CommonKmer c3(readLength, refLength, rotationAngle);
    CommonKmer c4(readLength, 0, rotationAngle);
    float c1Y = c1.m_rotatedVPosition;
    float c3Y = c3.m_rotatedVPosition;
    float c1BandLength = getLineLength(c1.m_hPosition, c1.m_vPosition,
                                       expectedSlope, readLength, refLength);
    float c3BandLength = getLineLength(c3.m_hPosition, c3.m_vPosition,
                                       expectedSlope, readLength, refLength);

    // Now we loop through the CommonKmer points, calculating their k-mer density (score) along the way.
    int commonKmerCount = m_commonKmers.size();
    for (int i = 0; i < commonKmerCount; ++i) {
        int bandStartIndex = std::max(i - halfBandSize, 0);
        int bandEndIndex = std::min(i + halfBandSize, commonKmerCount - 1);
        int thisBandSize = bandEndIndex - bandStartIndex;

        // Get the Y coordinates for the start and end of the band.
        float bandStartY, bandEndY;
        if (bandStartIndex < 0)
            bandStartY = c4.m_rotatedVPosition;
        else
            bandStartY = m_commonKmers[bandStartIndex].m_rotatedVPosition;
        if (bandEndIndex >= commonKmerCount)
            bandEndY = c2.m_rotatedVPosition;
        else
            bandEndY = m_commonKmers[bandEndIndex].m_rotatedVPosition;

        // Now we need to calculate the area of the band.
        float bandArea;

        // We'll need the starting point band length for all area possibilities, so we calculate
        // that now.
        float bandStartLength;
        if (bandStartIndex < 0)
            bandStartLength = 0.0;
        else
            bandStartLength = getLineLength(m_commonKmers[bandStartIndex].m_hPosition,
                                            m_commonKmers[bandStartIndex].m_vPosition,
                                            expectedSlope, readLength, refLength);

        // If both the start and end are in the middle of the rotated rectangle, then the area is a
        // parallelogram and calculating its area is easy.
        if (bandStartY >= c1Y && bandEndY <= c3Y)
            bandArea = (bandEndY - bandStartY) * bandStartLength;

        // Other cases are more complex, and we'll need the ending point band length too.
        else {
            float bandEndLength;
            if (bandEndIndex >= commonKmerCount)
                bandEndLength = 0.0;
            else
                bandEndLength = getLineLength(m_commonKmers[bandEndIndex].m_hPosition,
                                              m_commonKmers[bandEndIndex].m_vPosition,
                                              expectedSlope, readLength, refLength);

            // If both the start and end are in the bottom or top of the rotated rectangle, then the
            // area is a triangle/trapezoid.
            if ( (bandStartY <= c1Y && bandEndY <= c1Y) || (bandStartY >= c3Y && bandEndY >= c3Y) )
                bandArea = (bandEndY - bandStartY) * ((bandStartLength + bandEndLength) / 2.0);

            // If the start and end span a rectangle's corner, then the area is more complex. We need
            // to add both the parallelogram and trapezoid components.
            else if (bandStartY <= c1Y && bandEndY >= c1Y && bandEndY <= c3Y) {
                float trapezoidArea = (c1Y - bandStartY) * ((bandStartLength + c1BandLength) / 2.0);
                float parallelogramArea = (bandEndY - c1Y) * bandEndLength;
                bandArea = trapezoidArea + parallelogramArea;
            }
            else if (bandStartY >= c1Y && bandStartY <= c3Y && bandEndY >= c3Y) {
                float trapezoidArea = (bandEndY - c3Y) * ((c3BandLength + bandEndLength) / 2.0);
                float parallelogramArea = (c3Y - bandStartY) * bandStartLength;
                bandArea = trapezoidArea + parallelogramArea;
            }

            // The final, most complex scenario is when the band start and end span both C1 and C3.
            // This would be unusual, as it would require either a very large band or a very sparse
            // set of CommonKmers.
            else {
                float trapezoidArea1 = (c1Y - bandStartY) * ((bandStartLength + c1BandLength) / 2.0);
                float parallelogramArea = (c3Y - c1Y) * c1BandLength;
                float trapezoidArea2 = (bandEndY - c3Y) * ((c3BandLength + bandEndLength) / 2.0);
                bandArea = trapezoidArea1 + parallelogramArea + trapezoidArea2;
            }
        }

        // Now that we have the band area, we can get the density of CommonKmers in the band. Also,
        // we'll scale this to the expected level of CommonKmers (given a random sequence).
        float kmerDensity = thisBandSize / bandArea;
        float score = (kmerDensity / expectedDensity) - 1.0;
        m_commonKmers[i].m_score = score;
        m_maxScore = std::max(m_maxScore, score);
    }
}



// Given a point, a slope and rectangle bounds, this function returns the length of the line
// segment which goes through that point with that slope and is in those bounds.
float getLineLength(float x, float y, float slope, float xSize, float ySize) {
    float xStart, yStart, xEnd, yEnd;

    float yIntercept = y - (slope * x);
    if (yIntercept >= 0.0) {
        xStart = 0.0;
        yStart = yIntercept;
    }
    else {
        xStart = -yIntercept / slope;
        yStart = 0.0;
    }

    float yAtXEnd = (slope * xSize) + yIntercept;
    if (yAtXEnd <= ySize) {
        xEnd = xSize;
        yEnd = yAtXEnd;
    }
    else {
        xEnd = (ySize - yIntercept) / slope;
        yEnd = ySize;
    }

    float xLength = xEnd - xStart;
    float yLength = yEnd - yStart;
    return sqrt((xLength * xLength) + (yLength * yLength));
}


KmerPositions * newKmerPositions() {
    return new KmerPositions();
}

void addKmerPositions(KmerPositions * kmerPositions, char * nameC, char * sequenceC) {
    std::string name(nameC);
    std::string sequence(sequenceC);
    kmerPositions->addPositions(name, sequence);
}

// void deleteKmerPositions(KmerPositions * kmerPositions, char * name) {
//     kmerPositions->deletePositions(name);
// }

void deleteAllKmerPositions(KmerPositions * kmerPositions) {
    delete kmerPositions;
}


