
#include "commonkmerset.h"
#include "settings.h"


CommonKmerSet::CommonKmerSet(std::string & readName, std::string & refName, int readLength, int refLength,
                             float expectedSlope, KmerPositions * kmerPositions) :
    m_readName(readName), m_refName(refName),
    m_readLength(readLength), m_refLength(refLength),
    m_expectedSlope(expectedSlope), m_maxScore(0), m_maxScoreIndex(0)
{
    float rotationAngle = CommonKmer::getRotationAngle(m_expectedSlope);

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

    // Sort by rotated vertical position so lines should be roughly horizontal.
    std::sort(m_commonKmers.begin(), m_commonKmers.end(), [](const CommonKmer & a, const CommonKmer & b) {
        return a.m_rotatedVPosition < b.m_rotatedVPosition;   
    });

    scorePoints();
}


// Gives a score to each CommonKmer based on the density of CommonKmer points in the horizontal
// band (in the rotated frame).
void CommonKmerSet::scorePoints() {
    int halfBandPoints = COMMON_KMER_BAND_SIZE / 2;
    double halfBandThickness = COMMON_KMER_BAND_THICKNESS / 2.0;

    // There are four corners of the alignment rectangle which we also need to rotate.
    float rotationAngle = CommonKmer::getRotationAngle(m_expectedSlope);
    CommonKmer c1(0, 0, rotationAngle);
    CommonKmer c2(0, m_refLength, rotationAngle);
    CommonKmer c3(m_readLength, m_refLength, rotationAngle);
    CommonKmer c4(m_readLength, 0, rotationAngle);
    float c1BandLength = getLineLength(c1.m_hPosition, c1.m_vPosition,
                                       m_expectedSlope, m_readLength, m_refLength);
    float c3BandLength = getLineLength(c3.m_hPosition, c3.m_vPosition,
                                       m_expectedSlope, m_readLength, m_refLength);

    // If the reference is longer than the read, the C3 corner will be higher. If the read is
    // longer than the reference, the C1 corner will be higher.
    float lowerCornerY = c1.m_rotatedVPosition;
    float higherCornerY = c3.m_rotatedVPosition;
    if (lowerCornerY > higherCornerY)
        std::swap(lowerCornerY, higherCornerY);

    // We will scale the scores relative to the expected k-mer density.
    float expectedDensity = 1.0 / pow(4.0f, KMER_SIZE);

    // Now we loop through the CommonKmer points, calculating their k-mer density (score) along the
    // way. The density is gotten from a window which spans either some number of steps away from
    // the current point or some distance, whichever is larger.
    int bandStartIndex = -1, bandEndIndex = -1;
    int commonKmerCount = m_commonKmers.size();
    for (int i = 0; i < commonKmerCount; ++i) {

        // Advance the window bounds, if appropriate.
        double currentY = getY(m_commonKmers, i, commonKmerCount, c2, c4);
        while (bandStartIndex < i - halfBandPoints &&
               getY(m_commonKmers, bandStartIndex, commonKmerCount, c2, c4) < currentY - halfBandThickness)
            ++bandStartIndex;
        while (bandEndIndex < commonKmerCount &&
               (bandEndIndex < i + halfBandPoints || getY(m_commonKmers, bandEndIndex, commonKmerCount, c2, c4) < currentY + halfBandThickness) )
            ++bandEndIndex;

        int thisBandCount = bandEndIndex - bandStartIndex;
        float bandStartY = getY(m_commonKmers, bandStartIndex, commonKmerCount, c2, c4);
        float bandEndY = getY(m_commonKmers, bandEndIndex, commonKmerCount, c2, c4);

        // Now we need to calculate the area of the band, which is a little complex because it is a
        // horizontal band through a rotated rectangle.
        float bandArea;

        // We'll need the starting point band length for all area possibilities, so we calculate
        // that now.
        float bandStartLength;
        if (bandStartIndex < 0)
            bandStartLength = 0.0;
        else
            bandStartLength = getLineLength(m_commonKmers[bandStartIndex].m_hPosition,
                                            m_commonKmers[bandStartIndex].m_vPosition,
                                            m_expectedSlope, m_readLength, m_refLength);

        // If both the start and end are in the middle of the rotated rectangle, then the area is a
        // parallelogram and calculating its area is easy.
        if (bandStartY >= lowerCornerY && bandEndY <= higherCornerY) {
            float parallelogramArea = (bandEndY - bandStartY) * bandStartLength;
            bandArea = parallelogramArea;
        }

        // Other cases are more complex, and we'll need the ending point band length too.
        else {
            float bandEndLength;
            if (bandEndIndex >= commonKmerCount)
                bandEndLength = 0.0;
            else
                bandEndLength = getLineLength(m_commonKmers[bandEndIndex].m_hPosition,
                                              m_commonKmers[bandEndIndex].m_vPosition,
                                              m_expectedSlope, m_readLength, m_refLength);

            // If both the start and end are in the bottom or top of the rotated rectangle, then the
            // area is a triangle/trapezoid.
            if ( (bandStartY <= lowerCornerY && bandEndY <= lowerCornerY) || (bandStartY >= higherCornerY && bandEndY >= higherCornerY) ) {
                float trapezoidArea = (bandEndY - bandStartY) * ((bandStartLength + bandEndLength) / 2.0);
                bandArea = trapezoidArea;
            }

            // If the start and end span a rectangle's corner, then the area is more complex. We need
            // to add both the parallelogram and trapezoid components.
            else if (bandStartY <= lowerCornerY && bandEndY >= lowerCornerY && bandEndY <= higherCornerY) {
                float trapezoidArea = (lowerCornerY - bandStartY) * ((bandStartLength + c1BandLength) / 2.0);
                float parallelogramArea = (bandEndY - lowerCornerY) * bandEndLength;
                bandArea = trapezoidArea + parallelogramArea;
            }
            else if (bandStartY >= lowerCornerY && bandStartY <= higherCornerY && bandEndY >= higherCornerY) {
                float trapezoidArea = (bandEndY - higherCornerY) * ((c3BandLength + bandEndLength) / 2.0);
                float parallelogramArea = (higherCornerY - bandStartY) * bandStartLength;
                bandArea = trapezoidArea + parallelogramArea;
            }

            // The final, most complex scenario is when the band start and end span both C1 and C3.
            // This is less common, as it would require either a very large band or a very sparse
            // set of CommonKmers.
            else {
                float trapezoidArea1 = (lowerCornerY - bandStartY) * ((bandStartLength + c1BandLength) / 2.0);
                float parallelogramArea = (higherCornerY - lowerCornerY) * c1BandLength;
                float trapezoidArea2 = (bandEndY - higherCornerY) * ((c3BandLength + bandEndLength) / 2.0);
                bandArea = trapezoidArea1 + parallelogramArea + trapezoidArea2;
            }
        }

        // Now that we have the band area, we can get the density of CommonKmers in the band. Also,
        // we'll scale this to the expected level of CommonKmers (given a random sequence).
        float kmerDensity = thisBandCount / bandArea;
        float score = (kmerDensity / expectedDensity) - 1.0;
        m_commonKmers[i].m_score = score;
    }
    findMaxScore();
}


void CommonKmerSet::findMaxScore() {
    m_maxScore = std::numeric_limits<float>::min();
    m_maxScoreIndex = -1;
    for (size_t i = 0; i < m_commonKmers.size(); ++i) {
        double score = m_commonKmers[i].m_score;
        if (score > m_maxScore) {
            m_maxScore = score;
            m_maxScoreIndex = i;
        }
    }
}


// This function builds an alignment line around the highest scoring point in the set. It removes
// the points from the set that go into the line.
AlignmentLine * CommonKmerSet::extractAlignmentLine() {

    // If a max score isn't ready, then we can't proceed.
    if (m_maxScoreIndex == -1)
        return 0;

    // Start the line with the max scoring point.
    int lineStartI = m_maxScoreIndex;
    int lineEndI = m_maxScoreIndex;
    AlignmentLine * line = new AlignmentLine(m_commonKmers[lineStartI],
                                             m_readName, m_refName, m_readLength, m_refLength);

    // Extend the line in both directions (always choosing the higher scoring side) until either
    // the score drops below a threshold or the line error gets too high.
    int commonKmerCount = m_commonKmers.size();
    while (true) {

        // We want to add a single point to the line, so we assess the two candidate points (one
        // above the current points and one below).
        int i1 = lineStartI - 1;
        int i2 = lineEndI + 1;
        bool i1Good = (i1 >= 0);
        bool i2Good = (i2 < commonKmerCount);
        float p1Score, p2Score;
        if (i1Good)
            p1Score = m_commonKmers[i1].m_score;
        else
            p1Score = std::numeric_limits<float>::min();
        if (i2Good)
            p2Score = m_commonKmers[i2].m_score;
        else
            p2Score = std::numeric_limits<float>::min();

        // If both points are bad (have left the vector bounds) then we're done.
        if (!i1Good && !i2Good)
            break;

        // Determine which of the two options is better.
        int bestI;
        float bestScore;
        if (p1Score >= p2Score) {
            bestI = i1;
            bestScore = p1Score;
        }
        else { // p2Score > p1Score
            bestI = i2;
            bestScore = p2Score;
        }

        // If the best score is below the threshold, we're done.
        if (bestScore < MIN_LINE_SCORE)
            break;

        if (bestI == i1)
            --lineStartI;
        else // bestI == i2
            ++lineEndI;

        // Add the point to the line, which will cause the line to update its linear regression.
        line->addPoint(m_commonKmers[bestI]);

        // If the line error has gotten too large, we're done. But we only check this after the
        // line has gotten a few points, otherwise it could fail here too early.
        int pointCount = lineEndI - lineStartI + 1;
        if (pointCount >= MIN_POINT_COUNT) {
            if (line->getRelativeLineError() > MAX_ALLOWED_LINE_RELATIVE_ERROR ||
                line->m_error > MAX_ALLOWED_LINE_ABSOLUTE_ERROR)
                break;
        }
    }

    // We should now have an alignment line with at least one (hopefully more) points. We remove
    // those points from the set and recalculate the set's scores.
    m_commonKmers.erase(m_commonKmers.begin() + lineStartI, m_commonKmers.begin() + lineEndI + 1);
    scorePoints(); // TO DO: THIS IS PRETTY INEFFICIENT, AS I'M RECALCULATING EVERY SCORE. I COULD INSTEAD JUST RECALCUATE THE SCORES FOR POTENTIALLY AFFECTED POINTS (THOSE NEAR THE LINE).

    return line;
}


// This function is used in the CommonKmer scoring loop. It uses an index to return a Y value
// (a.k.a. a rotated V) from the CommonKmer. But if the index is out of range, it gives the
// appropriate Y value from the rectangle corner instead.
float getY(std::vector<CommonKmer> & commonKmers, int i, int count, CommonKmer & c2, CommonKmer & c4) {
    if (i < 0)
       return c4.m_rotatedVPosition;
    else if (i >= count)
        return c2.m_rotatedVPosition;
    else
        return commonKmers[i].m_rotatedVPosition;
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

