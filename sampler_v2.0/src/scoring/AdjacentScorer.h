#ifndef _SCORING_ADJACENTSCORING_H
#define _SCORING_ADJACENTSCORING_H

#include "../Config.h"
#include "ScoringScheme.h"

class AdjacentScorer : public ScoringScheme
{
public:
	AdjacentScorer(WindowContainer * winContRef);
	double configAdjustedWeight(Config & config);

private:
	double regionPairScore(const Region & x, const Region & y);
	double bestPairingRegions(const vector<Region> & regions);

};


#endif
