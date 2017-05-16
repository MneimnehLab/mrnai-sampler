#ifndef _SCORING_SCORINGSCHEME_H
#define _SCORING_SCORINGSCHEME_H


#include "../WindowContainer.h"
#include "../Config.h"


class ScoringScheme
{
public:
	ScoringScheme(WindowContainer * winContRef_) : 
		winCont(winContRef_)
	{
	};
	virtual double configAdjustedWeight(Config &) = 0;

protected:
	WindowContainer * winCont;
	vector<int> evenRNAs;
	vector<int> oddRNAs;
};

#endif