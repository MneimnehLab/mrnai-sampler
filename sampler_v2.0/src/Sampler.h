#ifndef _SAMPLER_H_
#define _SAMPLER_H_

#include <iostream>
#include <vector>
#include <string>
#include <tuple>
#include <unordered_map>
#include "Window.h"
#include "WindowContainer.h"
#include "Config.h"
#include "scoring/ScoringScheme.h"



using std::vector;
using std::string;

//typedef unordered_map<string, double>;

class Sampler
{
public:
	struct Params
	{
		bool 	equalWins;
		int 	maxGoodWins;
		int 	maxWinSize;
		double 	threshold;
		int 	intervalSteps;
		int 	burnSteps;
		int 	randStartSteps;
		int 	numEven;
		int 	numOdd;
		vector<int> even_levels;
		vector<int> odd_levels;

		enum AType { FAST, SLOW, SLOW_HALF, INTERVAL };
		AType atype;
	};
private:
	WindowContainer * winContainer;
	unordered_map<string, double> weights;

	int numOfRNAs;
	int curr_level_pair_id;
	int curr_level_even;
	int curr_level_odd;
	int intervalSteps;
	int maxWinSize;
	int randStartSteps;

	

	Config lastState;
	double lastStateProb;

	ScoringScheme * scorer;


public:
	vector<std::pair<int, int>> level_id_to_pair;
	
	Sampler(int is, int bs, string dataFilename);
	// Sampler(int is, int bs, int maxWinSize);
	Sampler(const Sampler::Params & param);
	~Sampler();

	const vector<Window> getAddable(Config s) const;
	
	std::tuple<Config, double, double> sampleOneNeighbor(Config S);
	int getNeighborhoodCountAndMatch(Config S, Config compare = Config());

	std::tuple<Config, double, double> getAllNbrsExp(Config S, int depth);
	std::tuple<Config, double> makeRandomState();
	void getNext();
	Config getSample();
	void setLevel(int l);
	WindowContainer* getWinContainer();

	double getWeight(Config c);
	int getNumOfRNAs();

	int accept_count;
	int reject_count;

	static const int MAX_DEPTH_NBRS = 2;	//Should be at least 1, to allow going to first neighborhood

	int burnSteps;


	void read_joint_prob_file(int rna_num);

	
};


#endif