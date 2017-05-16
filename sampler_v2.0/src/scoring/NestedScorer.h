#include <vector>
#include "../Config.h"
#include "ScoringScheme.h"

class NestedScorer : public ScoringScheme
{
public:
	NestedScorer(WindowContainer * winContRef);
	double configAdjustedWeight(Config & config);
	
private:
	double best_pairing(const std::vector<Region> & S, int i, int j);
	void read_joint_prob_file(int rna_num);

	unordered_map<string, double> pair_up_energy;
};
