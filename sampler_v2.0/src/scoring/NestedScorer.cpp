#include <algorithm>
#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include "NestedScorer.h"
#include "../WindowContainer.h"
#include "../Region.h"

#define MIN(x, y) (((x) < (y)) ? (x) : (y))
#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#define INFTY 999999.0

#define DEBUG 0
#define DEBUG_ADJ_SCORE 0

// #define SCORING_MATCHING 1

#define INDEX(i, j, n) ((i)*(n)+(j))

using std::vector;
using std::string;
using std::cout;
using std::endl;
using std::pair;


NestedScorer::NestedScorer(WindowContainer * winCont) : 
    ScoringScheme(winCont)
{    
	// Read joint probability/energy files:
    // for now, fix file names at joint_prob_rna_i...
    // for now, only for 2 rnas
    read_joint_prob_file(0);
    read_joint_prob_file(1);
}


double NestedScorer::configAdjustedWeight(Config & vec)
{
    try
    {
    vec.sort();
    if(DEBUG_ADJ_SCORE) cout << "Finding score of " << vec << endl;
    

if(DEBUG_ADJ_SCORE) cout << "numOfLevels = " << winCont->numOfLevels << endl;
if(DEBUG_ADJ_SCORE) cout << "numOfRNAs = " << winCont->numOfRNAs << endl;
    
    vector< vector<Region> > regionsByLevel(winCont->numOfRNAs);
    
    double totalScore = 0;

    for(Window win : vec)
    {
        
        int odd  = win.l2;
        int even = win.l1;

        regionsByLevel[even].push_back( Region (even, win.i, win.w1) );
        regionsByLevel[odd] .push_back( Region (odd,  win.j, win.w2) );

        // Compute add the interE for each window
        if(DEBUG_ADJ_SCORE) cout << win << "'s inter = " << winCont->interEs[win.toString()] << endl;
        totalScore += winCont->interEs[win.toString()];

    }
    if(DEBUG_ADJ_SCORE) cout << "Interaction scores = " << totalScore << endl;

    // for(int i=0; i<numOfLevels; i++)
    for(int i=0; i<winCont->numOfRNAs; i++)
    {
        sort(regionsByLevel[i].begin(), regionsByLevel[i].end(), WindowContainer::regionComparator);
        double rScore = best_pairing(regionsByLevel[i], 0, regionsByLevel[i].size()-1);  
        if(DEBUG_ADJ_SCORE) cout << "region score of level " << i << "  =  " << rScore << endl; 
        totalScore += rScore;
    }
    
    if(DEBUG_ADJ_SCORE) cout << "totalScore of " << vec << " = " << totalScore << endl << endl << endl; 



    return totalScore;
    }
    catch (const std::exception& e) { // caught by reference to base
        std::cerr << " !! a $tandard exception was caught, with message '"
                  << e.what() << "'\n";
                  throw e;
    }
}




double NestedScorer::best_pairing(const vector<Region> & S, int i, int j)
{
    const int MAX_VAL = 999999;
    if(DEBUG_ADJ_SCORE) cout << "start = " << i << ",  end = " << j << endl;

    vector<Region>::const_iterator first = S.begin() + i;
    vector<Region>::const_iterator last = S.begin() + j + 1; // + 1 because this is inclusive
    vector<Region> subVec(first, last);

    if(DEBUG_ADJ_SCORE) {cout << "in new_regionScore, list subVec = " << endl;
        for(auto r:subVec)
            cout << r << "  ";
        cout << endl;}
    


    if(j < 0 or i >= S.size())
    {
        return 0;
    }
    else if(i == j)
    {
        Region r = S[i];
        unordered_map<string, double>::const_iterator got = winCont->regionUpEs.find (r.toString());
        if ( got == winCont->regionUpEs.end())
        {
            return MAX_VAL;
        }
        else
        {
            return got->second;
        }
    }
    else if(i > j)
    {
        return 0;
    }

    Region r1 = S[i];
    Region r2 = S[j];

    vector<double> vals;
    unordered_map<string, double>::const_iterator got0 = pair_up_energy.find (r1.toString() + r2.toString());
    unordered_map<string, double>::const_iterator got1 = winCont->regionUpEs.find (r1.toString());
    unordered_map<string, double>::const_iterator got2 = winCont->regionUpEs.find (r2.toString());

    double val1, val2, val0;
    // x1 = opt_pairing_bt(S, i+1, j-1, memo, depth+1)
    // v1 = x1[0] + joint_data[(a,b,c,d)]
    if (got0 == pair_up_energy.end())
    {
    	if(DEBUG_ADJ_SCORE)  cout << ">>>>> didn't find pair for " << r1 << " and " << r2 << endl;
        val0 = MAX_VAL;
    }
    else
    {
    	// val0 = MAX_VAL;
        val0 = got0->second;
		if(DEBUG_ADJ_SCORE) cout << ">>>>> found pair for " << r1 << " and " << r2 << " with energy = " << val0 << endl;
	}


    val0 += best_pairing(S, i+1, j-1);

    if (got1 == winCont->regionUpEs.end())
        val1 = MAX_VAL;
    else
        val1 = got1->second;
    val1 += best_pairing(S, i+1, j);

    if (got2 == winCont->regionUpEs.end())
        val2 = MAX_VAL;
    else
        val2 = got2->second;
    val2 += best_pairing(S, i, j-1);

    vals.push_back(val0);
    vals.push_back(val1);
    vals.push_back(val2);

    for(int k=i+1; k<j-1; k++)
    {
        double val = best_pairing(S, i, k) + best_pairing(S, k+1, j);
        vals.push_back(val);
    }

    if(DEBUG_ADJ_SCORE) {
    cout << "*********" << endl;
    for(auto e: vals)
        cout << "-> " << e << endl;
    cout << "*********" << endl;
    }

    return *(std::min_element(std::begin(vals), std::end(vals)));
   
}


void NestedScorer::read_joint_prob_file(int rna_num)
{   
    std::stringstream fname_stream;
    fname_stream << "data/joint_human/joint_prob_rna_" << rna_num;
    string fname = fname_stream.str();
    std::ifstream joint_prob_file(fname);

    if (joint_prob_file.fail())
    {
        throw std::runtime_error("couldn't open file " + fname + "\n"); 
    }


    int i1, j1, i2, j2;
    double energy;
    while(joint_prob_file >> i1 >> j1 >> i2 >> j2 >> energy)
    {
        // data in joint prob file is zero indexed...
		i1++; j1++; i2++; j2++;

        int w1 = j1 - i1;
        int w2 = j2 - i2;
        Region r1(rna_num, j1+1, w1);
        Region r2(rna_num, j2+1, w2);

        string key1 = r1.toString() + r2.toString();
        string key2 = r2.toString() + r1.toString();
        pair_up_energy[key1] = energy;
        pair_up_energy[key2] = energy;
    }
}
