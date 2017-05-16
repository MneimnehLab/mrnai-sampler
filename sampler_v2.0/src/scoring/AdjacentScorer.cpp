#include <iostream>
#include <vector>
#include "AdjacentScorer.h"
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


AdjacentScorer::AdjacentScorer(WindowContainer * winCont) : 
    ScoringScheme(winCont)
{    
}

double AdjacentScorer::configAdjustedWeight(Config & vec)
{
    try
    {
    vec.sort();
    if(DEBUG_ADJ_SCORE) cout << "Finding score of " << vec << endl;
    

if(DEBUG_ADJ_SCORE) cout << "numOfLevels = " << winCont->numOfLevels << endl;
if(DEBUG_ADJ_SCORE) cout << "numOfRNAs = " << winCont->numOfRNAs << endl;
    // vector< vector<Region> > regionsByLevel(numOfLevels);
    vector< vector<Region> > regionsByLevel(winCont->numOfRNAs);
    // for(int i=0; i<numOfLevels; i++)
    for(int i=0; i<winCont->numOfRNAs; i++)
    {
        Region r(-1,-1,-1);
        regionsByLevel[i].push_back(r);
    }
    // Make region vectors:
    // vector<Region> topRegions(1), botRegions(1);
    
    double totalScore = 0;

    for(Window win : vec)
    {
        
        // topRegions.push_back( Region (win.l,   win.i, win.w1) );
        // botRegions.push_back( Region (win.l+1, win.j, win.w2) );

        int odd  = win.l2;
        int even = win.l1;

        regionsByLevel[even].push_back( Region (even,   win.i, win.w1) );
        regionsByLevel[odd].push_back( Region (odd, win.j, win.w2) );

        // Compute add the interE for each window
        if(DEBUG_ADJ_SCORE) cout << win << "'s inter = " << winCont->interEs[win.toString()] << endl;
        totalScore += winCont->interEs[win.toString()];

    }
    if(DEBUG_ADJ_SCORE) cout << "Interaction scores = " << totalScore << endl;

    // for(int i=0; i<numOfLevels; i++)
    for(int i=0; i<winCont->numOfRNAs; i++)
    {
        sort(regionsByLevel[i].begin(), regionsByLevel[i].end(), WindowContainer::regionComparator);
        double rScore = bestPairingRegions(regionsByLevel[i]);  
        if(DEBUG_ADJ_SCORE) cout << "region score of level " << i << "  =  " << rScore << endl; 
        totalScore += rScore;
    }
    // double t = bestPairingRegions(topRegions);
    // double b = bestPairingRegions(botRegions);
    // cout << "t = " << t << endl;
    // cout << "b = " << b << endl;
    
    // totalScore += t + b;
    if(DEBUG_ADJ_SCORE) cout << "totalScore of " << vec << " = " << totalScore << endl << endl << endl; 



    return totalScore;
    }
    catch (const std::exception& e) { // caught by reference to base
        std::cerr << " !! a $tandard exception was caught, with message '"
                  << e.what() << "'\n";
                  throw e;
    }
}

double AdjacentScorer::regionPairScore(const Region & x, const Region & y) 
{
    const double MAX_VAL = 9999999.0;
    // Compute up of union (assume x < y):
    // int gap = y.i - y.w - x.i;
    int unionWidth = y.i-(x.i-x.w);
    // return MAX_VAL;
    
    // if(gap > 6 || unionWidth > 25)
    //     return MAX_VAL;

    
    Region unionRegion(x.l, y.i, unionWidth);
    if(DEBUG_ADJ_SCORE)cout << "union region = " << unionRegion << endl;

    // Check if we have score of union
    unordered_map<string, double>::const_iterator got = winCont->regionUpEs.find (unionRegion.toString());
    if ( got == winCont->regionUpEs.end())
    {
        if(DEBUG_ADJ_SCORE)std::cout << "union region " << unionRegion << " not found" << endl;
        return MAX_VAL;
    }
    else
    {
        // return MAX_VAL;
        if(DEBUG_ADJ_SCORE)cout << "found union region " << unionRegion << endl;
        return winCont->regionUpEs[unionRegion.toString()];
    }


}

// For now don't worry about returning actual pairing, just give the score
double AdjacentScorer::bestPairingRegions(const vector<Region> & regions) 
{
    if(DEBUG_ADJ_SCORE) {cout << "in bestPairingRegions, list of regions = " << endl;
        for(auto r:regions)
            cout << r << "  ";
        cout << endl;}

    vector<double> M(regions.size());
    M[0] = 0.0;
    M[1] = winCont->regionUpEs[regions[1].toString()];

    for(int i=2; i<regions.size(); i++)
    {
        if(DEBUG_ADJ_SCORE) cout << "At region " << regions[i] << endl;
        double singletonScore = M[i-1] + winCont->regionUpEs[regions[i].toString()];
        if(DEBUG_ADJ_SCORE) cout << "calling pair func" << endl;
        double pairScore =      M[i-2] + regionPairScore(regions[i-1], regions[i]);
        if(DEBUG_ADJ_SCORE) cout << i << ". singletonScore = " << singletonScore <<",  b/c M[i-1] = " << M[i-1] << " and rs = " << winCont->regionUpEs[regions[i].toString()]  << endl;
        if(DEBUG_ADJ_SCORE) cout << i << ". pairScore = " << pairScore << endl;

        M[i] = MIN(singletonScore, pairScore);
        if(DEBUG_ADJ_SCORE) cout << "M[" << i << "] = " << M[i] << endl;
        if(DEBUG_ADJ_SCORE) cout << endl;
    }
    
    if(DEBUG_ADJ_SCORE) {for(auto x : M)
            cout << "-> " << x << endl;
        cout << "up score = " << M.back() << endl;
        cout << endl;}

    return M.back();
}


