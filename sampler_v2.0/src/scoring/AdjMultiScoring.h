#include "../Config.h"
#include "ScoringScheme.h"

class AdjacentScorer : public ScoringScheme
{
public:
	double configAdjustedWeight(Config & config);

private:

};


/*
double WindowContainer::new_unionScore(const vector<Region> & allRegions)
{
    const double MAX_VAL = 9999999.0;

    if(allRegions.size() == 0)
        return MAX_VAL;

    // Compute up of union (assume x < y):
    const Region & x = allRegions.front();
    const Region & y = allRegions.back();

    int gap = y.i - y.w - x.i;
    int unionWidth = y.i-(x.i-x.w);

    // if(gap > 6 || unionWidth > 25)
    //     return MAX_VAL;

    
    Region unionRegion(x.l, y.i, unionWidth);
    // cout << "union region = " << unionRegion << endl;

    // Check if we have score of union
    unordered_map<string, double>::const_iterator got = regionUpEs.find (unionRegion.toString());
    if ( got == regionUpEs.end())
    {
        // std::cout << "union region " << unionRegion << " not found" << endl;
        return MAX_VAL;
    }
    else
    {
        // return MAX_VAL;
        // cout << "found union region " << unionRegion << endl;
        return regionUpEs[unionRegion.toString()];
    }
} 


// For now don't worry about returning actual pairing, just give the score
double WindowContainer::new_regionScore(const vector<Region> & regions, int start, int end) 
{
    // cout << "start = " << start << ",  end = " << end << endl;

    vector<Region>::const_iterator first = regions.begin() + start;
    vector<Region>::const_iterator last = regions.begin() + end + 1; // + 1 because this is inclusive
    vector<Region> subVec(first, last);

    // if(DEBUG_ADJ_SCORE) {cout << "in new_regionScore, list subVec = " << endl;
    //     for(auto r:subVec)
    //         cout << r << "  ";
    //     cout << endl;}

    

    double unionScore = new_unionScore(subVec);
    double minScore = unionScore;

    if(DEBUG_ADJ_SCORE) cout << "Min for Union " << start << " - " << end << "  :  " << minScore << endl;



    for(int k=start; k<end; k++)
    {
        double left = new_regionScore(regions, start, k);
        double right = new_regionScore(regions, k+1, end);
        double total = left+right;

        if(total < minScore)
            minScore = total;
    }

    if(DEBUG_ADJ_SCORE) cout << "Min for range " << start << " - " << end << "  :  " << minScore << endl;

    return minScore;
}


double WindowContainer::configAdjustedWeight(Config & vec)
{
    try
    {
    vec.sort();
    if(DEBUG_ADJ_SCORE) cout << "Finding score of " << vec << endl;
    

if(DEBUG_ADJ_SCORE) cout << "numOfLevels = " << numOfLevels << endl;
if(DEBUG_ADJ_SCORE) cout << "numOfRNAs = " << numOfRNAs << endl;
    // vector< vector<Region> > regionsByLevel(numOfLevels);
    vector< vector<Region> > regionsByLevel(numOfRNAs);
    // for(int i=0; i<numOfLevels; i++)
    for(int i=0; i<numOfRNAs; i++)
    {
        // Region r(-1,-1,-1);
        // regionsByLevel[i].push_back(r);
    }
    // Make region vectors:
    // vector<Region> topRegions(1), botRegions(1);
    
    double totalScore = 0;

    for(Window win : vec)
    {
        
        // topRegions.push_back( Region (win.l,   win.i, win.w1) );
        // botRegions.push_back( Region (win.l+1, win.j, win.w2) );

        int odd  = win.l % numOdd;
        int even = win.l / numOdd;

        regionsByLevel[even].push_back( Region (even,   win.i, win.w1) );
        regionsByLevel[odd+numEven].push_back( Region (odd+numEven, win.j, win.w2) );

        // Compute add the interE for each window
        if(DEBUG_ADJ_SCORE) cout << win << "'s inter = " << interEs[win.toString()] << endl;
        totalScore += interEs[win.toString()];

    }
    if(DEBUG_ADJ_SCORE) cout << "Interaction scores = " << totalScore << endl;

    // for(int i=0; i<numOfLevels; i++)
    for(int i=0; i<numOfRNAs; i++)
    {
        sort(regionsByLevel[i].begin(), regionsByLevel[i].end(), regionComparator);
        double rScore = new_regionScore(regionsByLevel[i], 0, regionsByLevel[i].size()-1);  
        if(DEBUG_ADJ_SCORE) cout << "region score of level " << i << "  =  " << rScore << endl; 
        totalScore += rScore;

        if(DEBUG_ADJ_SCORE) cout << endl;
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
*/