#include <algorithm>    // std::sort
#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <fstream>
#include <unordered_map>
#include <cstdlib>
#include <cmath>
#include <tuple>
#include <chrono>
#include <utility>
#include "Window.h"
#include "WindowContainer.h"
#include "Sampler.h"
#include "scoring/AdjacentScorer.h"
#include "scoring/NestedScorer.h"

#define DEBUG_MODE  0
#define VERBOSE 0

#define MIN(x, y) (((x) < (y)) ? (x) : (y))
#define MAX(x, y) (((x) > (y)) ? (x) : (y))

#define LOG(args) {if(VERBOSE) std::cout << args;}

extern bool cmdOptionExists(char** begin, char** end, const std::string& option);




Sampler::Sampler(const Sampler::Params & params)
{
    cout << "X" << endl;
    intervalSteps   = params.intervalSteps;
    burnSteps       = params.burnSteps;
    maxWinSize      = params.maxWinSize;
    randStartSteps  = params.randStartSteps;

    int numEven     = params.numEven;
    int numOdd      = params.numOdd;
    bool equalWins  = params.equalWins;
    double threshold = params.threshold;

    curr_level_pair_id = 0;
    curr_level_even = params.even_levels.at(0);
    curr_level_odd  = params.odd_levels.at(0);

    numOfRNAs = numEven + numOdd;


    cout << "numOfRNAs = " << numOfRNAs << endl;
    level_pair_to_id.resize(numOfRNAs+1);
    for(auto& e : level_pair_to_id)
        e.resize(numOfRNAs+1);


    int id = 0;
    for(int even : params.even_levels)
        for(int odd : params.odd_levels)
        {
            cout << "Even = " << even << " , odd = " << odd << endl;
            cout << "id = " << id << endl;
            level_pair_to_id[even][odd] = id;
            level_id_to_pair.push_back(std::make_pair(even, odd));
            id++;
        }


    cerr << "* Num Even : " << numEven << endl;
    cerr << "* Num Odd : " << numOdd << endl;
    
    // winContainer = new WindowContainer(numEven, numOdd);
    winContainer = new WindowContainer(params.even_levels, params.odd_levels);
    
    vector<std::tuple<Window, double, double> > tempContainerToSort;

    int maxLevel = 0;
    for (string line; getline(cin, line);) 
    {

        if(DEBUG_MODE) cout << ">" << line << endl;
        
        /*
        Window win = winContainer->addWindow(line);
        //cout << win.toStringW() << endl;

        maxLevel = win.l > maxLevel ? win.l : maxLevel;
        */

        // Parse line here:
        stringstream ss;
        ss << line;

        // Read indices
        // level is essentially "pair id"
        int level1, level2, i_, i, j_, j;
        

        // Added on 1/20/2016: fpg file now includes level # in first column
        ss >> level1 >> level2 >> i_ >> i >> j_ >> j;

        // note, level is essentially equivalent to id = even * NUM_ODD + odd
        // int odd = level % numOdd;
        // int even = level / numOdd;

        // Only for use when storing in N-dimensional matrix
        // int eIndex = even;
        // int oIndex = numEven + odd;



        // Read energies
        // double intlnZ, totlUpE, delta, p1, p2, upe1, upe2, E, Edelta;
        double intlnZ, totlUpE, delta, upe1, upe2;
        ss >> delta >> totlUpE >> intlnZ >> upe1 >> upe2;;

        // cout << "level1 = " << level1 << ",  level2 = " << level2 << endl;
        Window win(level1, level2, i, j, i-i_, j-j_, delta, intlnZ);


        // if(DEBUG_MODE) cout << "win = " << win << endl;
        // cout << "** win = " << win << endl;
        // cout << upe1 << "    " << upe2 << endl;
        winContainer->addEvenRegion(win, upe1);
        winContainer->addOddRegion(win, upe2);



        // if(delta > 0)
        if(intlnZ > 0)
            continue;

        if(equalWins && win.w1 != win.w2)
            continue;

        if(win.w1 > maxWinSize or win.w2 > maxWinSize)
        // if(win.w1 > maxWinSize)
            continue;

        // Now create region information

        //cout << win.toStringW() << endl;
        
        // winContainer->addWindow(win, upe1, upe2);
        // maxLevel = win.l > maxLevel ? win.l : maxLevel;

        // if(win.weight/(win.w1 + win.w2) > threshold)
  //        continue;


        tempContainerToSort.push_back(std::make_tuple(win, upe1, upe2));

    }

    
    // Sort container based on energy: indx 0 => min energy
    std::sort(tempContainerToSort.begin(), tempContainerToSort.end(), 
        [](tuple<Window, double, double> a, tuple<Window, double, double> b) 
    {
        return std::get<0>(a).intrEnergy < std::get<0>(b).intrEnergy;   
    });

    // for(auto elem: tempContainerToSort)
    // {
 //     Window win(0,0,0,0,0);
 //     double upe1, upe2;
 //     std::tie(win, upe1, upe2) = elem;

    //  cout << win << "\t" << win.weight << "/" << (win.w1 + win.w2) 
    //       << "\t" << win.weight/(win.w1 + win.w2) 
    //       << endl;
    // }







    // Now add good windows to winContainer
    int maxGoodWindows = params.maxGoodWins;

    for(int i=0; i<maxGoodWindows && i<tempContainerToSort.size(); i++)
    {
        // tuple<Window, double, double> winpk = tempContainerToSort[i];

        Window win(0,0,0,0,0,0);
        double upe1, upe2;
        std::tie(win, upe1, upe2) = tempContainerToSort[i];

        // if(win.weight/(win.w1 + win.w2) > -0.2)
        //  continue;

        if(DEBUG_MODE) cout << "Adding window " << win << "   with intr = " << win.intrEnergy 
        //   << "  and ratio = " << win.weight/(win.w1 + win.w2) 
             << endl;
        // winContainer->addWindow(win, upe1, upe2);

        win.id = i;
        winContainer->addWindow(win);


        maxLevel = MAX( MAX(win.l1, win.l2), maxLevel);
    }


    


    
    if(DEBUG_MODE) cout << "Number of RNAS: " << numOfRNAs << endl;

 //    if(params.atype == Params::FAST)
 //     winContainer->makeOverlaps();
 //    else if(params.atype == Params::SLOW)
 //     winContainer->makeOverlaps_Slow();
 //    else if(params.atype == Params::SLOW_HALF)
    //  winContainer->makeOverlaps_Slow_Half();
    // else
    // {
    //  cout << "intervals" << endl;
    //  winContainer->makeOverlaps_with_intervals();
    // }

    winContainer->makeOverlaps_with_intervals();    

    scorer = new AdjacentScorer(winContainer);
    // scorer = new NestedScorer(winContainer);

    

    cerr << "**All wins size : "      << winContainer->allWins.size() << endl;
    // cerr << "** with intlnZ <= 0"      << endl;
    cerr << "**** Rand start steps: " << randStartSteps << endl;
    cerr << "**** Burn steps: "       << burnSteps << endl;
    cerr << "**** Intervals steps: "  << intervalSteps << endl;
    cerr << "**** Max win size: "     << maxWinSize << endl;
    cerr << "**** Equal win size: "   << equalWins << endl;
    
    accept_count = 0;
    reject_count = 0;




    // Config S1 = Config::fromString("[(0, 1, 4, 4, 3, 3), (0, 1, 23, 29, 15, 17), (0, 1, 58, 38, 8, 8)]");
    // cout << "ans: " << scorer->configAdjustedWeight(S1) << endl << endl;

    // Config S2 = Config::fromString("[(0, 1, 4, 4, 3, 3), (0, 1, 14, 15, 3, 3), (0, 1, 58, 38, 8, 8)]");
    // cout << "ans: " << scorer->configAdjustedWeight(S2) << endl << endl;

    // Config S3 = Config::fromString("[(0, 1, 4, 4, 3, 3), (0, 1, 14, 15, 3, 3), (0, 1, 17, 20, 2, 2), (0, 1, 58, 38, 8, 8)]");
    // cout << "ans: " << scorer->configAdjustedWeight(S3) << endl << endl;

    // Config S4 = Config::fromString("[(0, 1, 4, 4, 3, 3), (0, 1, 17, 20, 6, 8), (0, 1, 58, 38, 8, 8)]");
    // cout << "ans: " << scorer->configAdjustedWeight(S4) << endl << endl;


    // exit(0);

    // Config S2 = Config::fromString("[(1, 0, 20, 24, 9, 9), (1, 2, 10, 10, 9, 9), (3, 0, 9, 8, 7, 6), (3, 0, 14, 14, 3, 3), (3, 2, 19, 15, 4, 4), (3, 2, 32, 30, 5, 5)]");
    // cout << "ans: " << winContainer->configAdjustedWeight(S2) << endl << endl;

    // Config S3 = Config::fromString("[(0, 4, 4, 3, 3), (0, 17, 20, 6, 8), (0, 58, 38, 8, 8)]");
    // cout << "ans: " << winContainer->configAdjustedWeight(S3) << endl << endl;


    // Config S1 = Config::fromString("[(0, 30, 9, 5, 5), (0, 58, 38, 9, 9)]");
    // cout << "ans: " << winContainer->configAdjustedWeight(S1) << endl << endl;

    // Config S2 = Config::fromString("[(0, 4, 4, 3, 3), (0, 14, 15, 3, 3), (), (0, 58, 38, 8, 8)]");
    // cout << "ans: " << winContainer->configAdjustedWeight(S2) << endl << endl;

    // Config S3 = Config::fromString("[(0, 4, 4, 3, 3), (0, 17, 20, 6, 8), (0, 58, 38, 8, 8)]");
    // cout << "ans: " << winContainer->configAdjustedWeight(S3) << endl << endl;


    // exit(0);



    // Config S1 = Config::fromString("[(0, 4, 4, 3, 3), (0, 12, 13, 1, 1), (0, 17, 20, 4, 4), (0, 23, 29, 3, 3), (0, 58, 38, 8, 8)]");
    // cout << winContainer->configAdjustedWeight(S1) << endl << endl << endl;

    // Config S2 = Config::fromString("[(0, 4, 4, 3, 3), (0, 13, 14, 2, 2), (0, 17, 20, 4, 4), (0, 23, 29, 3, 3), (0, 58, 38, 8, 8)]");
    // cout << winContainer->configAdjustedWeight(S2) << endl << endl << endl;


    // [(0, 4, 4, 3, 3), (0, 17, 20, 4, 4), (0, 23, 29, 3, 3), (0, 58, 38, 8, 8)]
/*
    Config S1 = Config::fromString("[(0, 4, 4, 3, 3), (0, 14, 15, 3, 3), (0, 17, 20, 2, 2), (0, 58, 38, 8, 8)]");
    cout << winContainer->configAdjustedWeight(S1) << endl ;

    Config S2 = Config::fromString("[(0, 4, 4, 3, 3), (0, 14, 15, 3, 3), (0, 17, 20, 2, 2), (0, 23, 29, 3, 3), (0, 58, 38, 8, 8)]");
    cout << winContainer->configAdjustedWeight(S2) << endl;

    Config S3 = Config::fromString("[(0, 4, 4, 3, 3), (0, 17, 20, 6, 8), (0, 23, 29, 3, 3), (0, 58, 38, 8, 8)]");
    cout << winContainer->configAdjustedWeight(S3) << endl << endl;

    Config S4 = Config::fromString("[(0, 23, 29, 3, 3), (0, 58, 38, 8, 8)]");
    cout << winContainer->configAdjustedWeight(S4) << endl;

    Config S5 = Config::fromString("[(0, 23, 29, 3, 3)]");
    cout << winContainer->configAdjustedWeight(S5) << endl;

    Config S6 = Config::fromString("[(0, 58, 38, 8, 8)]");
    cout << winContainer->configAdjustedWeight(S6) << endl << endl;

*/
        
    // Config S2 = Config::fromString("[(0, 5, 5, 4, 4), (0, 16, 18, 3, 3), (1, 30, 6, 5, 5), (1, 36, 14, 3, 3)]");
    // cout << winContainer->configAdjustedWeight(S2) << endl << endl << endl;

    // Config S3 = Config::fromString("[(0, 5, 5, 4, 4), (0, 16, 18, 3, 3), (1, 28, 6, 5, 5), (1, 36, 14, 3, 3)]");
    // cout << winContainer->configAdjustedWeight(S3) << endl << endl << endl;

    // Config S4 = Config::fromString("[(0, 7, 7, 2, 2), (1, 16, 18, 4, 4), (1, 19, 23, 2, 2), (2, 7, 10, 2, 2), (3, 13, 12, 3, 3)]");
    // cout << winContainer->configAdjustedWeight(S4) << endl << endl << endl;
    
    // cout << "is valid = " << S4.isValid() << endl;

    // exit(0);

    /*
    auto begin = chrono::high_resolution_clock::now();  

    for(auto win : winContainer->allWins)
    {
        // cout << win << endl;
        vector<Window> vec {win};
        int replacing, totalSize;
        bool staying;

        winContainer->sampleAddOrReplace(vec, 0, replacing, totalSize, staying);
    }

    auto end = chrono::high_resolution_clock::now();    
    auto dur = end - begin;
    auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(dur).count();

    int NN = (winContainer->maxRNALen) * (winContainer->maxRNALen) * (winContainer->maxRNALen) * (winContainer->maxRNALen) ;
    // double avg = (double)ms / winContainer->allWins.size();
    // cout << "Avg Time: " << (double)ms << " / " << winContainer->allWins.size() << " = " << avg << endl;
    double avg = (double)ms / NN;
    // cout << "Avg Time: " << (double)ms << " / " << NN << " = " << avg << endl;
        
    */
}


Sampler::~Sampler()
{
    delete winContainer;
}

Config Sampler::getSample()
{
    auto begin = chrono::high_resolution_clock::now();

    static int sC = 0;
    sC++;

    // cout << "Restarting..." << endl;
    // Config clean;
    // lastState = clean;
    // cout << "lastState = " << lastState << endl;



    // New idea: 
    // first, random walk for 50
    // Second, burn MCMC for 50
    // Third, get 50 samples from MCMC and choose one uniformly


    // step 1, random walk for 50
if(VERBOSE)
    cout << "Random Restarting... level " << curr_level_even << ", " << curr_level_odd << endl;

    Config orig;
    for(auto win: lastState)
    {
        orig.add(win);
    }

    
    Config T;

    // Restart only this level
    for(auto win: lastState)
    {
        if(win.l1 != curr_level_even or win.l2 != curr_level_odd)
            T.add(win);
    }




if(VERBOSE)
    cout << "post-clearing: " << T << endl;



    // At this point, config T contains no windows in level "currLevel"
    // During sampling of currLevel, we will initially add currLevel windows
    // to T, and then add or replace or delete during "randStartSteps" 
    // many steps.
    // 
    // Cycle detection:
    // We don't want to add or replace windows that will create cycles
    // Approach: consider the (original) T. Any window from currLevel that, when
    // added to T, creates a cycle, should not be considered for addition or 
    // replacement. 
    // We will look for all such windows and "delete" them from the set of allowed
    // windows.

    // We will incorporate the approach used in the previous version of this 
    // program: first find all paths in T. i.e., for any two regions A,B in T
    // we should be able to in O(1) whether A~>B or B~>A


    // Create graph from T, and get all paths
    unordered_map<string, vector<Region>> adjList = 
            winContainer->createGraphFromConfig(T);
    unordered_map<pair<string, string>, bool, pairhash> isPath =  
            winContainer->findAllPaths(T, adjList);

    // Now go through all windows in currLevel, and see which ones create a cycle
    // WindowContainer uses a single bit array to keep track of this, so we need
    // to reset it each time.
    winContainer->clearWindowCycleBits();
    for(auto win : *(winContainer->getLevel(curr_level_pair_id)))
    {
        // Not sure why this function requires T: it doesn't use it
        bool canAddWin = winContainer->addWinAndTest(T, adjList, win, isPath);
        if(!canAddWin)
            winContainer->setWindowCycleBit(win);
    }
    



    double _;
    for(int i=0;i<randStartSteps;i++)
    {
        std::tie(T, _, _) = sampleOneNeighbor(T);
        if(VERBOSE)
        cout << endl;
    }
if(VERBOSE)
    cout << "After " << randStartSteps << " steps of random walk, T = " << T << endl;
    // cout << "lastState = " << lastState << endl;


    // cout << endl << endl;
if(VERBOSE)
    cout << endl;

    // return T;
    lastState = T;
if(VERBOSE)
    cout << "lastState = " << lastState << endl;

    // Step 2, burn MCMC for 50
if(VERBOSE)
    cout << "Now burning in MCMC" << endl;
    // for(int i=0; i<intervalSteps; i++)
    for(int i=0; i<burnSteps; i++)
    {
        getNext();
    }
if(VERBOSE)
    cout << "Done burn in in MCMC (for " << burnSteps  << ")" << endl;
    // Step 3, get 50 samles from MCMC

if(VERBOSE)
    cout << "Generating " << burnSteps << " from MCMC" << endl;
    vector<Config> mcmcSamples;
    for(int i=0; i<intervalSteps; i++)
    {
        getNext();
        mcmcSamples.push_back(lastState);
if(VERBOSE)
        cout << "    for random >>  " << lastState << endl;
    }
if(VERBOSE){
    cout << "Done  from MCMC" << endl;

    cout << "Orig: " << orig << endl;
    cout << "For uniform: " << endl;
    for(auto x : mcmcSamples)
        cout << x << endl;
    // cout << endl;
    }
    
    int random = rand() % intervalSteps;
    lastState =  mcmcSamples[random];
if(VERBOSE)
    cout << "Uniformly Picked " << lastState << " as sample" << endl;
if(VERBOSE)
    cout << endl << endl;


    auto end = chrono::high_resolution_clock::now();    
    auto dur = end - begin;
    auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(dur).count();
    if(DEBUG_MODE)cout << "**Time for one sample (# " << sC << ") :" << (double)ms << "ms" << endl;



    return lastState;
}



void Sampler::getNext()
{
    static int sc = 0;
    sc++;
    Config S = lastState;

    double Qf, Qb;  //, Qb = lastStateProb;
    Config T;

    auto begin = chrono::high_resolution_clock::now();  
    

    std::tie(T, Qf, Qb) = sampleOneNeighbor(S);


    auto end = chrono::high_resolution_clock::now();    
    auto dur = end - begin;
    auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(dur).count();
    if(VERBOSE)cout << "Time to make back,forward :" << (double)ms << "ms" << endl;

    S.sort();
    T.sort();
if(VERBOSE)
    cout << "S = " << (S) << endl;
if(VERBOSE)
    cout << "T = " << (T) << endl;
    
    // double S_weight = winContainer->configAdjustedWeight(S);    //S.totalWeight();
    // double T_weight = winContainer->configAdjustedWeight(T);    //T.totalWeight();
    double S_weight = scorer->configAdjustedWeight(S);    //S.totalWeight();
    double T_weight = scorer->configAdjustedWeight(T);    //T.totalWeight();

if(VERBOSE)
    cout << "Weight(S) = " << S_weight << endl;
if(VERBOSE)
    cout << "Weight(T) = " << T_weight << endl;

    //Qf = Qb = 1;
    double RT = 0.616321;
    // double RT = 1;
    double Z_ratio = exp((-1*T_weight - (-1)*S_weight)/RT );
    ////cout << "Z_T / Z_S = " << Z_T / Z_S << endl;
if(VERBOSE)
    cout << "Z_T / Z_S = " << Z_ratio << endl;
if(VERBOSE)
    cout << "Qf = " << Qf << endl;
if(VERBOSE)
    cout << "Qb = " << Qb << endl;

    //cout << "Qb / Qf = " << (Qb / Qf) << endl;

    double ratio = (Qb / Qf) * Z_ratio;
if(VERBOSE)
    cout << "ratio = " << ratio << endl;

    double r = ((double) rand() / (RAND_MAX));
if(VERBOSE)
    cout << "r = " << r << endl;

    //Accept with prob min(1, ratio)
    if(r < MIN(1, ratio))
    {
    if(VERBOSE)
        cout << "Accepted" << endl;
        accept_count++;
    
        lastState = T;
        lastStateProb = Qf;

        /*
        if(T.size() < S.size())
            cout << "**case 1" << endl;
        else if(T.size() > S.size())
            cout << "**case 2" << endl;
        else
            cout << "**case 3" << endl;
        */
    }
    else
    {
        /*
        if(T.size() < S.size())
            cout << "**case 4" << endl;
        else if(T.size() > S.size())
            cout << "**case 5" << endl;
        else
            cout << "**case 6" << endl;
        */

        reject_count++;
    if(VERBOSE)
        cout << "Rejected" << endl;
    }
    if(VERBOSE)
    cout << "Sample count: " << sc << endl;
if(VERBOSE)
    cout << endl;

}


int sampleFromDistribution(vector<double> dist)
{
    double r = ((double) rand() / (RAND_MAX));
    double cumm = 0;

    // for(int i=0; i<dist.size(); i++)
    //  cout << dist[i] << " ";
    // cout << endl;

    for(int i=0; i<dist.size(); i++)
    {
        cumm += dist[i];
        if(r <= cumm)
        {
            return i;
        }
    }

    throw std::runtime_error("couldn't sample from probability distribution"); 
    // return -1;   //to keep compiler happy
}




















































        




        
int Sampler::getNeighborhoodCountAndMatch(Config S, Config compare)
{
    return winContainer->getNbrSizeOnly(S.getVector(), curr_level_even, curr_level_odd, curr_level_pair_id);
}


std::tuple<Config, double, double> Sampler::sampleOneNeighbor(Config S)
{
    if(VERBOSE)
    cout << "Start with " << S << endl;
    

    // getAddable(S);


    vector<Window> vec = S.getVector();
    int replacing = -1, totalNbrs = 0;
    bool staying_here = false;
    int choice = winContainer->sampleAddOrReplace(vec, curr_level_even, curr_level_odd, curr_level_pair_id, replacing, totalNbrs, staying_here);
    if(DEBUG_MODE)cout << ">> choice = " << choice << endl;
    if(DEBUG_MODE)cout << "<< totalNbrs = " << totalNbrs << endl;

    
    Config newConfig;

    if(staying_here)
    {
        for(int i=0;i<vec.size();i++)
                newConfig.add(vec[i]);
    }
    else if(replacing != -1)
    {
        if(DEBUG_MODE)cout << "** BIT: replacing window # " << replacing << ": " << endl;
        if(DEBUG_MODE)cout << vec[replacing] << endl;
        if(DEBUG_MODE)cout << "with" << winContainer->allWins[choice] << endl;
        if(DEBUG_MODE)cout << "where choice = " << choice << endl;
        // Replacing... replace window "replacing" with "choice"
        for(int i=0;i<vec.size();i++)
            if(i != replacing)
                newConfig.add(vec[i]);

        newConfig.add(winContainer->allWins[choice]);
    }
    else if(choice < 0)
    {
        if(DEBUG_MODE)cout << "** BIT: deleting window # " << choice*-1 << ": " << endl;
        if(DEBUG_MODE)cout << vec[(-1*choice)-1] << endl;
        for(int i=0;i<vec.size();i++)
            if(i != (-1*choice)-1 )
                newConfig.add(vec[i]);
    }
    else
    {
        if(DEBUG_MODE)cout << "** BIT: ADding window # " << choice << ": " << endl;
        if(DEBUG_MODE)cout << "this " << winContainer->allWins[choice] << endl;
        
        for(int i=0;i<vec.size();i++)
            newConfig.add(vec[i]);
        newConfig.add(winContainer->allWins[choice]);
    }
    // cout << "XXXXX: " <<  << endl;
    newConfig.sort();

// cout << "YYYYY" << endl;




if(DEBUG_MODE)cout << endl << "** Final Forward = " << newConfig << endl << endl;
if(DEBUG_MODE)cout << endl << "** Final Forward totalNbrs = " << totalNbrs << endl << endl;
    double forwardProb = 1.0/totalNbrs;
    // Now compute backward prob
    // For this, get nbrs of newConfig
    int reverseCount = getNeighborhoodCountAndMatch(newConfig, S);
    double backProb = 1.0/reverseCount;
    

    // cout << "will stop here" << endl;
    //backProb = forwardProb = 1;
    std::tuple<Config, double, double> val = std::make_tuple(newConfig, forwardProb, backProb);
    return val;
    

}

void Sampler::setLevel(int l)
{
    curr_level_pair_id = l;

    std::pair<int,int> p = level_id_to_pair[l];
    curr_level_even = std::get<0>(p);
    curr_level_odd = std::get<1>(p);
}

double Sampler::getWeight(Config c)
{
    // return winContainer->configAdjustedWeight(c);
    return scorer->configAdjustedWeight(c);
}

int Sampler::getNumOfRNAs()
{
    return numOfRNAs;
}


WindowContainer* Sampler::getWinContainer()
{
    return winContainer;
}