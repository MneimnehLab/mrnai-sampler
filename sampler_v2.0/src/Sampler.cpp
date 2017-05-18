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

#define LOG(args) {if(VERBOSE) {std::cout << (args);}

Sampler::Sampler(const Sampler::Params & params)
{
    intervalSteps   = params.intervalSteps;
    burnSteps       = params.burnSteps;
    maxWinSize      = params.maxWinSize;
    randStartSteps  = params.randStartSteps;

    int numEven     = params.numEven;
    int numOdd      = params.numOdd;
    bool equalWins  = params.equalWins;
    // double threshold = params.threshold;

    curr_level_pair_id = 0;
    curr_level_even = params.even_levels.at(0);
    curr_level_odd  = params.odd_levels.at(0);

    numOfRNAs = numEven + numOdd;


    // level_pair_to_id.resize((numOfRNAs+1)*(numOfRNAs+1));
    // for(auto& e : level_pair_to_id)
    //     e.resize(numOfRNAs+1);


    int id = 0;
    for(int even : params.even_levels)
        for(int odd : params.odd_levels)
        {
            // level_pair_to_id[even][odd] = id;
            level_id_to_pair.push_back(std::make_pair(even, odd));
            id++;
        }


    cerr << "* Num Even : " << numEven << endl;
    cerr << "* Num Odd : " << numOdd << endl;
    
    winContainer = new WindowContainer(params.even_levels, params.odd_levels);
    
    vector<std::tuple<Window, double, double> > tempContainerToSort;

    int maxLevel = 0;
    for (string line; getline(cin, line);) 
    {
        // Parse line here:
        stringstream ss;
        ss << line;
        int level1, level2, i_, i, j_, j;
        
        ss >> level1 >> level2 >> i_ >> i >> j_ >> j;

        // Read energies
        double intlnZ, totlUpE, delta, upe1, upe2;
        ss >> delta >> totlUpE >> intlnZ >> upe1 >> upe2;;

        Window win(level1, level2, i, j, i-i_, j-j_, delta, intlnZ);

        // Register all regions' weights, regardless of whether these 
        // windows get filtered out
        winContainer->addEvenRegion(win, upe1);
        winContainer->addOddRegion (win, upe2);

        // Filter windows, based on program runtime params
        // if(delta > 0)
        if(intlnZ > 0)
            continue;

        if(equalWins && win.w1 != win.w2)
            continue;

        if(win.w1 > maxWinSize or win.w2 > maxWinSize)
            continue;

        tempContainerToSort.push_back(std::make_tuple(win, upe1, upe2));
    }

    
    // Sort container based on energy: indx 0 => min energy
    std::sort(tempContainerToSort.begin(), tempContainerToSort.end(), 
        [](tuple<Window, double, double> a, tuple<Window, double, double> b) 
    {
        return std::get<0>(a).intrEnergy < std::get<0>(b).intrEnergy;   
    });


    // Now add "good" windows to winContainer
    int maxGoodWindows = params.maxGoodWins;

    for(int i=0; i<maxGoodWindows && i<tempContainerToSort.size(); i++)
    {
        Window win(0,0,0,0,0,0);
        double upe1, upe2;
        std::tie(win, upe1, upe2) = tempContainerToSort[i];

        win.id = i;
        winContainer->addWindow(win);

        maxLevel = MAX( MAX(win.l1, win.l2), maxLevel);
    }

    // After adding windows, create overlap matrices
    winContainer->makeOverlaps_with_intervals();    

    // Choose one scoring scheme (TODO: make program param)
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
}


Sampler::~Sampler()
{
    delete winContainer;
}

// Note: this gets one sample from the Gibbs Sampler
// so here we deal with many samples of Metropolis Hastings
Config Sampler::getSample()
{
    auto begin = chrono::high_resolution_clock::now();

    static int sC = 0;
    sC++;

    // The algorithm has three phases: 
    // first, random walk for 50
    // Second, burn MCMC for 50
    // Third, get 50 samples from MCMC and choose one uniformly


    // step 1, random walk for 50

    Config orig;
    for(auto win: lastState)
    {
        orig.add(win);
    }

    
    // Restart only this level, by cleaning out all wins in lastState that are from this level
    Config T; // T will contain all those windows that do not have even,odd == l1,l2
    for(auto win: lastState)
    {
        if(win.l1 != curr_level_even or win.l2 != curr_level_odd)
            T.add(win);
    }


    // At this point, config T contains no windows in level "currLevel*"
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
    

    // Now perform random walk for "randStartSteps" steps
    double _; // two returned vals of sampleOneNeighbor are probabilities,  
              // which we don't care about for now
    for(int i=0;i<randStartSteps;i++)
        std::tie(T, _, _) = sampleOneNeighbor(T);
    

    // lastState is the last state to by Metropolis Hastings
    // Each call to getNext() updates lastState, since it is hte last state visited by MH
    // initialize it to T.
    lastState = T;

    // Step 2, burn MCMC for 50
    
    for(int i=0; i<burnSteps; i++)
        // note: this will update lastState
        // states are accepted/rejected  per MH criterion
        getNext();


    // Step 3, get 50 samples from MCMC

    vector<Config> mcmcSamples;
    for(int i=0; i<intervalSteps; i++)
    {
        // again, getNext() updates lastState, so we can can add it to our vector
        getNext();
        mcmcSamples.push_back(lastState);
    }
    
    int random = rand() % intervalSteps;
    lastState =  mcmcSamples[random];

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

    double Qf, Qb;
    Config T;

    auto begin = chrono::high_resolution_clock::now();  
    
    std::tie(T, Qf, Qb) = sampleOneNeighbor(S);

    auto end = chrono::high_resolution_clock::now();    
    auto dur = end - begin;
    auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(dur).count();
    if(VERBOSE)cout << "Time to make back,forward :" << (double)ms << "ms" << endl;

    // Commenting out sort(), since I don't think anything depends on sorted configs...
    // S.sort();
    // T.sort();

    if(VERBOSE) cout << "S = " << (S) << endl;
    if(VERBOSE) cout << "T = " << (T) << endl;
    
    // get weights for both current and candidate states, 
    // using the specified dependency model
    double S_weight = scorer->configAdjustedWeight(S);
    double T_weight = scorer->configAdjustedWeight(T);

    if(VERBOSE) cout << "Weight(S) = " << S_weight << endl;
    if(VERBOSE) cout << "Weight(T) = " << T_weight << endl;

    double RT = 0.616321;
    double Z_ratio = exp((-1*T_weight - (-1)*S_weight)/RT );
    
    if(VERBOSE) cout << "Z_T / Z_S = " << Z_ratio << endl;
    if(VERBOSE) cout << "Qf = " << Qf << endl;
    if(VERBOSE) cout << "Qb = " << Qb << endl;

    double ratio = (Qb / Qf) * Z_ratio;

    if(VERBOSE) cout << "ratio = " << ratio << endl;

    double r = ((double) rand() / (RAND_MAX));
    if(VERBOSE) cout << "r = " << r << endl;

    //Accept with prob min(1, ratio)
    if(r < ratio)
    {
        if(VERBOSE) cout << "Accepted" << endl;
        accept_count++; // for testing...
    
        lastState = T;
    }
    else
    {
        reject_count++; // for testing...
        if(VERBOSE) cout << "Rejected" << endl;

        // if rejected, then lastState remains lastState
    }
    
    if(VERBOSE) cout << "Sample count: " << sc << endl;
    if(VERBOSE) cout << endl;
}


int Sampler::getNeighborhoodCountAndMatch(Config S, Config compare)
{
    return winContainer->getNbrSizeOnly(S.getVector(), curr_level_even, curr_level_odd, curr_level_pair_id);
}


// This performs one step of the random walk, by generating one neighbor randomly
// If current state is S, it generates a candidate C, and also finds the probabilites
// of going in both directions
// The Probabaility P(S|C) = 1 / |neighbors of C|, so it computes |neighbors of C|
std::tuple<Config, double, double> Sampler::sampleOneNeighbor(Config S)
{
    vector<Window> vec = S.getVector();
    int replacing = -1, totalNbrs = 0;
    bool staying_here = false;
    
    int choice = winContainer->sampleAddOrReplace(vec, curr_level_even, curr_level_odd, curr_level_pair_id, 
                                                    replacing, totalNbrs, staying_here);
    
    Config newConfig;

    // Using rules specified by "sampleAddOrReplace", figure out which action
    // we are taking here.
    if(staying_here)
    {
        for(int i=0;i<vec.size();i++)
                newConfig.add(vec[i]);
    }
    else if(replacing != -1)
    {
        // Replacing... replace window "replacing" with "choice"
        for(int i=0;i<vec.size();i++)
            if(i != replacing)
                newConfig.add(vec[i]);

        newConfig.add(winContainer->allWins[choice]);
    }
    else if(choice < 0)
    {
        for(int i=0;i<vec.size();i++)
            if(i != (-1*choice)-1 )
                newConfig.add(vec[i]);
    }
    else
    {
        for(int i=0;i<vec.size();i++)
            newConfig.add(vec[i]);
        newConfig.add(winContainer->allWins[choice]);
    }
    newConfig.sort();

    double forwardProb = 1.0/totalNbrs;
    
    // Now compute backward prob
    // For this, get nbrs of newConfig
    int reverseCount = getNeighborhoodCountAndMatch(newConfig, S);
    double backProb = 1.0/reverseCount;    

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
