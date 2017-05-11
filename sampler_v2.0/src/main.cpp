#include <chrono>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <unordered_map>
#include <cstdlib>
#include <ctime>
#include <algorithm>
#include <chrono>
#include "Window.h"
#include "WindowContainer.h"
#include "Sampler.h"
#include "Config.h"
using namespace std;

#define DOTIME 0

char* getCmdOption(char ** begin, char ** end, const std::string & option);
bool cmdOptionExists(char** begin, char** end, const std::string& option);
bool compar (std::tuple<string, int, double> i, std::tuple<string, int, double> j) ;


int main(int argc, char *argv[])
{
    srand(time(0));


    if(cmdOptionExists(argv, argv+argc, "-h"))
    {
        cout << "Usage: ./main  [-f output_filename]";
        cout << "               [-s number of samples]";

        return 0;
    }

    string structFilename = "mystructs";
    
    char * filename = getCmdOption(argv, argv + argc, "-f");

    if (filename)
    {
        structFilename = string(filename);
    }

    char * analysis_type = getCmdOption(argv, argv + argc, "-t");
    Sampler::Params::AType a_type = Sampler::Params::INTERVAL;
    if (analysis_type)
    {
        if(string(analysis_type) == "fast")
            a_type = Sampler::Params::FAST;
        else if(string(analysis_type) == "slow")
            a_type = Sampler::Params::SLOW;
        else if(string(analysis_type) == "slow_half")
            a_type = Sampler::Params::SLOW_HALF;
    }

    char * numSamplesChr = getCmdOption(argv, argv + argc, "-s");
    int numSamples = 100;
    if (numSamplesChr)
    {
        //cout << "numSamplesChr = " << numSamplesChr << endl;
        numSamples = stoi(string(numSamplesChr));
    }

    char * intervalStepsChr = getCmdOption(argv, argv + argc, "-i");
    int intervalSteps = 50;
    if(intervalStepsChr)
        intervalSteps = stoi(string(intervalStepsChr));

    char * burnStepsChr = getCmdOption(argv, argv + argc, "-b");
    int burnSteps = 50;
    if(burnStepsChr)
        burnSteps = stoi(string(burnStepsChr));


    char * randStartStepsChr = getCmdOption(argv, argv + argc, "-r");
    int randStartSteps = 20;
    if(randStartStepsChr)
        randStartSteps = stoi(string(randStartStepsChr));


    char * maxWinSizeChr = getCmdOption(argv, argv + argc, "-w");
    int maxWinSize = 25;
    if(maxWinSizeChr)
        maxWinSize = stoi(string(maxWinSizeChr));

    char * maxGoodWinChr = getCmdOption(argv, argv + argc, "-m");
    int maxGoodWin = 10000000;
    if(maxGoodWinChr)
        maxGoodWin = stoi(string(maxGoodWinChr));

    // cout << "Pair Sampling with Split and Merge" << endl;
    // cout << "----------------------------------" << endl;
    // cout << "Generating " << numSamples << " samples" << endl;
    // cout << "Writing to " << structFilename << " samples" << endl;

    // Sampler sampler(intervalSteps, burnSteps);


    int numEven = -1, numOdd = -1;
    vector<vector<int>> even_odd_levels(2);

    for(int i=0; i<2; i++)
    {
        string line;
        getline(cin, line);
        stringstream ss(line);
        
        int rnaNum = -1;
        while(ss >> rnaNum)
        {
            even_odd_levels[i].push_back(rnaNum);
        }

        if(rnaNum == -1)
            throw std::runtime_error("\nInvalid RNA Number in input\n\n"); 
    }

    numEven = even_odd_levels[0].size();
    numOdd = even_odd_levels[1].size();

    cout << "Even rnas: "<< endl;
    for(auto e : even_odd_levels[0])
        cout << e << endl;
    cout << endl;

    cout << "Odd rnas: "<< endl;
    for(auto e : even_odd_levels[1])
        cout << e << endl;
    cout << endl;

    
    // Read bipartite graph of the optimal solution
    vector<vector<int>> adj_matrix(numEven+numOdd, vector<int>(numEven+numOdd, 0));
    {
        string line;
        getline(cin, line);
        stringstream ss(line);
        int even, odd;
        while(ss >> even >> odd)
            adj_matrix[even][odd] = 1;
    }
    
    Sampler::Params params;

    params.equalWins = false;
    if(cmdOptionExists(argv, argv+argc, "-e"))
        params.equalWins = true;
    params.maxGoodWins = maxGoodWin;
    params.maxWinSize = maxWinSize;
    params.threshold = 99999;
    // params.threshold = 0.2;
    params.intervalSteps = intervalSteps;
    params.burnSteps = burnSteps;
    params.randStartSteps = randStartSteps;
    params.numEven = numEven;
    params.numOdd = numOdd;
    params.even_levels = even_odd_levels[0];
    params.odd_levels = even_odd_levels[1];
    params.atype = a_type;

    Sampler sampler(params);


    // exit(0);


    // return 0;

    // int numRNA = sampler.getNumOfRNAs();

    unordered_map<string, int> structCount;

    int T = numSamples;
    Config S;
    int pstep = 10;

    auto begin_trials = chrono::high_resolution_clock::now();  


    
    for(int i=0;i<T;i++)
    {
        if((i+1) % (int)( (double)(T) / ( (double)(100) / pstep))  == 0)
        {
            int percent = (i+1) / (double)(T) * 100;
            cerr << percent << "% iterations done" << endl;
        }

        #if DOTIME == 1
        auto begin = chrono::high_resolution_clock::now();  
        #endif

        int level_pair_id = 0;
        for(int evenId : even_odd_levels[0])
        {
            for(int oddId : even_odd_levels[1])
            {
                // Disable edge
                // if(evenId == 1 and oddId == 0)
                //     continue;

                if(adj_matrix[evenId][oddId] == 1)
                {
            
                    // auto begin = chrono::high_resolution_clock::now();  

                    sampler.setLevel(level_pair_id);
                    S = sampler.getSample();

                    // auto end = chrono::high_resolution_clock::now();    
                    // auto dur = end - begin;
                    // auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(dur).count();
                    
                    // cout << "Time  for level " << level << " sample:" << (double)ms << "ms" << endl;
                }
                
                level_pair_id++;
            
            }
        }

        Config S_ = sorted(S);
        structCount[S_.toString()] += 1;
        
        #if DOTIME == 1
        auto end = chrono::high_resolution_clock::now();    
        auto dur = end - begin;
        auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(dur).count();
        
        cout << "Time  after 3 levels:" << (double)ms << "ms" << endl;
        #endif

    }


    auto end_trials = chrono::high_resolution_clock::now();    
    auto dur_trials = end_trials - begin_trials;
    auto ms_trials = std::chrono::duration_cast<std::chrono::milliseconds>(dur_trials).count();
    
    cout << "Sampling duration:" << (double)ms_trials << "ms" << endl;

    vector<std::tuple<string, int, double> > vals;
    // cout << endl << "Struct count:: " << endl << endl;

    //for (unordered_map<string, int, double>::iterator it = structCount.begin(); it != structCount.end(); ++it) 
    for (unordered_map<string, int>::iterator it = structCount.begin(); it != structCount.end(); ++it) 
    {
        //std::cout << it->first << ": " << it->second << endl;;

        Config c = Config::fromString(it->first);
         

        vals.push_back(std::make_tuple(it->first, it->second, sampler.getWeight(c)));
    }

    sort(vals.begin(), vals.end(), compar);

    for(int i=0; i<10; i++)
    {
        std::tuple<string, int, double> val = vals[i];
        
        int count;
        string str;
        double wt;
        std::tie(str, count, wt) = val;
        Config c = Config::fromString(str);
        std::cout << c << ": " << count << "   (weight = " << wt << ")" << endl;;
    }

    std::cout << std::endl; 

    cout << "accept_count = " << sampler.accept_count << endl;
    cout << "reject_count = " << sampler.reject_count << endl;

    
    ofstream structsOut(structFilename.c_str());

    int n = 1;
    for(std::tuple<string, int, double> val : vals)
    {
        int count;
        string str;
        double wt;
        std::tie(str, count, wt) = val;
        Config c = Config::fromString(str);

        structsOut << n << "\t"
                   << c << "\t"
                   << wt << "\t"
                   << count << endl;
        n++;
    }

    structsOut.close();

    return 0;
    
}






char* getCmdOption(char ** begin, char ** end, const std::string & option)
{
    char ** itr = std::find(begin, end, option);
    if (itr != end && ++itr != end)
    {
        return *itr;
    }
    return 0;
}

bool cmdOptionExists(char** begin, char** end, const std::string& option)
{
    return std::find(begin, end, option) != end;
}

bool compar (std::tuple<string, int, double> i, std::tuple<string, int, double> j) 
{ 
    return (std::get<2>(i) < std::get<2>(j)); 
}
