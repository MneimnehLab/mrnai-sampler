#ifndef _WINDOW_CONTAINER_H_
#define _WINDOW_CONTAINER_H_

#include <iostream>
#include <vector>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include "Window.h"
#include "Region.h"
#include "Config.h"
#include "BitSetTools.h"
// #include "scoring/AdjacentScoring.h"



using std::vector;
using std::string;
using std::pair;


struct pairhash {
public:
  template <typename T, typename U>
  std::size_t operator()(const std::pair<T, U> &x) const
  {
    return std::hash<T>()(x.first) ^ std::hash<U>()(x.second);
  }
};

struct triplehash {
public:
  template <typename T, typename U, typename V>
  std::size_t operator()(const std::tuple<T, U, V> &x) const
  {
    return std::hash<T>()(std::get<0>(x)) ^ std::hash<U>()(std::get<1>(x)) ^ std::hash<U>()(std::get<2>(x));
  }
};


class abc {
  std::unordered_map<std::pair<int,int>, int, pairhash> rules;
};


typedef unordered_map<string, vector<Window> *> WinMap;
typedef unordered_map<string, mrnai_tools::BitSetTools::BIT_STORE*> WinMapLong;
typedef unordered_map<pair<int,int>, mrnai_tools::BitSetTools::BIT_STORE *, pairhash> pair_bitmap;


class WindowContainer
{
	friend class AdjacentScoring;
private:
	vector<vector<Window> *> windowsByLevel;
	
	// WinMap overlapsHM;
	// WinMap nonOverlapsHM;
	

	
	bool hasMadeOverlaps;

	double bestPairingRegions(const vector<Region> & regions) ;
	double regionPairScore(const Region & x, const Region & y) ;

	void sweepRightToLeft(int level, pair_bitmap & winBitsOnRightOfEdge);
	void sweepLeftToRight(int level, pair_bitmap & winBitsOnLeftOfEdge);

	void map_corners_to_windows(vector<Window> * vec);

	int bitVecLongs;

	WinMapLong nonOverlapsHM_Long;
	WinMapLong overlapsHM_Long;


	unordered_map<int, mrnai_tools::BitSetTools::BIT_STORE*> levelMapper;
	
	
	mrnai_tools::BitSetTools::BIT_STORE* winThatCreateCycles;




	
public:
int maxRNALen;
	int numEven;
	int numOdd;

	vector<int> even_levels;
	vector<int> odd_levels;

	vector<int> level_pair_to_id; // use as a 2D array
	vector<std::pair<int, int>> level_id_to_pair;

	// We now have to distinguish between number of RNAs and number of levels
	// number of RNAs   = even + odd
	// number of Levels = even * odd
	int numOfLevels;
	int numOfRNAs;

	int numOfWins;


	// unordered_map<std::tuple<int,int,int>, vector<Window> *, triplehash> windowsByEndingEdge;


	vector<Window> allWins;
	int sampleAddOrReplace(const vector<Window> & vec, int level1, int level2, int level_pair_id, int & replacing, int & totalSize, bool & staying);
	int getNbrSizeOnly(const vector<Window> & vec, int level1, int level2, int level_pair_id);

	
	void map_corners_to_regions(vector<Window> * vec);
	// intervals_by_right_corner[level] -> hashtable of interval_str -> 
	vector<vector<vector<int>>> regions_by_right_corner;
	vector<vector<vector<int>>> regions_by_left_corner;

	unordered_map<std::pair<int,int>, vector<int>, pairhash> map_windows_by_right_corner;
	unordered_map<std::pair<int,int>, vector<int>, pairhash> map_windows_by_left_corner;



	
	unordered_map<string, double> interEs;
	// Need to explictly look up upE of a region, so this should be in a hashtable... don't care about per level, I guess
	unordered_map<string, double> regionUpEs;

	
	
	WindowContainer(vector<int>, vector<int>) ;
	~WindowContainer();

	void addWindow(Window win);
	void addWindow(Window win, double topUpE, double botUpE);
	Window addWindow(string);

	void addEvenRegion(Window win, double topUpE);
	void addOddRegion(Window win, double botUpE);

	const vector<Window> * getLevel(int i);

	void makeOverlaps();
	void makeOverlaps_Slow();
	void makeOverlaps_Slow_Half();
	
	void makeOverlaps_with_intervals();

	
	
	int getNumOfLevels();
	
	double configAdjustedWeight(Config & vec);

	unordered_map<string, vector<Region>> createGraphFromConfig(Config S);
	unordered_set<string> simpleBFS(string start, unordered_map<string, vector<Region>>&  adjList);
	unordered_map<pair<string, string>, bool, pairhash> findAllPaths(Config S, unordered_map<string, vector<Region>>&  adjList);
	bool addWinAndTest(Config S, 
			unordered_map<string, vector<Region>>&  adjList, 
			Window newWin, 
			unordered_map<pair<string, string>, bool, pairhash> & isPath);

	void setWindowCycleBit(Window win);
	void clearWindowCycleBits();


	double new_regionScore(const vector<Region> & regions, int start, int end) ;
	double new_unionScore(const vector<Region> & allRegions);

	double best_pairing(const vector<Region> & S, int i, int j);

	static bool regionComparator (Region i, Region j);

};

#endif