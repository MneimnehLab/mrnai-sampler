#include <iostream>
#include <bitset>
#include <algorithm>
#include <chrono>
#include <vector>
#include <string>
#include <set>
#include <unordered_map>
#include <fstream>
#include <queue>
#include <unordered_set>
#include "WindowContainer.h"
#include "Window.h"
#include "Region.h"
#include "util.h"
#include "Config.h"
#include "BitSetTools.h"

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


bool WindowContainer::regionComparator (Region i, Region j) { return (i<j); }


WindowContainer::WindowContainer(vector<int> even_levels, vector<int> odd_levels) 
    : hasMadeOverlaps(false), 
      even_levels(even_levels),
      odd_levels(odd_levels)
{
    numEven = even_levels.size();
    numOdd = odd_levels.size(); 
    numOfLevels = numEven*numOdd; 
    numOfRNAs = numEven+numOdd;
  
    Window::setNumEvenOdd(numEven, numOdd);

    level_pair_to_id.resize(numOfRNAs*numOfRNAs);
    
    // for(int i=0; i<numOfLevels; i++)
    int id = 0;
    for(int even : even_levels)
        for(int odd : odd_levels)
        {
            // level_pair_to_id[even][odd] = id;
            level_pair_to_id[INDEX(even, odd, numOfRNAs)] = id;
            level_id_to_pair.push_back(std::make_pair(even, odd));

            windowsByLevel.push_back(new vector<Window>());
            id++;
        }

    maxRNALen = 0;
}



WindowContainer::~WindowContainer()
{
    for(int i=0;i<numOfLevels;i++)
    {
        if(windowsByLevel[i] != nullptr)
        {
            delete windowsByLevel[i];
            // cout << "deleted level" << endl;
        }
    }
    // cout << "deleting" << endl;

    delete[] winThatCreateCycles;
}



inline void setBitTo1(int winNum, mrnai_tools::BitSetTools::BIT_STORE array[])
{
    int pos = winNum / (sizeof(mrnai_tools::BitSetTools::BIT_STORE)*8);
    int idx = winNum % (sizeof(mrnai_tools::BitSetTools::BIT_STORE)*8);

    // have to set a 1 at bit 'idx' (from right) in integer 'pos'
    mrnai_tools::BitSetTools::BIT_STORE bits = array[pos];
    // BIT_STORE one = 1; one = one << idx;
    mrnai_tools::BitSetTools::BIT_STORE mask = 0; 
    mask = ~((mask-1)>> 1);
    mask = mask >> idx;
    bits = bits | mask;
    array[pos] = bits;
}

inline void setBitTo0(int winNum, mrnai_tools::BitSetTools::BIT_STORE array[])
{
    int pos = winNum / (sizeof(mrnai_tools::BitSetTools::BIT_STORE)*8);
    int idx = winNum % (sizeof(mrnai_tools::BitSetTools::BIT_STORE)*8);

    mrnai_tools::BitSetTools::BIT_STORE bits = array[pos];
    mrnai_tools::BitSetTools::BIT_STORE mask = -1; 
    mask = ~(mask>> 1);
    mask = mask >> idx;
    mask = ~mask;
    bits = bits & mask;
    array[pos] = bits;
}






void WindowContainer::addEvenRegion(Window win, double topUpE)
{
    int even_level = win.l1;
    
    Region evenR(even_level,  win.i, win.w1, topUpE);
    regionUpEs[evenR.toString()] = topUpE;
}

void WindowContainer::addOddRegion(Window win, double botUpE)
{
    int odd_level = win.l2;
    
    Region oddR(odd_level,   win.j, win.w2, botUpE);
    regionUpEs[oddR.toString()] = botUpE;
}


void WindowContainer::addWindow(Window win)
{
    // cout << "Adding win: " << win << endl;
    int level_id = level_pair_to_id[INDEX(win.l1, win.l2, numOfRNAs)];  //[win.l1][win.l2];
    windowsByLevel[level_id]->push_back(win);
    interEs[win.toString()] = win.intrEnergy;
    allWins.push_back(win);

    maxRNALen = MAX(maxRNALen, MAX(win.i, win.j));
}


// Allows us to look up all windows that have an interval
// beginning or ending at a given point.
void WindowContainer::map_corners_to_windows(vector<Window> * vec)
{
    int counter = 0;
    for(const Window& win : *vec)
    {
        // for each window "win", add it to the four hashtables that correspond
        // to its four corners
        int even = win.l1;
        int odd  = win.l2;

        int even_right_corner = win.i;
        int even_left_corner  = win.i - win.w1;
        int odd_right_corner  = win.j;
        int odd_left_corner   = win.j - win.w2;

        pair<int, int> even_right = std::make_pair(even, even_right_corner);
        pair<int, int> even_left  = std::make_pair(even, even_left_corner);
        pair<int, int> odd_right  = std::make_pair(odd,  odd_right_corner);
        pair<int, int> odd_left   = std::make_pair(odd,  odd_left_corner);

        map_windows_by_right_corner[even_right].push_back(counter);
        map_windows_by_left_corner [even_left ].push_back(counter);
        map_windows_by_right_corner[odd_right ].push_back(counter);
        map_windows_by_left_corner [odd_left  ].push_back(counter);

        counter++;
    }
}

void WindowContainer::sweepRightToLeft(int level, pair_bitmap & winBitsOnRightOfEdge)
{
    // use iterators everywhere to avoid looking up the same (key,value) pair more than once
    for(int i=maxRNALen-1; i>0; i--) // n -> 0 direction crucial for correctness
    {
        // Idea: winBitsOnRightOfEdge[(level,i)] = winBitsOnRightOfEdge[(level, i+1)]  -> A
        //                                   union wins_starting_from[(level,i+1)]  -> B
        
        pair<int, int> lev_i   = std::make_pair(level,i);
        pair<int, int> lev_i_1 = std::make_pair(level,i+1);

        mrnai_tools::BitSetTools bst(numOfWins);
        mrnai_tools::BitSetTools::Container temp_bitset = bst.makeBitSet();

        // Init with A
        auto iterA = winBitsOnRightOfEdge.find(lev_i_1);
        if(iterA != winBitsOnRightOfEdge.end())
            bst.s_union(temp_bitset, iterA->second);

        // union with B
        auto gotL1 = map_windows_by_left_corner.find (lev_i_1);
        if ( gotL1 != map_windows_by_left_corner.end() )
        {
            vector<int> intervals = gotL1->second;
            for(int w_ : intervals)
                setBitTo1(w_, temp_bitset);
        }

        winBitsOnRightOfEdge[lev_i] = temp_bitset;
    }
}



void WindowContainer::sweepLeftToRight(int level, pair_bitmap & winBitsOnLeftOfEdge)
{
    for(int i=1; i<=maxRNALen; i++)
    {
        // Idea: winBitsOnLeftOfEdge[(level,i)] = winBitsOnLeftOfEdge[(level, i-1)]  -> A
        //                                       union intervals_endint_at[(level,i-1)]  -> B
        
        pair<int, int> lev_i   = std::make_pair(level,i);
        pair<int, int> lev_i_1   = std::make_pair(level,i-1);

        mrnai_tools::BitSetTools bst(numOfWins);
        mrnai_tools::BitSetTools::Container temp_bitset = bst.makeBitSet();

        // Init with A
        auto iterA = winBitsOnLeftOfEdge.find(lev_i_1);
        if(iterA != winBitsOnLeftOfEdge.end())
            bst.copy(temp_bitset, iterA->second);
                
        // union with B
        auto gotL1 = map_windows_by_right_corner.find (lev_i_1);
        if ( gotL1 != map_windows_by_left_corner.end() )
        {
            vector<int> intervals = gotL1->second;
            for(int w_ : intervals)
                setBitTo1(w_, temp_bitset);
        }

        winBitsOnLeftOfEdge[lev_i] = temp_bitset;
    }
}







void WindowContainer::makeOverlaps_with_intervals()
{
    if(hasMadeOverlaps)
    {
        cerr << "Already computed overlaps" << endl;
        return;
    }

    cout << "Max RNA Len = " << maxRNALen << endl;
    auto begin = chrono::high_resolution_clock::now();    

    numOfWins   = allWins.size();

    // create tools for this universe size
    mrnai_tools::BitSetTools bst(numOfWins);
    bitVecLongs = bst.num_words();
    
    // we will reset this everytime
    winThatCreateCycles = bst.makeBitSet();

    // Create set of windows per level
    for(int i=0; i<numOfLevels; i++)
        levelMapper[i] = bst.makeBitSet();

    // create mappings from Corner point to list of windows that have that
    // point as start or end
    map_corners_to_windows(&allWins);

    vector<pair_bitmap> winBitsOnRightOfPeg(numEven + numOdd);
    vector<pair_bitmap> winBitsOnLeftOfPeg(numEven + numOdd);

    for(int num=0; num<numEven+numOdd; num++)
    {
        sweepRightToLeft(num, winBitsOnRightOfPeg[num]);
        sweepLeftToRight(num, winBitsOnLeftOfPeg[num]);
    }

    // levelMapper needs to be filled in before main loop, not in it
    // due to bst.minus(leftright, levelMapper[level_id]);
    int counter = 0;
    for(const Window& win : allWins)
    {
        int even = win.l1;
        int odd  = win.l2;
    
        int level_id = level_pair_to_id[INDEX(even, odd, numOfRNAs)];    // [even][odd];
        setBitTo1(counter++, levelMapper[level_id]);
    }

    counter = 0;

    // create once, then clean and reuse everytime, b/c even if we delete 
    // and realloc everytime, we still have to set set every byte to zero
    // so better off only resetting to 0 instead of reallocing
    mrnai_tools::BitSetTools::Container tempNonOverlapArray_left_A  = bst.makeBitSet();
    mrnai_tools::BitSetTools::Container tempNonOverlapArray_right_A = bst.makeBitSet();
    mrnai_tools::BitSetTools::Container tempNonOverlapArray_left_B  = bst.makeBitSet();
    mrnai_tools::BitSetTools::Container tempNonOverlapArray_right_B = bst.makeBitSet();

    // for windows with one interval on level != l1,l2
    // reuse array.
    mrnai_tools::BitSetTools::Container leftright = bst.makeBitSet();   


    for(const Window& win : allWins)
    {
        int even = win.l1;
        int odd  = win.l2;
    
        mrnai_tools::BitSetTools::Container tempOverlapArray    = bst.makeBitSet();
        mrnai_tools::BitSetTools::Container tempNonOverlapArray = bst.makeBitSet();
        
        // Reset bit arrays so that we may use them again 
        // without reallocating
        bst.reset(tempNonOverlapArray);
        bst.reset(tempNonOverlapArray_left_A);
        bst.reset(tempNonOverlapArray_right_A);
        bst.reset(tempNonOverlapArray_left_B);
        bst.reset(tempNonOverlapArray_right_B);

        int even_right_corner = win.i;
        int even_left_corner  = win.i - win.w1;
        int odd_right_corner  = win.j;
        int odd_left_corner   = win.j - win.w2;

        pair<int, int> even_right = std::make_pair(even, even_right_corner);
        pair<int, int> even_left  = std::make_pair(even, even_left_corner);

        pair<int, int> odd_right  = std::make_pair(odd,  odd_right_corner);
        pair<int, int> odd_left   = std::make_pair(odd,  odd_left_corner);
        
        pair<int, int> even_right_plus1 = std::make_pair(even, even_right_corner+1);
        pair<int, int> even_left_minus1 = std::make_pair(even, even_left_corner-1);

        pair<int, int> odd_right_plus1  = std::make_pair(odd,  odd_right_corner+1);
        pair<int, int> odd_left_minus1  = std::make_pair(odd,  odd_left_corner-1);



        pair_bitmap::const_iterator gotL0 = winBitsOnLeftOfPeg [even].find(even_left);
        pair_bitmap::const_iterator gotR0 = winBitsOnRightOfPeg[even].find(even_right);

        pair_bitmap::const_iterator gotL1 = winBitsOnLeftOfPeg [odd].find(odd_left);
        pair_bitmap::const_iterator gotR1 = winBitsOnRightOfPeg[odd].find(odd_right);

    
        pair_bitmap::const_iterator gotL0_minus1 = winBitsOnLeftOfPeg [even].find(even_left_minus1);
        pair_bitmap::const_iterator gotR0_plus1  = winBitsOnRightOfPeg[even].find(even_right_plus1);

        pair_bitmap::const_iterator gotL1_minus1 = winBitsOnLeftOfPeg [odd].find(odd_left_minus1);
        pair_bitmap::const_iterator gotR1_plus1  = winBitsOnRightOfPeg[odd].find(odd_right_plus1);


        // Non overlaps: union of 4 sets
        // where each set is the intersection of two sets of windows

        // windows on left: gap on bottom
        if(gotL0 != winBitsOnLeftOfPeg[even].end())
            bst.copy(tempNonOverlapArray_left_A, gotL0->second);

        if(gotL1_minus1 != winBitsOnLeftOfPeg[odd].end())
            bst.s_intersection(tempNonOverlapArray_left_A, gotL1_minus1->second);
        else
            bst.reset(tempNonOverlapArray_left_A);


        // windows on left: gap on top
        if(gotL0_minus1 != winBitsOnLeftOfPeg[even].end())
            bst.copy(tempNonOverlapArray_left_B, (gotL0_minus1->second));
                
        if(gotL1 != winBitsOnLeftOfPeg[odd].end())
            bst.s_intersection(tempNonOverlapArray_left_B, (gotL1->second));
        else
            bst.reset(tempNonOverlapArray_left_B);


        // windows on right: gap on bottom
        if(gotR0 != winBitsOnRightOfPeg[even].end())
            bst.copy(tempNonOverlapArray_right_A, (gotR0->second));
        
        if(gotR1_plus1 != winBitsOnRightOfPeg[odd].end())
            bst.s_intersection(tempNonOverlapArray_right_A, (gotR1_plus1->second));
        else
            bst.reset(tempNonOverlapArray_right_A);


        // windows on right: gap on top
        if(gotR0_plus1 != winBitsOnRightOfPeg[even].end())
            bst.copy(tempNonOverlapArray_right_B, (gotR0_plus1->second));
        
        if(gotR1 != winBitsOnRightOfPeg[odd].end())
            bst.s_intersection(tempNonOverlapArray_right_B, (gotR1->second));
        else
            bst.reset(tempNonOverlapArray_right_B);

        // Union of 4 sets...
        bst.copy   (tempNonOverlapArray, tempNonOverlapArray_right_A);
        bst.s_union(tempNonOverlapArray, tempNonOverlapArray_right_B);
        bst.s_union(tempNonOverlapArray, tempNonOverlapArray_left_A );
        bst.s_union(tempNonOverlapArray, tempNonOverlapArray_left_B );
        

        // Now for all other levels:
        // get non overlapping itervals on both even and odd, and 
        // remove windows that are in current level
        bst.reset(leftright);
        

        int level_id = level_pair_to_id[INDEX(win.l1, win.l2, numOfRNAs)];    // [win.l1][win.l2];
            
        if(gotL0 != winBitsOnLeftOfPeg[even].end())
            bst.s_union(leftright, gotL0->second);
        
        if(gotR0 != winBitsOnRightOfPeg[even].end())
            bst.s_union(leftright, gotR0->second);

        if(gotL1 != winBitsOnLeftOfPeg[odd].end())
            bst.s_union(leftright, gotL1->second);
        
        if(gotR1 != winBitsOnRightOfPeg[odd].end())
            bst.s_union(leftright, gotR1->second);

        bst.minus(leftright, levelMapper[level_id]);
        bst.s_union(tempNonOverlapArray, leftright);
    

        // when both levels are different, we have no overlap
        for(int e_even : even_levels)
        {
            if(e_even == even) continue;
            for(int o_odd : odd_levels)
            {
                if(o_odd == odd) continue;
                
                int level_id = level_pair_to_id[INDEX(e_even, o_odd, numOfRNAs)]; // [e_even][o_odd];
                
                bst.s_union(tempNonOverlapArray, levelMapper[level_id]);
            }
        }   
    

        bst.complement(tempOverlapArray, tempNonOverlapArray);

        setBitTo0(counter, tempOverlapArray);

        string winStr = win.toString();
        nonOverlapsHM_Long[winStr] = tempNonOverlapArray;
        overlapsHM_Long[winStr] = tempOverlapArray;

        counter++;
    }

    // Clean up
    bst.deleteBitSet(tempNonOverlapArray_left_A);
    bst.deleteBitSet(tempNonOverlapArray_right_A);
    bst.deleteBitSet(tempNonOverlapArray_left_B);
    bst.deleteBitSet(tempNonOverlapArray_right_B);
    bst.deleteBitSet(leftright);

    for(int num=0; num<numEven+numOdd; num++)
    {
        for(const auto& kv : winBitsOnLeftOfPeg[num]) 
            bst.deleteBitSet(kv.second);

        for(const auto& kv : winBitsOnRightOfPeg[num]) 
            bst.deleteBitSet(kv.second);
    }

    auto end = chrono::high_resolution_clock::now();    
    auto dur = end - begin;
    auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(dur).count();
    cout << "Time to make overlaps (intervals):" << (double)ms << "ms" << endl;


    // ofstream outfile("interval_out");
    // for(std::vector<Window>::iterator it = vec->begin(); it != vec->end(); ++it, counter++)
    // {
    //     Window win = *it;
    //     string winStr = win.toString();
    //     for(int k=0;k<bitVecLongs;k++) 
    //     {
    //         std::bitset<64> bb(nonOverlapsHM_Long[winStr][k]);
    //         outfile << bb;
    //         // outfile << nonOverlapsHM_Long[winStr] << " ";
    //         // outfile << (unsigned long long)nonOverlapsHM_Long[winStr] << " ";
    //     }
    //     outfile << "\n";
    // }
    // outfile.close();
    // cout << "allWins.size() = " << allWins.size() << endl;
}








/* Return rules:
If returning 'add' or 'replace', return 0-based index of window that will be added from allWins.
For 'replace', set replacing to 0-based index of window that will be replaced in S
If deleting, return -1*(1-based index) of window that will be deleted from S
*/
int WindowContainer::sampleAddOrReplace(const vector<Window> & vec, int level1, int level2, 
                    int level_pair_id, int & replacing, int & totalSize, bool & staying)
{

    vector<mrnai_tools::BitSetTools::BIT_STORE*> bitVecs, nolBitVecs;
    for(auto elem : vec) 
    {
        string s = elem.toString();
        bitVecs.push_back(overlapsHM_Long[s]) ;
        nolBitVecs.push_back(nonOverlapsHM_Long[s]) ;
    }


    // Add window: get all wins that don't overlap vec
    int c = 0, choice = 99, windowNum = 0;
    
    for(int i=0;i<bitVecLongs;i++) 
    {
        mrnai_tools::BitSetTools::BIT_STORE x = -1;

        // make sure no overlaps with all windows in the current set:
        for(auto word : nolBitVecs)
            x = x & word[i];

        // reject windows that create a cycle
        // this bit is 1 if the window creates a cycle
        x = x & ~winThatCreateCycles[i];

        // &ing here set those window bits to 0 which are not on this level
        // thereby preventing us from choosing them
        x = x & levelMapper[level_pair_id][i];

        mrnai_tools::BitSetTools::BIT_STORE one = -1; 
        one = ~(one >> 1);
        for(; one; one >>= 1, windowNum++)
            if(one & x) 
                if(rand() % ++c == 0) 
                    choice = windowNum;
    }

    // Replace window
    replacing = -1;
    unsigned int sz = vec.size();//, addReplCount = 0;
    for(int wi=0; wi<sz; wi++) 
    {
        if(vec[wi].l1 != level1 or vec[wi].l2 != level2)
            continue;

        windowNum = 0;
        for(int i=0;i<bitVecLongs;i++) 
        {
            // Initially x contains 1 for those windows that overlap window-to-be-replaced
            // so all 1's here represent candidate windows (that will replace)
            mrnai_tools::BitSetTools::BIT_STORE x = bitVecs[wi][i];

            // Now & with bitvector of no-overlaps, so if a candidate overlaps existing
            // windows in the set, we set it to 0
            for(int j=0; j<nolBitVecs.size(); j++)
                if(j != wi)
                    x = x & nolBitVecs[j][i];
            
            // &ing here set those candidate window bits to 0 which are not on this level
            x = x & levelMapper[level_pair_id][i];

            // reject windows that create a cycle
            // this bit is 1 if the window creates a cycle
            x = x & ~winThatCreateCycles[i];
            
            mrnai_tools::BitSetTools::BIT_STORE one = -1; 
            one = ~(one >> 1);
            for(; one; one >>= 1, windowNum++)
                if(one & x) 
                    if(rand() % ++c == 0) 
                    {
                        choice = windowNum; 
                        replacing = wi;
                    }
        }
    }

    // Sample neighbors via deleted
    int delC = -1;
    for(int i=0;i<sz;i++)
    {
        if(vec[i].l1 != level1 or vec[i].l2 != level2)
            continue;
    
        if(rand() % ++c == 0) 
        {
            choice = c;
            delC = i;
        }
    }
    
    if(delC >= 0) 
    {
        choice = -1*delC - 1;
        replacing = -1;
    }

    // randomly decide to stay
    staying = false;
    if(rand() % ++c == 0) 
        staying = true;

    totalSize = c;

    return choice;
}



int WindowContainer::getNbrSizeOnly(const vector<Window> & vec, int level1, int level2, int level_pair_id)
{
    vector<mrnai_tools::BitSetTools::BIT_STORE*> bitVecs, nolBitVecs;
    for(auto elem : vec) 
    {
        string s = elem.toString();
        bitVecs.push_back(overlapsHM_Long[s]) ;
        nolBitVecs.push_back(nonOverlapsHM_Long[s]) ;
    }

    // Add window: get all wins that don't nonoverlap vec
    int sum = 0;
    for(int i=0;i<bitVecLongs;i++) 
    {
        mrnai_tools::BitSetTools::BIT_STORE x = nolBitVecs.size() == 0 ? 0 : -1;
        for(auto word : nolBitVecs)
            x = x & word[i];
    
        x = x & levelMapper[level_pair_id][i];

        // reject windows that create a cycle
        // this bit is 1 if the window creates a cycle
        x = x & ~winThatCreateCycles[i];

        
        for(; x; x >>= 1)
            sum += x & 1;
    }
    
    // Replace window
    unsigned int sz = vec.size();
    for(int wi=0; wi<sz; wi++) 
    {
        if(vec[wi].l1 != level1 or vec[wi].l2 != level2)
            continue;

        for(int i=0;i<bitVecLongs;i++) 
        {
            mrnai_tools::BitSetTools::BIT_STORE x = bitVecs[wi][i];

            for(int j=0; j<nolBitVecs.size(); j++)
                if(j != wi)
                    x = x & nolBitVecs[j][i];

            x = x & levelMapper[level_pair_id][i];

            // reject windows that create a cycle
            // this bit is 1 if the window creates a cycle
            x = x & ~winThatCreateCycles[i];
            
            for(; x; x >>= 1)
                sum += x & 1;
        }
    }
    

    // For neighbors via delete
    for(Window win : vec)
        if(win.l1 == level1 and win.l2 == level2)
            sum++;
    
    // for staying
    sum += 1;
    
    return sum;
}



const vector<Window> * WindowContainer::getLevel(int i)
{
    return windowsByLevel[i];
}





int WindowContainer::getNumOfLevels()
{
    return numOfLevels;
}



unordered_map<string, vector<Region>> WindowContainer::createGraphFromConfig(Config S)
{
    if(DEBUG) cout << endl << "*** create graph for " << S << endl << endl;
    vector<Window> vec = S.getVector();
    unordered_map<string, vector<Region>> adjList;
    vector<vector<Region>> intervalsByLevel(numOfRNAs);
    
    for(auto win : vec)
    {
        int odd  = win.l2;
        int even = win.l1;

        Region interval1(even, win.i, win.w1);
        Region interval2(odd,  win.j, win.w2);

        adjList[interval1.toString()].push_back(interval2);
        adjList[interval2.toString()].push_back(interval1);

        intervalsByLevel[even].push_back(interval1);
        intervalsByLevel[odd].push_back(interval2);
    }

    for(int i=0; i<numOfLevels;i++)
    {
        std::sort(intervalsByLevel[i].begin(), intervalsByLevel[i].end(), regionComparator); 
    }

    
    for(int i=0; i<numOfLevels;i++)
    {
        if(intervalsByLevel[i].size() >= 2)
        {
            for(int j=0; j<intervalsByLevel[i].size()-1; j++)
            {
                Region p = intervalsByLevel[i][j];
                Region n = intervalsByLevel[i][j+1];

                adjList[p.toString()].push_back(n);
            }
        }
    }

    // Print adjacency list

    // cout << endl << "adjacency list: " << endl;
    // for(auto it : adjList)
    // {
    //  cout << it.first << " : ";
    //  for(auto win : it.second)
    //      cout << win << "  ";
    //  cout << endl;
    // }

    return adjList;
}


unordered_set<string> WindowContainer::simpleBFS(string start, unordered_map<string, vector<Region>>&  adjList)
{
    unordered_set<string> visited;
    queue<string> Q;

    Q.push(start);

    while(Q.size() > 0)
    {
        string v = Q.front();
        Q.pop();
        visited.insert(v);

        for(Region nbr : adjList[v])
        {
            if(visited.find(nbr.toString()) == visited.end())
            {
                Q.push(nbr.toString());
            }
        }
    }

    return visited;
}

unordered_map<pair<string, string>, bool, pairhash> WindowContainer::
    findAllPaths(Config S, unordered_map<string, vector<Region>>&  adjList)
{
    unordered_map<pair<string, string>, bool, pairhash> isPath;

    for(auto it : adjList)
    {
        string v = it.first;
        unordered_set<string> pathsTo = simpleBFS(v, adjList);
        for(auto v2 : pathsTo)
        {
            isPath[make_pair(v, v2)] = true;
        }
    }

    return isPath;
}


template<typename T>
void printVec(vector<T> & v)
{
    for(auto e : v)
        cout << e << "  ";
    cout << endl;
}


bool WindowContainer::addWinAndTest(Config S, 
            unordered_map<string, vector<Region>>&  adjList, 
            Window win, 
            unordered_map<pair<string, string>, bool, pairhash> & isPath)
{
    int odd  = win.l2;
    int even = win.l1;

    // Each undirected edge representing a window is (X,Y), 
    // where X is the even interval and Y is the odd interval
    Region X(even, win.i, win.w1);
    Region Y(odd, win.j, win.w2); 

    vector<Region> vec {X, Y};

    string beforeXY[2], afterXY[2];

    for(int xy=0; xy<2; xy++)
    {
        Region interval = vec[xy];
        int level = interval.l;
        int i_    = interval.i;

        vector<Region> before;
        vector<Region> after;

        for(auto it : adjList)
        {
            Region r(it.first);

            if(level == r.l and r.i < i_)
                before.push_back(r);
            else if(level == r.l and r.i > i_)
                after.push_back(r);
        }

        sort(before.begin(), before.end(), regionComparator);
        sort(after.begin(),  after.end(),  regionComparator);

        if(before.size() == 0)
            beforeXY[xy] = "";
        else
            beforeXY[xy] = before.back().toString();

        if(after.size() == 0)
            afterXY[xy] = "";
        else 
            afterXY[xy] = after[0].toString();
    }

    string beforeX = beforeXY[0], 
           afterX  = afterXY[0], 
           beforeY = beforeXY[1], 
           afterY  = afterXY[1];    

    // case 1: ->X , X-Y, Y->
    string Xnbr = beforeX;
    string Ynbr = afterY;


    bool a = false;
    if(Xnbr != "" and Ynbr != "")
        a = isPath[make_pair(Ynbr, Xnbr)];
    
    // case 2: X->, X-Y, ->Y
    Xnbr = afterX;
    Ynbr = beforeY;
    bool b = false;
    if(Xnbr != "" and Ynbr != "")
        b = isPath[make_pair(Xnbr, Ynbr)];
    
    return !(a or b);
}


void WindowContainer::clearWindowCycleBits()
{
    for(int k=0; k<bitVecLongs; k++) 
        winThatCreateCycles[k] = 0;
}

void WindowContainer::setWindowCycleBit(Window win)
{
    setBitTo1(win.id, winThatCreateCycles);
}
