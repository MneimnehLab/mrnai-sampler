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
#include "BitSetTools.h"
using namespace std;


int main(int argc, char *argv[])
{
    mrnai_tools::BitSetTools bst(100);

    mrnai_tools::BitSetTools::Container blah = bst.makeBitSet();
    bst.setBitTo1(62, blah);

    bst.printInts(blah);



    return 0;
}
