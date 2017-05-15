#include "BitSetTools.h"

namespace mrnai_tools
{

BitSetTools::BitSetTools(int universe_size) : universe_size(universe_size)
{
	array_size = universe_size / (sizeof(BitSetTools::BIT_STORE)*8) + 1;
}

inline BitSetTools::BIT_STORE* BitSetTools::makeBitSet()
{
    BIT_STORE* array = new BIT_STORE[array_size];
    for(int k=0; k<array_size; k++) 
        array[k] = 0;
	return array;
}

// Modifies A to (A intersection B)
inline void BitSetTools::s_intersection(BitSetTools::BIT_STORE * A, BitSetTools::BIT_STORE * B)
{
	for(int k=0; k<array_size; k++) 
        A[k] &= B[k];
}

// Modifies A to (A union B)
inline void BitSetTools::s_union(BitSetTools::BIT_STORE * A, BitSetTools::BIT_STORE * B)
{
	for(int k=0; k<array_size; k++) 
        A[k] |= B[k];
}

inline void BitSetTools::reset(BitSetTools::BIT_STORE * array)
{
	for(int k=0; k<array_size; k++) 
        array[k] = 0;
}


inline void BitSetTools::setBitTo1(int winNum, BitSetTools::BIT_STORE array[])
{
    int pos = winNum / (sizeof(BIT_STORE)*8);
    int idx = winNum % (sizeof(BIT_STORE)*8);

    // have to set a 1 at bit 'idx' (from right) in integer 'pos'
    BIT_STORE bits = array[pos];
    // BIT_STORE one = 1; one = one << idx;
    BIT_STORE mask = 0; 
    mask = ~((mask-1)>> 1);
    mask = mask >> idx;
    bits = bits | mask;
    array[pos] = bits;
}

inline void BitSetTools::setBitTo0(int winNum, BitSetTools::BIT_STORE array[])
{
    int pos = winNum / (sizeof(BIT_STORE)*8);
    int idx = winNum % (sizeof(BIT_STORE)*8);

    BIT_STORE bits = array[pos];
    BIT_STORE mask = -1; 
    mask = ~(mask>> 1);
    mask = mask >> idx;
    mask = ~mask;
    bits = bits & mask;
    array[pos] = bits;
}


}
