#ifndef __BITSET_H__
#define __BITSET_H__
#include <bitset>

namespace mrnai_tools
{
	
class BitSetTools
{
public:
	typedef uint64_t BIT_STORE;
	typedef BIT_STORE* Container;
private:
	size_t universe_size;
	size_t array_size;
public:
	BitSetTools(int universe_size) : universe_size(universe_size)
	{
		array_size = universe_size / (sizeof(BIT_STORE)*8) + 1;
	}

	inline BIT_STORE* makeBitSet()
	{
	    BIT_STORE* array = new BIT_STORE[array_size];
	    for(int k=0; k<array_size; k++) 
	        array[k] = 0;
		return array;
	}

	inline void deleteBitSet(BIT_STORE* array)
	{
	    delete[] array;
	}

	inline void copy(BIT_STORE * A, BIT_STORE * B)
	{
		for(int k=0; k<array_size; k++) 
	        A[k] = B[k];
	}

	// Modifies A to (A intersection B)
	inline void s_intersection(BIT_STORE * A, BIT_STORE * B)
	{
		for(int k=0; k<array_size; k++) 
	        A[k] &= B[k];
	}

	// Modifies A to (A union B)
	inline void s_union(BIT_STORE * A, BIT_STORE * B)
	{
		for(int k=0; k<array_size; k++) 
	        A[k] |= B[k];
	}

	inline void reset(BIT_STORE * array)
	{
		for(int k=0; k<array_size; k++) 
	        array[k] = 0;
	}

	// sets A = ~B
	inline void complement(BIT_STORE * A, BIT_STORE * B)
	{
		for(int k=0; k<array_size; k++) 
	        A[k] = ~B[k];
	}

	// Modifies A to be A - B
	inline void minus(Container A, Container B)
	{
		for(int k=0; k<array_size; k++) 
			A[k] &= ~B[k];
	}


	inline void setBitTo1(int winNum, BIT_STORE array[])
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

	inline void setBitTo0(int winNum, BIT_STORE array[])
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

	inline void printInts(BIT_STORE array[])
	{
		for(int k=0; k<array_size; k++) 
			std::cout << array[k] << " ";
		std::cout << std::endl;
	}

	inline size_t num_words()
	{
	    return array_size;
	}

	inline size_t get_universe_size()
	{
		return universe_size;
	}



};
};

#endif