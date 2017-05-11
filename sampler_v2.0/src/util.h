#ifndef __UTIL_H__
#define __UTIL_H__
#include <vector>
#include <string>

using std::vector;

class Util
{
public:
	template <class T>
	static vector<T> fastIntersection(const vector<vector<T>* > & list);
};

#include "util.cpp"
#endif