#include <chrono>
#include <iostream>
#include <string>
#include <sstream>
#include <vector>
#include <algorithm>
//#include "Window.h"
using namespace std;


template <class T>
vector<T> Util::fastIntersection(const vector<vector<T>* > & list)
{
	vector<T> result;
	vector<int> indices(list.size(), 0);
	int sameCount = 0;

	bool done = false;
	while(!done)
	{
		for(int i=0;i<indices.size();i++)
		{
			//cout << i << ",";
			if(indices[i] == (*list[i]).size())
			{
				done = true;
				break;
			}
		}
		//cout << endl;
		if(done)
			break;

		bool firstAreSame = true;
		T first = (*list[0])[indices[0]];
		for(int i=1;i<indices.size();i++)
		{
			T first_ = (*list[i])[indices[i]];
			if(first != first_)
			{
				firstAreSame = false;
				break;
			}
			else
				first = first_;

		}

		if(firstAreSame)
		{
			//cout << "Same = " << first << endl;
			sameCount++;
			result.push_back(first);

			for(int i=0;i<indices.size();i++)
			{
				indices[i] += 1;
			}
			continue;
		}

		T minVal = (*list[0])[indices[0]];
		int minList = 0;
		for(int i=1;i<indices.size();i++)
		{
			if((*list[i])[indices[i]] < minVal)
			{
				minVal = (*list[i])[indices[i]];
				minList = i;
			}
		}
		indices[minList]++;

	}

	//cout << "Same: " << sameCount << endl;

	return result;
}






