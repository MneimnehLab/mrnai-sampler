#ifndef _REGION_H_
#define _REGION_H_

#include <iostream>
#include <vector>
#include <string>

using namespace std;

struct Region
{
	int l,i,w;
  double upEnergy;

  Region() {};
	Region(int, int, int) ;
  Region(int, int, int, double) ;

	Region(std::string) ;
	
	bool operator< (Region &I2);
	bool operator== (Region &I2);
	bool operator!= (Region &I2);	
	friend bool operator== (const Region& wA, const Region& wB);
  friend bool operator< (const Region& wA, const Region& wB);

	bool overlaps(Region);
	string toString() const;
  string toStringW() const;

	static std::vector<Region>* createRegionSetFromPyString(std::string line);
  static bool equal(const Region& r1, const Region& r2);
  static bool lessThan(const Region& r1, const Region& r2);
};

std::ostream& operator<<(std::ostream& os, const Region& I);

namespace std {
template <>
  struct hash<Region>
  {
    std::size_t operator()(const Region& k) const
    {
      using std::size_t;
      using std::hash;
      using std::string;

      // Compute individual hash values for first,
      // second and third and combine them using XOR
      // and bit shifting:
      //cout << "here = " << k << endl;
      cout << (hash<string>()(k.toString())) << endl;
      return (hash<string>()(k.toString()));
    }
  };
}

#endif