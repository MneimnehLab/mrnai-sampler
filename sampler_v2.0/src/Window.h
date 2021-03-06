#ifndef _WINDOW_H_
#define _WINDOW_H_

#include <iostream>
#include <vector>
#include <string>
#include <stdexcept>

using namespace std;

struct Window
{
  int l1,l2,i,j,w1,w2;
  double weight;
  double intrEnergy;

  int id;

  // static int numRNAs;

  Window(int, int, int, int, int, int) ;
  Window(int, int, int, int, int, int, double, double) ;

  Window(std::string) ;
  
  bool operator< (Window &I2);
  bool operator== (Window &I2);
  bool operator!= (Window &I2); 
  friend bool operator== (const Window& wA, const Window& wB);
  friend bool operator!= (const Window& wA, const Window& wB);
  friend bool operator< (const Window& wA, const Window& wB);

  bool overlaps(Window);
  string toString() const;
  string toStringW() const;

  static std::vector<Window>* createWindowSetFromPyString(std::string line);

  static bool lessThan(const Window&, const Window&);

  static void setNumEvenOdd(int e, int o)
  {
    numEven = e;
    numOdd = o;
  }
private:
  static int numEven;
  static int numOdd;
  
};

std::ostream& operator<<(std::ostream& os, const Window& I);

namespace std {
template <>
  struct hash<Window>
  {
    std::size_t operator()(const Window& k) const
    {
      using std::size_t;
      using std::hash;
      using std::string;

      // Compute individual hash values for first,
      // second and third and combine them using XOR
      // and bit shifting:
      cout << "here = " << k << endl;
      cout << (hash<string>()(k.toString())) << endl;
      return (hash<string>()(k.toString()));
    }
  };
}

#endif