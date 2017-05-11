#include <vector>
#include <string>
#include <sstream>
#include "Region.h"

Region::Region(int l_, int i_, int w_) 
		: l(l_), i(i_), w(w_), upEnergy(0) {}

Region::Region(int l_, int i_, int w_, double upEnergy_) 
        : l(l_), i(i_), w(w_), upEnergy(upEnergy_) {}

Region::Region(std::string line) 
{
    for(int i=0;i<line.length();i++)
        if ( !isdigit(line.at(i)) )
            line[i] = ' ';

    std::stringstream ss;
    ss.str(line);
    
    ss >> l >> i >> w;
    

    if(ss.rdbuf()->in_avail() == 0)
        upEnergy = 0;
    else
        ss >> upEnergy;

}

bool Region::lessThan(const Region& r1, const Region& r2)
{
    //Order matters!!!
    if(r1.l < r2.l)
        return true;
    if(r1.l > r2.l)
        return false;
    
    if(r1.i < r2.i)
        return true;
    if(r1.i > r2.i)
        return false;

    if(r1.w < r2.w)
        return true;
    if(r1.w > r2.w)
        return false;
    
    return false;
}

bool Region::operator< (Region &win)
{
    return Region::lessThan(*this, win);
}

bool Region::equal(const Region& I1, const Region& I2)
{
    if(I1.l == I2.l && I1.i == I2.i && I1.w == I2.w )
        return true;
    
    return false;
}

bool Region::operator== (Region &I2)
{
    return Region::equal(*this, I2);
}

bool operator== (const Region& wA, const Region& wB)
{
    return Region::equal(wA, wB);
}

bool operator< (const Region& wA, const Region& wB)
{
    // return wA < wB;
    return Region::lessThan(wA, wB);
}

bool Region::operator!= (Region &I2)
{
    return !(*this == I2);
}

std::ostream& operator<<(std::ostream& os, const Region& I)
{
	os << "(" << I.l << ", " << I.i << ", " << I.w << ", " << ")";
	return os;
}

string Region::toString() const
{
    stringstream ss;
    ss.str("");
    ss << "(" << l << ", " << i << ", " << w << ", " << ")";
    return ss.str();
}

string Region::toStringW() const
{
    stringstream ss;
    ss.str("");
    ss << "(" << l << ", " << i << ", " << w << ", " << ") : " << upEnergy;
    return ss.str();
}


std::vector<Region>* Region::createRegionSetFromPyString(std::string line)
{
    std::vector<Region> * ints = new std::vector<Region>();

    for(int i=0;i<line.length();i++)
        if ( !isdigit(line.at(i)) )
            line[i] = ' ';

    std::stringstream ss;
    ss.str(line);
    int l,i,w;
    //int count = 0;
    while(ss >> l >> i >> w)
    {
        Region Region(l, i, w );

    }

    return ints;
}

bool Region::overlaps(Region win)
{
    return false;

    /*    
    int ATopLeft = i - w;
    int ATopRight = i;
    int ABotLeft = j - w2;
    int ABotRight = j;

    int BTopLeft = win.i - win.w;
    int BTopRight = win.i;
    int BBotLeft = win.j - win.w2;
    int BBotRight = win.j;

    if(l == win.l)
    {
        //return !( (ATopLeft > BTopRight and ABotLeft > BBotRight) or (ATopRight < BTopLeft and ABotRight < BBotLeft) );

        // We will consider it an over lap if both Regions are immediately adjacent, i.e., 
        // (ATopLeft == BTopRight+1 and ABotLeft == BBotRight+1) or (ATopRight == BTopLeft-1 and ABotRight == BBothLeft-1)
        if((ATopLeft == BTopRight+1 and ABotLeft == BBotRight+1) or (ATopRight == BTopLeft-1 and ABotRight == BBotLeft-1))
            return true;

        return !( (ATopLeft > BTopRight and ABotLeft > BBotRight) or (ATopRight < BTopLeft and ABotRight < BBotLeft) );
    }
    
    else if(l == win.l - 1)
        return !(ABotLeft > BTopRight or ABotRight < BTopLeft)  ;      
    
    else if(l == win.l + 1)
        return !(ATopLeft > BBotRight or ATopRight < BBotLeft) ;

    */

}

