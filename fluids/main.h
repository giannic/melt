#ifndef MAIN_H_
#define MAIN_H_ 1

#include <string>
using std::string;

#include <iostream>
using std::cout;
using std::cerr;
using std::cin;
using std::endl;
using std::ostream;
using std::istream;
using std::ifstream;
using std::ofstream;

#include <iomanip>
using std::setprecision;
using std::setw;

#include <algorithm>
using std::min;
using std::max;
using std::swap;
using std::sort;
using std::fill;
using std::copy;

#include <vector>
using std::vector;
#include <list>
using std::list;
#include <set>
using std::set;
#include <map>
using std::map;
#include <deque>
using std::deque;

#include <cmath>
using std::sqrt;
using std::abs;
using std::exp;
using std::cos;
using std::sin;
using std::atan;
using std::atan2;
using std::log;
using std::exp;
using std::ceil;
using std::floor;

#include <ctime>
using std::clock;
using std::time;
using std::difftime;

#include <cstdlib>
#include <cassert>
#include <fstream>
using std::ostringstream;

#ifndef M_PI
#define M_PI 3.1415926535898
#endif

#include <sstream>
#include <iomanip>
#include <cstdio>

//------------------------------CONSTANTS---------------------------------------------
typedef double Float;
typedef double Double;
typedef unsigned int u_int;

#ifndef INFINITY
#define INFINITY HUGE_VAL
#endif

#ifndef HALF_PI
const Float HALF_PI = 2.0 * atan(1.0);
#endif
//const double M_PI = 4.0 * atan(1.0);
#ifndef TWO_PI
const Float TWO_PI = 8.0 * atan(1.0);
#endif

const Float INV_PI = 1.0/M_PI; 
const Float INV_TWO_PI = 1.0/TWO_PI;
const Float PI_INV180 = M_PI   / 180.;
const Float INV_PI180 = INV_PI * 180;

const Float ONE_THIRD = 1./3.;
const Float INV_255 = 1. / 255;
const Float SQRT_TWO = sqrt(2.);

#define HH					0.5       // CELL SIZE
#define PARTICLES_PER_NODE	32
#define NX                  10
#define NY                  100
#define NZ                  100

const Float MAX_U               = NX * 0.005;
const Float MAX_V               = NY * 0.005;
const Float MAX_W               = NZ * 0.005;
const Float RADIUS_MIN			= 0.1; // Don't Scale these by H
const Float RADIUS_MAX			= 0.5;
const Float RESEED_THRESHOLD	= 2.0 * HH;  
//const Float PARTICLE_DELETE		= 1.5 * RADIUS_MIN;
const Float PARTICLE_DELETE		= 100 * RADIUS_MIN;
const Float DT					= 4.9 / ((MAX_U + MAX_V + MAX_W) / HH);
const Float FASTMARCH_LIMIT		= 10.0 * HH; // extent of influence of fast marching
const Float SEMILAGRA_LIMIT		= 5.0 * HH; // extent of influence of semi-lagrangian
//const Float SEMILAGRA_LIMIT		= 100.0 * HH; // extent of influence of semi-lagrangian
const int PARTICLES_PER_INTERFACE_NODE = 2 * PARTICLES_PER_NODE;

#define SAMPLEPHI                  LinearSample      
//#define SAMPLEPHI                  CubicSample

#define FOR_LS     for(int k=1; k<=theDim[2]; k++) { for(int j=1; j<=theDim[1]; j++) { for(int i=1; i<=theDim[0]; i++) {
#define FOR_ALL_LS for(int k=0; k<(theDim[2]+2); k++) { for(int j=0; j<(theDim[1]+2); j++) { for(int i=0; i<(theDim[0]+2); i++) {
#define END_FOR_THREE }}}
#define END_FOR_TWO }}
#define FOR_GRID   for(int i=0; i < size; i++)
#define FOR_GRIDZY for(int k=0; k <= theDim[2]; k++) { for(int j=0; j <= theDim[1]; j++) { 
#define FOR_GRIDZX for(int k=0; k <= theDim[2]; k++) { for(int i=0; i <= theDim[0]; i++) { 
#define FOR_GRIDYX for(int j=0; j <= theDim[1]; j++) { for(int i=0; i <= theDim[0]; i++) { 

class FastMarch;
class Grid;
class LevelSet;
class Particle;
class ParticleSet;
class Timer;
class Vector;
class Velocity;
class Water;

inline bool CmpFtoZero(const Float &compare) {
    if((compare < 0.00001) && (compare > -0.00001)) return true;
    else return false;
}

inline Float Lerp(const Float &alpha, const Float &a, const Float &b) 
    { return (1.0 - alpha) * a + alpha * b; }

// Monotonic Cubic Interpolation
// interpolation between fk1 and fk2. alpha = position - k1;
inline Float MCerp(const Float &alpha, const Float &fk0, const Float &fk1, 
                                       const Float &fk2, const Float &fk3) 
{
    static Float dk02, dk13, delk;
    delk = fk2 - fk1;
    if(CmpFtoZero(delk)) dk02 = dk13 = 0;
    else {
        dk02 = 0.5 * (fk2 - fk0);
        dk13 = 0.5 * (fk3 - fk1);
        if(dk02 * delk < 0.) dk02 = 0;
        if(dk13 * delk < 0.) dk13 = 0;
    }
    return (dk02 + dk13 - 2 * delk) * alpha * alpha * alpha 
            + (3. * delk - 2. * dk02 - dk13) * alpha * alpha + dk02 * alpha + fk1;
}

inline int Step(const Float &x, const Float &alpha) 
    { if(alpha < x) return 0; return 1; }

inline int Round(const Float &alpha) { return (int)floor(alpha+0.5); }

inline Double square(const Double &x) { return x*x; }

inline Float Radians(const Float &deg) {return PI_INV180 * deg; }
inline Float Degrees(const Float &rad) {return INV_PI180 * rad; }

inline Float Clamp(const Float &alpha, const Float &a, const Float &b) { 
	if(alpha < a) return a;
	if(alpha > b) return b;
	return alpha;
}

inline int Clamp(const int &alpha, const int &a, const int &b) {
	if(alpha < a) return a;
	if(alpha > b) return b;
	return alpha;
}

inline Float SmoothStep(const Float &min, const Float &max, const Float &alpha) {
	Float temp = Clamp((alpha - min) / (max - min), 0.0, 1.0);
	return -2.0 * temp * temp * temp + 3.0 * temp * temp;
}

inline int Mod(int a, int b, const Double &bInv) {
	int n = int(a*bInv);
	a -= n*b;
	if (a < 0)
		a += b;
	return a;
}

//----------------need random code for this - else take it out------------------------------
extern Float genrand();
extern unsigned long genrandlong();

inline Float RandomFloat(Float min = 0.0, Float max = 1.0) {
	return Lerp(genrand(), min, max);
}

inline unsigned long RandomInt() {
	return genrandlong();
}

#endif