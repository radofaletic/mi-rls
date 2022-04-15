/*
 contains the Angle namespace for dealing with angle conversions
 
 Rado Faletic
 Department of Physics
 Faculty of Science
 Australian National University  ACT  0200
 Australia
 
 Rado.Faletic@anu.edu.au
 18th April 2004
 */


#ifndef _ANGLES_
#define _ANGLES_


/* ---------- standard header files ---------- */
#include <cmath>
#include <valarray>
/* ------------ user header files ------------ */
#include "extra_math.h"
/* ------------------------------------------- */

/* ------------------------------------------ */
/* ---------- function definitions ---------- */
/* ------------------------------------------ */


namespace Angle
{
// the number PI
// 3.141592653589793238462643383279502884197169399375105820974944592
const long double pi = (long double)(4) * std::atan((long double)(1));

// convert radians to degrees
template<class T> inline T r2d(const T& radians) { return ( is_zero(radians) ) ? T(0) : radians * T( (long double)(180) / Angle::pi ); };

// convert degrees to radians
template<class T> inline T d2r(const T& degrees) { return ( is_zero(degrees) ) ? T(0) : degrees * T( Angle::pi / (long double)(180) ); };

// convert Cartesian to polar coordinates (x,y) to (r,theta)
template<class T> inline void C2p(std::valarray<T>& point)
{
	T x = point[0]; T y = point[1]; point[0] = std::sqrt(x*x+y*y);
	if ( is_zero(x) && is_zero(y) ) { point[1] = T(0); }
	else if ( is_zero(x) && y < T(0) ) { point[1] = T(270); }
	else if ( is_zero(x) && T(0) < y ) { point[1] = T(90); }
	else if ( x < T(0) && is_zero(y) ) { point[1] = T(180); }
	else if ( T(0) < x && is_zero(y) ) { point[1] = T(0); }
	else
	{
		point[1] = std::atan2(std::abs(y),std::abs(x)); point[1] = Angle::r2d(point[1]);
		if ( x < T(0) && y < T(0) ) { point[1] += T(180); }
		else if ( x < T(0) ) { point[1] = T(180) - point[1]; }
		else if ( y < T(0) ) { point[1] = T(360) - point[1]; };
	};
};

enum axes { X, Y, Z, XY, YZ, ZX, XYZ };
}


#endif /* _ANGLES_ */
