/*
 plane
 
 A collection of C++ routines for dealing with planes.
 
 Rado Faletic
 Department of Physics
 Faculty of Science
 Australian National University  ACT  0200
 Australia
 
 Rado.Faletic@anu.edu.au
 4th May 2004
 */


#ifndef _PLANE_
#define _PLANE_


/* ---------- standard header files ---------- */
#include <cmath>
#include <sstream>
#include <string>
#include <valarray>
/* ------- user header files -------- */
#include "extra_math.h"
#include "line.h"
#include "shape.h"
/* ---------------------------------- */


/* --------------------------------------- */
/* ---------- class declaration ---------- */
/* --------------------------------------- */


// `plane' class
template<class T>
class plane : public shape<T>
{
private:
	std::valarray<T> zero_;
	std::valarray<T> normal_;
public:
	void init(const std::valarray<T>&, const std::valarray<T>&);
	void init(const std::valarray<T>*, const std::valarray<T>*, const std::valarray<T>*);
	plane(const std::valarray<T>&, const std::valarray<T>&);
	plane(const std::valarray<T>*, const std::valarray<T>*, const std::valarray<T>*);
	plane() { };
	unsigned short dim() const;
	std::string print() const;
	const std::valarray<T>* normal_p() const;
	std::valarray<T> zero() const;
	std::valarray<T> drop(const std::valarray<T>&) const;
	T distance_p(const std::valarray<T>*) const;
	std::valarray<T> closest(const std::valarray<T>&) const;
	bool contains_p(const std::valarray<T>*) const;
	T intersect(const line<T>&, std::valarray<T>&, bool&) const;
};


/* ------------------------------------------ */
/* ---------- function definitions ---------- */
/* ------------------------------------------ */


/* ---------- init ---------- */
template<class T> void
plane<T>::init(const std::valarray<T>& nm, const std::valarray<T>& ze)
{
	if ( nm.size() != ze.size() || nm.size() < 2 )
	{
		throw; return;
	}
	
	// allocate plane properties
	this->normal_.resize(nm.size(), T(0));
	this->zero_.resize(ze.size(), T(0));
	if ( nm.size() == 3 )
	{
		this->normal_ = nm;
		normalise(this->normal_);
		this->zero_ = ze;
	}
}
/* -------------------------- */

/* ---------- init ---------- */
template<class T> inline void
plane<T>::init(const std::valarray<T>* pt1, const std::valarray<T>* pt2, const std::valarray<T>* pt3)
{
	this->init(cross(pt1, pt2, pt3), *pt1);
}
/* -------------------------- */

/* ---------- plane ---------- */
template<class T> inline
plane<T>::plane(const std::valarray<T>& nm, const std::valarray<T>& ze)
{
	this->init(nm, ze);
}
/* --------------------------- */

/* ---------- plane ---------- */
template<class T> inline
plane<T>::plane(const std::valarray<T>* pt1, const std::valarray<T>* pt2, const std::valarray<T>* pt3)
{
	this->init(pt1, pt2, pt3);
}
/* --------------------------- */

/* ---------- dim ---------- */
template<class T> inline unsigned short
plane<T>::dim() const
{
	return this->zero_.size();
}
/* ------------------------- */

/* ---------- print ---------- */
template<class T> std::string
plane<T>::print() const
{
	std::ostringstream tmp;
	tmp << "< x - ";
	for (unsigned short i=0; i<(this->zero_).size(); i++)
	{
		if ( i )
		{
			tmp << ',';
		}
		tmp << (this->zero_)[i];
	}
	tmp << ", ";
	for (unsigned short i=0; i<(this->normal_).size(); i++)
	{
		if ( i )
		{
			tmp << ',';
		}
		tmp << (this->normal_)[i];
	}
	tmp << " > = 0";
	return tmp.str();
}
/* --------------------------- */

/* ---------- normal_p ---------- */
template<class T> inline const std::valarray<T>*
plane<T>::normal_p() const
{
	return &(this->normal_);
}
/* ------------------------------ */

/* ---------- zero ---------- */
template<class T> inline std::valarray<T>
plane<T>::zero() const
{
	return this->zero_;
}
/* -------------------------- */

/* ---------- drop ---------- */
template<class T> inline std::valarray<T>
plane<T>::drop(const std::valarray<T>& pt) const
{
	// drop = pt - c.normal_ (where c is an unknown constant)
	// solution: c = dot(pt - zero_, normal_) / norm2(normal_)
	//           but norm(normal_) = 1
	return pt - ( dot(std::valarray<T>(pt - this->zero_), this->normal_) * this->normal_ );
}
/* -------------------------- */

/* ---------- distance_p ---------- */
template<class T> inline T
plane<T>::distance_p(const std::valarray<T>* pt) const
{
	// Theorem 3.5.2, Anton pg 142
	T d = std::abs( dot(this->normal_, *pt) - dot(this->normal_, this->zero_) );
	return ( is_zero(d) ) ? T(0) : d;
}
/* -------------------------------- */

/* ---------- closest ---------- */
template<class T> inline std::valarray<T>
plane<T>::closest(const std::valarray<T>& pt) const
{
	return this->drop(pt);
}
/* ----------------------------- */

/* ---------- contains_p ---------- */
template<class T> inline bool
plane<T>::contains_p(const std::valarray<T>* pt) const
{
	return ( this->dim() < 3 || is_zero(this->distance_p(pt)) ) ? true : false;
}
/* -------------------------------- */

/* ---------- intersect ---------- */
template<class T> T
plane<T>::intersect(const line<T>& l, std::valarray<T>& intercept, bool& parallel) const
{
	intercept.resize(0);
	parallel = false;
	if ( this->dim() != l.dim() )
	{
		throw; return T(0);
	}
	else if ( this->dim() < 3 ) // a line always intersects a plane in 2d
	{
		parallel = true;
		return T(1);
	}
	
	// see if the line is parallel to the plane
	T dnl = dot(this->normal_, *(l.gradient_p()));
	if ( is_zero(dnl) )
	{
		if ( this->contains(l.zero()) )
		{
			parallel = true;
			return T(1);
		}
		return T(0);
	}
	
	// solve for the intercept
	//T t = ( dot(zero_, normal_) - dot(l.zero(), normal_) ) / dot(*(l.gradient_p()), normal_);
	//intercept = l.evaluate(t);
	intercept.resize(this->dim());
	intercept = l.evaluate( ( dot(this->zero_, this->normal_) - dot(l.zero(), this->normal_) ) / dnl );
	return T(1);
}
/* ------------------------------- */


#endif /* _PLANE_ */
