/**
 segment
 
 A collection of C++ routines for dealing with line segments.
 
 Rado Faletic
 Department of Physics
 Faculty of Science
 Australian National University  ACT  0200
 Australia
 
 Rado.Faletic@anu.edu.au
 21st June 2004
 22nd April 2022
 */





#ifndef _SEGMENT_
#define _SEGMENT_





/* ---------- standard header files ---------- */
#include <algorithm>
#include <valarray>





/* ------- user header files -------- */
#include "extra_math.h"
#include "line.h"
#include "shape.h"





/* ---------- class & function declaration ---------- */

// `segment' class
template<class T> class segment : public compact_shape<T>
{
private:
	line<T> container_;
public:
	void init(const std::valarray<T>*, const std::valarray<T>*);
	segment(const std::valarray<T>*, const std::valarray<T>*);
	segment() { };
	T make_measure() const;
	T scale() const;
	const std::valarray<T>* gradient_p() const;
	T distance_p(const std::valarray<T>*) const;
	bool contains_p(const std::valarray<T>*) const;
	T intersect(const line<T>&, std::valarray<T>&, bool&) const;
	const compact_shape<T>* face(const unsigned short&) const;
};





/* ---------- function definitions ---------- */





/* ---------- init ---------- */
template<class T> void segment<T>::init(const std::valarray<T>* pt1, const std::valarray<T>* pt2)
{
	// allocate segment properties
	this->assign_nodes(2, pt1, pt2);
	
	this->container_.init(this->p(0), this->p(1));
	//this->measure_ = this->make_measure();
}





/* ---------- segment ---------- */
template<class T> inline segment<T>::segment(const std::valarray<T>* pt1, const std::valarray<T>* pt2)
{
	this->init(pt1, pt2);
}





/* ---------- make_measure ---------- */
template<class T> inline T segment<T>::make_measure() const
{
	return norm(this->p(0), this->p(1));
}





/* ---------- scale ---------- */
template<class T> inline T segment<T>::scale() const
{
	return this->measure();
}





/* ---------- gradient_p ---------- */
template<class T> inline const std::valarray<T>* segment<T>::gradient_p() const
{
	return this->container_.gradient_p();
}





/* ---------- distance_p ---------- */
template<class T> T segment<T>::distance_p(const std::valarray<T>* pt) const
{
	if ( this->contains(this->container_.closest(*pt)) )
	{
		return this->container_.distance_p(pt);
	}
	return std::min(norm(this->p(0), pt), norm(this->p(1), pt));
}





/* ---------- contains_p ---------- */
template<class T> bool segment<T>::contains_p(const std::valarray<T>* pt) const
{
	if ( this->dim() != pt->size() )
	{
		throw;
        return false;
	}
	else if ( v_eq(pt, this->p(0)) || v_eq(pt, this->p(1)) )
	{
		return true;
	}
	for (unsigned short i=0; i<this->dim(); i++)
	{
		if ( (*pt)[i] < std::min((*(this->p(0)))[i],(*(this->p(1)))[i]) ||
			std::max((*(this->p(0)))[i],(*(this->p(1)))[i]) < (*pt)[i] )
		{
			return false;
		}
	}
	// see if the point lies in the line
	if ( !this->container_.contains_p(pt) )
	{
		return false;
	}
	return true;
}





/* ---------- intersect ---------- */
template<class T> T segment<T>::intersect(const line<T>& line2, std::valarray<T>& cut, bool& parallel) const
{
	cut.resize(0);
	parallel = false;
	if ( this->dim() != line2.dim() )
	{
		throw;
        return T(0);
	}
	// see if the segment lies in the line
	else if ( v_pm_eq(line2.gradient_p(), this->container_.gradient_p()) )
	{
		if ( line2.contains_p(this->p(0)) )
		{
			cut.resize(2 * (this->dim()));
			cut[std::slice(0, this->dim(), 1)] = *(this->p(0));
			cut[std::slice(this->dim(), this->dim(), 1)] = *(this->p(1));
			parallel = true;
			return this->measure();
		}
	}
	// check the intersection properties
	if ( !this->container_.intersect(line2, cut, parallel) )
	{
		cut.resize(0);
		parallel = false;
		return T(0);
	}
	parallel = false;
	for (unsigned short i=0; i<this->dim(); i++)
	{
		if ( cut[i] < std::min((*(this->p(0)))[i],(*(this->p(1)))[i]) ||
			std::max((*(this->p(0)))[i],(*(this->p(1)))[i]) < cut[i] )
		{
			cut.resize(0);
			return T(0);
		}
	}
	return T(1);
}





/* ---------- face ---------- */
template<class T> inline const compact_shape<T>* segment<T>::face(const unsigned short& en) const
{
	return 0;
}





#endif /* _SEGMENT_ */
