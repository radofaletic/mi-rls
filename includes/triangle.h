/**
 triangle
 
 A collection of C++ routines for dealing with generic triangles.
 
 Rado Faletic
 Department of Physics
 Faculty of Science
 Australian National University  ACT  0200
 Australia
 
 Rado.Faletic@anu.edu.au
 21st June 2004
 22nd April 2022
 */





#ifndef _TRIANGLE_
#define _TRIANGLE_





/* ---------- standard header files ---------- */
#include <algorithm>
#include <cmath>
#include <valarray>
#include <vector>





/* ------- user header files -------- */
#include "extra_math.h"
#include "line.h"
#include "plane.h"
#include "segment.h"
#include "shape.h"





/* ---------- class declaration ---------- */

// `triangle' class
template<class T> class triangle : public compact_shape<T>
{
private:
	plane<T> container_;
	segment<T> edge0_;
	segment<T> edge1_;
	segment<T> edge2_;
public:
	void init(const std::valarray<T>*, const std::valarray<T>*, const std::valarray<T>*);
	triangle(const std::valarray<T>*, const std::valarray<T>*, const std::valarray<T>*);
	triangle() { };
	T make_measure() const;
	T scale() const;
	T distance_p(const std::valarray<T>*) const;
	const std::valarray<T>* normal_p() const;
	bool contains_p(const std::valarray<T>*) const;
	bool contains_p(const std::valarray<T>*, const bool&) const;
	T intersect(const line<T>&, std::valarray<T>&, bool&) const;
	const compact_shape<T>* face(const unsigned short&) const;
};





/* ---------- function definitions ---------- */





/* ---------- init ---------- */
template<class T> void triangle<T>::init(const std::valarray<T>* pt1, const std::valarray<T>* pt2, const std::valarray<T>* pt3)
{
	if ( pt1->size() < 2 )
	{
		throw;
        return;
	}
	
	// allocate triangle nodes
	this->assign_nodes(3, pt1, pt2, pt3);
	
	this->container_.init(this->p(0), this->p(1), this->p(2));
	
	this->edge0_.init(this->p(1), this->p(2));
	this->edge1_.init(this->p(2), this->p(0));
	this->edge2_.init(this->p(0), this->p(1));
	
	// area
	//this->measure_ = this->make_measure();
}





/* ---------- triangle ---------- */
template<class T> inline triangle<T>::triangle(const std::valarray<T>* pt1, const std::valarray<T>* pt2, const std::valarray<T>* pt3)
{
	this->init(pt1, pt2, pt3);
}





/* ---------- make_measure ---------- */
template<class T> inline T triangle<T>::make_measure() const
{
	T a = norm(this->p(1), this->p(0));
	T b = norm(this->p(2), this->p(1));
	T c = norm(this->p(0), this->p(2));
	T s = ( a + b + c ) / T(2);
	return std::sqrt(s * ( s - a ) * ( s - b ) * ( s - c ));
}





/* ---------- scale ---------- */
template<class T> inline T triangle<T>::scale() const
{
	return ( this->edge0_.scale() + this->edge1_.scale() + this->edge2_.scale() ) / T(3);
}





/* ---------- distance_p ---------- */
template<class T> T triangle<T>::distance_p(const std::valarray<T>* pt) const
{
	if ( this->contains(this->container_.drop(*pt)) )
	{
		return this->container_.distance_p(pt);
	}
	T dist = this->edge0_.distance_p(pt);
	dist = std::min(dist, this->edge1_.distance_p(pt));
	dist = std::min(dist, this->edge2_.distance_p(pt));
	return dist;
}





/* ---------- normal_p ---------- */
template<class T> inline const std::valarray<T>* triangle<T>::normal_p() const
{
	return this->container_.normal_p();
}





/* ---------- contains_p ---------- */
template<class T> inline bool triangle<T>::contains_p(const std::valarray<T>* pt) const
{
	return this->contains_p(pt, false);
}





/* ---------- contains_p ---------- */
template<class T> bool triangle<T>::contains_p(const std::valarray<T>* pt, const bool& internal) const
{
	if ( this->dim() != pt->size() )
	{
		throw;
        return false;
	}
	// see if the point is contained within the plane of the triangle
	else if ( !internal && !this->container_.contains_p(pt) )
	{
		return false;
	}
	// see if the point lies on any of the edges
	for (unsigned short i=0; i<3; i++)
	{
		if ( (this->face(i))->contains_p(pt) )
		{
			return true;
		}
	}
	// projection onto 2-space to dispense with numerical errors when trying to calculate in dim()-space
	std::vector< std::valarray<T> > n(3, std::valarray<T>(2));;
	std::valarray<T> pt2(2);
	bool found = false;
	for (unsigned short i=0; i<this->dim()-1; i++)
	{
		for (unsigned short j=i; j<this->dim(); j++)
		{
			if (
				// we don't want a 1-space
				( i == j ) ||
				// we don't want a 2-space in which the triangle collapses into a segment, or point
				( ( eq( (*(this->p(0)))[i], (*(this->p(1)))[i] ) &&
				   eq( (*(this->p(1)))[i], (*(this->p(2)))[i] ) ) ||
				 ( eq( (*(this->p(0)))[j], (*(this->p(1)))[j] ) &&
				  eq( (*(this->p(1)))[j], (*(this->p(2)))[j] ) ) ||
				 ( eq( (*(this->p(0)))[i], (*(this->p(1)))[i] ) &&
				  eq( (*(this->p(0)))[j], (*(this->p(1)))[j] ) ) ||
				 ( eq( (*(this->p(1)))[i], (*(this->p(2)))[i] ) &&
				  eq( (*(this->p(1)))[j], (*(this->p(2)))[j] ) ) )
				)
			{
				continue;
			}
			n[0][0] = (*(this->p(0)))[i];
			n[0][1] = (*(this->p(0)))[j];
			n[1][0] = (*(this->p(1)))[i];
			n[1][1] = (*(this->p(1)))[j];
			n[2][0] = (*(this->p(2)))[i];
			n[2][1] = (*(this->p(2)))[j];
			pt2[0] = (*pt)[i];
			pt2[1] = (*pt)[j];
			found = true;
			break;
		}
		if ( found )
		{
			break;
		}
	}
	if ( !found )
	{
		return false;
	}
	
	// see if the point lies inside the triangle by seeing which side of each edge it lies on (nifty... huh!)
	for (unsigned short i=0; i<3; i++)
	{
		// ref is perpendicular to ( n_2 - n_1 ), ie. rotate by 90 degrees
		std::valarray<T> ref = n[(i+2)%3] - n[(i+1)%3];
		std::swap(ref[0], ref[1]);
		ref[0] *= T(-1);
		T dotref = dot(ref, std::valarray<T>(n[i] - n[(i+1)%3]));
		T dotrel = dot(ref, std::valarray<T>(pt2 - n[(i+1)%3]));
		if ( ( dotref >= T(0) && dotrel < T(0) ) ||
			( dotref <= T(0) && dotrel > T(0) ) )
		{
			return false;
		}
	}
	return true;
}





/* ---------- intersect ---------- */
/// in 2d parallel refers to whether or not l lies along any of the sides of the triangle
/// in 3d parallel refers to whether or not l lies in the plane of the triangle
template<class T> T triangle<T>::intersect(const line<T>& l, std::valarray<T>& intercept, bool& parallel) const
{
	intercept.resize(0);
	parallel = false;
	if ( this->dim() != l.dim() )
	{
		throw;
        return T(0);
	}
	
	// firstly, make sure the line actually intersects the plane containing the triangle
	if ( !this->container_.intersect(l, intercept, parallel) )
	{
		intercept.resize(0);
		parallel = false;
		return T(0);
	}
	std::vector< std::valarray<T> > itemp(0);
	if ( intercept.size() && !this->contains_p(&intercept, true) )
	{
		parallel = false;
		intercept.resize(0);
		return T(0);
	}
	else
	{
		for (unsigned short i=0; i<intercept.size()/this->dim(); i++)
		{
			itemp.push_back(intercept[std::slice(i*this->dim(), this->dim(), 1)]);
		}
	}
	// if the line is parallel to the triangle plane then we need to find the intersects along the edges
	if ( parallel )
	{
		if ( this->dim() < 3 )
		{
			parallel = false;
		}
		for (unsigned short i=0; i<3; i++)
		{
			bool on_edge = false;
			if ( !is_zero((this->face(i))->intersect(l, intercept, on_edge)) )
			{
				if ( on_edge && this->dim() < 3 )
				{
					parallel = true;
				}
				for (unsigned short j=0; j<intercept.size()/this->dim(); j++)
				{
					itemp.push_back(intercept[std::slice(j*this->dim(), this->dim(), 1)]);
				}
			}
		}
	}
	remove_duplicates(itemp);
	intercept.resize(itemp.size() * this->dim());
	for (unsigned short i=0; i<itemp.size(); i++)
	{
		intercept[std::slice(i * this->dim(), this->dim(), 1)] = itemp[i];
	}
	if ( !intercept.size() )
	{
		return T(0);
	}
	else if ( itemp.size() == 1 && this->dim() > 2 )
	{
		return T(1);
	}
	return maximise_distance(itemp);
}





/* ---------- face ---------- */
template<class T> inline const compact_shape<T>* triangle<T>::face(const unsigned short& en) const
{
	switch(en)
	{
		case 0:
			return &(this->edge0_);
			break;
		case 1:
			return &(this->edge1_);
			break;
		case 2:
			return &(this->edge2_);
			break;
	}
	return 0;
}





#endif /* _TRIANGLE_ */
