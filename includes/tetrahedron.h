/**
 tetrahedron
 
 A collection of C++ routines for dealing with generic tetrahedrons.
 
 Rado Faletic
 Department of Physics
 Faculty of Science
 Australian National University  ACT  0200
 Australia
 
 Rado.Faletic@anu.edu.au
 21st June 2004
 22nd April 2022
 */





#ifndef _TETRAHEDRON_
#define _TETRAHEDRON_





/* ---------- standard header files ---------- */
#include <algorithm>
#include <valarray>
#include <vector>





/* ------- user header files -------- */
#include "extra_math.h"
#include "line.h"
#include "plane.h"
#include "shape.h"
#include "triangle.h"





/* ---------- class declaration ---------- */

// `tetrahedron' class
template<class T> class tetrahedron : public compact_shape<T>
{
private:
	triangle<T> face0_;
	triangle<T> face1_;
	triangle<T> face2_;
	triangle<T> face3_;
public:
	void init(const std::valarray<T>*, const std::valarray<T>*, const std::valarray<T>*, const std::valarray<T>*);
	tetrahedron(const std::valarray<T>*, const std::valarray<T>*, const std::valarray<T>*, const std::valarray<T>*);
	tetrahedron() { };
	T make_measure() const;
	T scale() const;
	T distance_p(const std::valarray<T>*) const;
	bool contains_p(const std::valarray<T>*) const;
	T intersect(const line<T>&, std::valarray<T>&, bool&) const;
	const compact_shape<T>* face(const unsigned short&) const;
};





/* ---------- function definitions ---------- */





/* ---------- init ---------- */
template<class T> void tetrahedron<T>::init(const std::valarray<T>* pt1, const std::valarray<T>* pt2, const std::valarray<T>* pt3, const std::valarray<T>* pt4)
{
	if ( pt1->size() < 3 )
	{
		throw;
        return;
	}
	
	// just make sure that the vertices are not coplanar
	plane<T> faceplane(pt1, pt2, pt3);
	if ( faceplane.contains_p(pt4) )
	{
		throw;
        return;
	}
	
	// allocate tetrahedron nodes
	this->assign_nodes(4, pt1, pt2, pt3, pt4);
	
	this->face0_.init(this->p(1), this->p(2), this->p(3));
	this->face1_.init(this->p(2), this->p(3), this->p(0));
	this->face2_.init(this->p(3), this->p(0), this->p(1));
	this->face3_.init(this->p(0), this->p(1), this->p(2));
	
	// volume
	//this->measure_ = this->make_measure();
}





/* ---------- tetrahedron ---------- */
template<class T> inline tetrahedron<T>::tetrahedron(const std::valarray<T>* pt1, const std::valarray<T>* pt2, const std::valarray<T>* pt3, const std::valarray<T>* pt4)
{
	this->init(pt1, pt2, pt3, pt4);
}





/* ---------- make_measure ---------- */
template<class T> inline T tetrahedron<T>::make_measure() const
{
	T m = T(0);
	// take the determinant of the three
	// vectors defining the tetrahedron
	m += ( ((*this)[1])[0] - ((*this)[0])[0] )
	* ( ((*this)[2])[1] - ((*this)[0])[1] )
	* ( ((*this)[3])[2] - ((*this)[0])[2] );
	m -= ( ((*this)[1])[2] - ((*this)[0])[2] )
	* ( ((*this)[2])[1] - ((*this)[0])[1] )
	* ( ((*this)[3])[0] - ((*this)[0])[0] );
	m += ( ((*this)[1])[1] - ((*this)[0])[1] )
	* ( ((*this)[2])[2] - ((*this)[0])[2] )
	* ( ((*this)[3])[0] - ((*this)[0])[0] );
	m -= ( ((*this)[1])[0] - ((*this)[0])[0] )
	* ( ((*this)[2])[2] - ((*this)[0])[2] )
	* ( ((*this)[3])[1] - ((*this)[0])[1] );
	m += ( ((*this)[1])[2] - ((*this)[0])[2] )
	* ( ((*this)[2])[0] - ((*this)[0])[0] )
	* ( ((*this)[3])[1] - ((*this)[0])[1] );
	m -= ( ((*this)[1])[1] - ((*this)[0])[1] )
	* ( ((*this)[2])[0] - ((*this)[0])[0] )
	* ( ((*this)[3])[2] - ((*this)[0])[2] );
	return m / T(6); // and divide by 6
}





/* ---------- scale ---------- */
template<class T> inline T tetrahedron<T>::scale() const
{
	return ( this->face0_.scale() + this->face1_.scale() + this->face2_.scale() + this->face3_.scale() ) / T(4);
}





/* ---------- distance_p ---------- */
template<class T> inline T tetrahedron<T>::distance_p(const std::valarray<T>* pt) const
{
	if ( this->contains_p(pt) )
	{
		return T(0);
	}
	T dist = this->face0_.distance_p(pt);
	dist = std::min(dist, this->face1_.distance_p(pt));
	dist = std::min(dist, this->face2_.distance_p(pt));
	dist = std::min(dist, this->face3_.distance_p(pt));
	return dist;
}





/* ---------- contains_p ---------- */
template<class T> bool tetrahedron<T>::contains_p(const std::valarray<T>* pt) const
{
	if ( this->dim() != pt->size() )
	{
		throw;
        return false;
	}
	
	// see if the point lies inside the tetrahedron
	for (unsigned short i=0; i<4; i++)
	{
		const triangle<T>* tf;
		switch(i)
		{
			case 0:
				tf = &(this->face0_);
				break;
			case 1:
				tf = &(this->face1_);
				break;
			case 2:
				tf = &(this->face2_);
				break;
			case 3:
				tf = &(this->face3_);
				break;
		}
		if ( tf->contains_p(pt) )
		{
			return true;
		}
		T dotref = dot(*(tf->normal_p()), std::valarray<T>(*(this->p(i)) - *(this->p((i+1)%4))));
		T dotrel = dot(*(tf->normal_p()), std::valarray<T>((*pt) - *(this->p((i+1)%4))));
		if ( ( dotref >= T(0) && dotrel < T(0) ) ||
			( dotref <= T(0) && dotrel > T(0) ) )
		{
			return false;
		}
	}
	return true;
}





/* ---------- intersect ---------- */
/// 'on_edge' refers to whether or not l lies along any of the triangular faces of the tetrahedron
template<class T> T tetrahedron<T>::intersect(const line<T>& l, std::valarray<T>& intercepts, bool& on_edge) const
{
	intercepts.resize(0);
	on_edge = false;
	if ( this->dim() != l.dim() )
	{
		throw;
        return T(0);
	}
	
	std::vector< std::valarray<T> > itemp;
	
	// test each face
	bool parallel = false;
	for (unsigned short i=0; i<4; i++)
	{
		if ( !is_zero((this->face(i))->intersect(l, intercepts, parallel)) )
		{
			if ( parallel )
			{
				on_edge = true;
			}
			for (unsigned short j=0; j<intercepts.size()/this->dim(); j++)
			{
				itemp.push_back(intercepts[std::slice(j * this->dim(), this->dim(), 1)]);
			}
		}
	}
	
	remove_duplicates(itemp);
	intercepts.resize(itemp.size()*this->dim());
	for (unsigned short i=0; i<itemp.size(); i++)
	{
		intercepts[std::slice(i * this->dim(), this->dim(), 1)] = itemp[i];
	}
	return maximise_distance(itemp);
}





/* ---------- face ---------- */
template<class T> inline const compact_shape<T>* tetrahedron<T>::face(const unsigned short& fn) const
{
	switch(fn)
	{
		case 0:
			return &(this->face0_);
			break;
		case 1:
			return &(this->face1_);
			break;
		case 2:
			return &(this->face2_);
			break;
		case 3:
			return &(this->face3_);
			break;
	}
	return 0;
}





#endif /* _TETRAHEDRON_ */
