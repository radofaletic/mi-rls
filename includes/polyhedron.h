/*
 polyhedron
 
 A collection of C++ routines for dealing with generic polyhedra.
 
 Rado Faletic
 Department of Physics
 Faculty of Science
 Australian National University  ACT  0200
 Australia
 
 Rado.Faletic@anu.edu.au
 21st June 2004
 */


#ifndef _POLYHEDRON_
#define _POLYHEDRON_


/* ---------- standard header files ---------- */
#include <algorithm>
#include <valarray>
#include <vector>
/* ------- user header files -------- */
#include "extra_math.h"
#include "line.h"
#include "polygon.h"
#include "shape.h"
#include "tetrahedron.h"
/* ---------------------------------- */


/* --------------------------------------- */
/* ---------- class declaration ---------- */
/* --------------------------------------- */


// `polyhedron' class
template<class T>
class polyhedron : public compact_shape<T>
{
private:
	std::vector< polygon<T> >  face_;
	void init_();
public:
	void init(const std::valarray<T>*, const std::valarray<T>*, const std::valarray<T>*, const std::valarray<T>*);
	polyhedron(const std::valarray<T>*, const std::valarray<T>*, const std::valarray<T>*, const std::valarray<T>*);
	void init(const std::valarray<T>*, const std::valarray<T>*, const std::valarray<T>*, const std::valarray<T>*, const std::valarray<T>*, const std::valarray<T>*, const std::valarray<T>*, const std::valarray<T>*);
	polyhedron(const std::valarray<T>*, const std::valarray<T>*, const std::valarray<T>*, const std::valarray<T>*, const std::valarray<T>*, const std::valarray<T>*, const std::valarray<T>*, const std::valarray<T>*);
	polyhedron() { };
	T make_measure() const;
	T scale() const;
	T distance_p(const std::valarray<T>*) const;
	bool contains_p(const std::valarray<T>*) const;
	T intersect(const line<T>&, std::valarray<T>&, bool&) const;
	const compact_shape<T>* face(const unsigned short&) const;
};


/* ------------------------------------------ */
/* ---------- function definitions ---------- */
/* ------------------------------------------ */


/* ---------- init_ ---------- */
template<class T> void
polyhedron<T>::init_()
{
	// check for consistency
	if ( (this->p(0))->size() < 3 )
	{
		throw; return;
	}
	
	// volume
	//this->measure_ = this->make_measure();
}
/* --------------------------- */

/* ---------- init ---------- */
template<class T> void
polyhedron<T>::init(const std::valarray<T>* pt1, const std::valarray<T>* pt2, const std::valarray<T>* pt3, const std::valarray<T>* pt4)
// this version of 'init' is a tetrahedron
{
	// allocate polyhedron nodes
	this->assign_nodes(4, pt1, pt2, pt3, pt4);
	
	// create faces of the polyhedron
	this->face_.resize(4);
	this->face_[0].init(this->p(1),this->p(2),this->p(3));
	this->face_[1].init(this->p(2),this->p(3),this->p(0));
	this->face_[2].init(this->p(3),this->p(0),this->p(1));
	this->face_[3].init(this->p(0),this->p(1),this->p(2));
	
	this->init_();
}
/* -------------------------- */

/* ---------- polyhedron ---------- */
template<class T> inline
polyhedron<T>::polyhedron(const std::valarray<T>* pt1, const std::valarray<T>* pt2, const std::valarray<T>* pt3, const std::valarray<T>* pt4)
{
	this->init(pt1, pt2, pt3, pt4);
}
/* -------------------------------- */

/* ---------- init ---------- */
template<class T> void
polyhedron<T>::init(const std::valarray<T>* pt1, const std::valarray<T>* pt2, const std::valarray<T>* pt3, const std::valarray<T>* pt4, const std::valarray<T>* pt5, const std::valarray<T>* pt6, const std::valarray<T>* pt7, const std::valarray<T>* pt8)
// this version of 'init' assumes that the nodes are given in order as for a Plot3D grid cell
{
	// allocate polyhedron nodes
	this->assign_nodes(8, pt1, pt2, pt3, pt4, pt5, pt6, pt7, pt8);
	
	// create faces of the polyhedron
	this->face_.resize(6);
	this->face_[0].init(this->p(0),this->p(1),this->p(3),this->p(2)); // front
	this->face_[1].init(this->p(1),this->p(5),this->p(7),this->p(3)); // right
	this->face_[2].init(this->p(3),this->p(7),this->p(6),this->p(2)); // top
	this->face_[3].init(this->p(7),this->p(5),this->p(4),this->p(6)); // back
	this->face_[4].init(this->p(6),this->p(4),this->p(0),this->p(2)); // left
	this->face_[5].init(this->p(4),this->p(5),this->p(1),this->p(0)); // bottom
	
	this->init_();
}
/* -------------------------- */

/* ---------- polyhedron ---------- */
template<class T> inline
polyhedron<T>::polyhedron(const std::valarray<T>* pt1, const std::valarray<T>* pt2, const std::valarray<T>* pt3, const std::valarray<T>* pt4, const std::valarray<T>* pt5, const std::valarray<T>* pt6, const std::valarray<T>* pt7, const std::valarray<T>* pt8)
{
	this->init(pt1, pt2, pt3, pt4, pt5, pt6, pt7, pt8);
}
/* -------------------------------- */

/* ---------- make_measure ---------- */
template<class T> inline T
polyhedron<T>::make_measure() const
{
	T m = T(0);
	tetrahedron<T> piece;
	switch(this->size())
	{
		case 4: // tetrahedron
			piece.init(this->p(0), this->p(1), this->p(2), this->p(3));
			m = piece.measure();
			break;
		default:
			// calculate the centroid
			std::valarray<T> centroid(T(0), this->dim());
			for (unsigned short i=0; i<this->size(); i++)
			{
				centroid += *(this->p(i));
			}
			centroid /= T(this->size());
			// calculate tetrahedral volumetric pieces
			for (unsigned short i=0; i<this->face_.size(); i++)
			{
				switch((this->face_[i]).size())
				{
					case 3:
						piece.init(&centroid, (this->face_[i]).p(0), (this->face_[i]).p(1), (this->face_[i]).p(2));
						m += piece.measure();
						break;
					default:
						for (unsigned short j=0; j<(this->face_[i]).size(); j++)
						{
							piece.init(&centroid, (this->face_[i]).centroid_p(),
									   (this->face_[i]).p(j), (this->face_[i]).p((j+1)%(this->face_[i]).size()));
							m += piece.measure();
						}
						break;
				}
			}
			break;
	}
	return m;
}
/* ---------------------------------- */

/* ---------- scale ---------- */
template<class T> inline T
polyhedron<T>::scale() const
{
	T s = T(0);
	for (unsigned short i=0; i<this->face_.size(); i++)
	{
		s += face_[i].scale();
	}
	return s / T(this->face_.size());
}
/* --------------------------- */

/* ---------- distance_p ---------- */
template<class T> inline T
polyhedron<T>::distance_p(const std::valarray<T>* pt) const
{
	tetrahedron<T> piece;
	T dist;
	switch(this->size())
	{
		case 4: // tetrahedron
			piece.init(this->p(0), this->p(1), this->p(2), this->p(3));
			return piece.distance_p(pt);
			break;
		default:
			T dist = (this->face_[0]).distance_p(pt);
			// calculate the centroid
			std::valarray<T> centroid(T(0), this->dim());
			for (unsigned short i=0; i<this->size(); i++)
			{
				centroid += *(this->p(i));
			}
			centroid /= T(this->size());
			// calculate tetrahedral volumetric pieces
			for (unsigned short i=0; i<this->face_.size(); i++)
			{
				switch((this->face_[i]).size())
				{
					case 3:
						piece.init(&centroid, (this->face_[i]).p(0), (this->face_[i]).p(1), (this->face_[i]).p(2));
						dist = std::min(dist, piece.distance_p(pt));
						break;
					default:
						for (unsigned short j=0; j<(this->face_[i]).size(); j++)
						{
							piece.init(&centroid, (this->face_[i]).centroid_p(),
									   (this->face_[i]).p(j), (this->face_[i]).p((j+1)%(this->face_[i]).size()));
							dist = std::min(dist, piece.distance_p(pt));
						}
						break;
				}
			}
			break;
	}
	return dist;
}
/* -------------------------------- */

/* ---------- contains_p ---------- */
template<class T> bool
polyhedron<T>::contains_p(const std::valarray<T>* pt) const
{
	return ( is_zero(this->distance_p(pt)) ) ? true : false;
}
/* -------------------------------- */

/* ---------- intersect ---------- */
template<class T> T
polyhedron<T>::intersect(const line<T>& l, std::valarray<T>& intercepts, bool& on_edge) const
// on_edge refers to whether l lies along any of the faces of the polyhedron
{
	intercepts.resize(0);
	on_edge = false;
	if ( this->dim() != l.dim() )
	{
		throw; return T(0);
	}
	
	std::vector< std::valarray<T> > itemp;
	
	for (unsigned short i=0; i<this->face_.size(); i++)
	{
		bool parallel = false;
		if ( !is_zero(this->face_[i].intersect(l, intercepts, parallel)) )
		{
			if ( parallel )
			{
				on_edge = true;
			}
			for (size_t j=0; j<intercepts.size()/this->dim(); j++)
			{
				itemp.push_back(intercepts[std::slice(j * this->dim(), this->dim(), 1)]);
			}
		}
	}
	
	remove_duplicates(itemp);
	intercepts.resize(itemp.size()*this->dim());
	for (unsigned short i=0; i<itemp.size(); i++)
	{
		intercepts[std::slice(i*this->dim(), this->dim(), 1)] = itemp[i];
	}
	return maximise_distance(itemp);
}
/* ------------------------------- */

/* ---------- face ---------- */
template<class T> inline const compact_shape<T>*
polyhedron<T>::face(const unsigned short& fn) const
{
	return &(this->face_[fn]);
}
/* -------------------------------- */


#endif /* _POLYHEDRON_ */
