/*
  polygon
  
  A collection of C++ routines for dealing with generic polygons.
  
  Rado Faletic
  Department of Physics
  Faculty of Science
  Australian National University  ACT  0200
  Australia
  
  Rado.Faletic@anu.edu.au
  21st June 2004
*/


#ifndef _POLYGON_
#define _POLYGON_


/* ---------- standard header files ---------- */
#include <algorithm>
#include <valarray>
#include <vector>
/* ------- user header files -------- */
#include "extra_math.h"
#include "line.h"
#include "segment.h"
#include "shape.h"
#include "tetrahedron.h"
#include "triangle.h"
/* ---------------------------------- */


/* --------------------------------------- */
/* ---------- class declaration ---------- */
/* --------------------------------------- */


// `polygon' class
template<class T>
class polygon : public compact_shape<T>
{
 private:
  std::valarray<T> centroid_;
  std::vector< triangle<T> > panel_;
  void init_();
 public:
  void init(const std::valarray<T>*, const std::valarray<T>*, const std::valarray<T>*);
  polygon(const std::valarray<T>*, const std::valarray<T>*, const std::valarray<T>*);
  void init(const std::valarray<T>*, const std::valarray<T>*, const std::valarray<T>*, const std::valarray<T>*);
  polygon(const std::valarray<T>*, const std::valarray<T>*, const std::valarray<T>*, const std::valarray<T>*);
  polygon() { };
  T make_measure() const;
  T scale() const;
  T distance_p(const std::valarray<T>*) const;
  bool contains_p(const std::valarray<T>*) const;
  T intersect(const line<T>&, std::valarray<T>&, bool&) const;
  const compact_shape<T>* face(const unsigned short&) const;
  const std::valarray<T>* centroid_p() const;
};


/* ------------------------------------------ */
/* ---------- function definitions ---------- */
/* ------------------------------------------ */


/* ---------- init_ ---------- */
template<class T> void
polygon<T>::init_()
{
  // check for consistency
  if ( this->dim() < 2 )
    {
      throw; return;
    }

  // panels
  switch(this->size())
    {
    case 3: // triangle
      this->panel_.resize(1);
      this->panel_[0].init(this->p(0), this->p(1), this->p(2));
      break;
    default: // generic polygon
      this->centroid_.resize(this->dim(), T(0));
      for (unsigned short i=0; i<this->size(); i++)
	{
	  this->centroid_ += *(this->p(i));
	}
      this->centroid_ /= T(this->size());
      this->panel_.resize(this->size());
      for (unsigned short i=0; i<this->size(); i++)
	{
	  this->panel_[i].init(&(this->centroid_), this->p(i), this->p((i+1)%this->size()));
	}
      break;
    }

  // area
  //this->measure_ = this->make_measure();
}
/* --------------------------- */

/* ---------- init ---------- */
template<class T> void
polygon<T>::init(const std::valarray<T>* pt1, const std::valarray<T>* pt2, const std::valarray<T>* pt3)
  // this version of 'init' is a triangle
{
  this->assign_nodes(3, pt1, pt2, pt3);
  this->init_();
}
/* -------------------------- */

/* ---------- polygon ---------- */
template<class T> inline
polygon<T>::polygon(const std::valarray<T>* pt1, const std::valarray<T>* pt2, const std::valarray<T>* pt3)
{
  this->init(pt1, pt2, pt3);
}
/* -------------------------- */

/* ---------- init ---------- */
template<class T> void
polygon<T>::init(const std::valarray<T>* pt1, const std::valarray<T>* pt2, const std::valarray<T>* pt3, const std::valarray<T>* pt4)
  // this version of 'init' is a quadrilateral, with nodes in cyclic order
{
  this->assign_nodes(4, pt1, pt2, pt3, pt4);
  this->init_();
}
/* -------------------------- */

/* ---------- polygon ---------- */
template<class T> inline
polygon<T>::polygon(const std::valarray<T>* pt1, const std::valarray<T>* pt2, const std::valarray<T>* pt3, const std::valarray<T>* pt4)
{
  this->init(pt1, pt2, pt3, pt4);
}
/* -------------------------- */

/* ---------- make_measure ---------- */
template<class T> inline T
polygon<T>::make_measure() const
{
  T m = T(0);
  for (unsigned short i=0; i<this->panel_.size(); i++)
    {
      m += this->panel_[i].measure();
    }
  return m;
}
/* ---------------------------------- */

/* ---------- scale ---------- */
template<class T> inline T
polygon<T>::scale() const
{
  switch (this->size())
    {
    case 3:
      return (this->panel_[0]).scale();
      break;
    }
  T s = T(0);
  for (unsigned short i=0; i<this->panel_.size(); i++)
    {
      s += (this->panel_[i].face(0))->scale();
    }
  return s / T(this->panel_.size());
}
/* --------------------------- */

/* ---------- distance_p ---------- */
template<class T> inline T
polygon<T>::distance_p(const std::valarray<T>* pt) const
{
  T dist = this->panel_[0].distance_p(pt);
  for (unsigned short i=1; i<this->panel_.size(); i++)
    {
      dist = std::min(dist, this->panel_[i].distance_p(pt));
    }
  return dist;
}
/* -------------------------------- */

/* ---------- contains_p ---------- */
template<class T> inline bool
polygon<T>::contains_p(const std::valarray<T>* pt) const
{
  return is_zero(this->distance_p(pt));
}
/* -------------------------------- */

/* ---------- intersect ---------- */
template<class T> T
polygon<T>::intersect(const line<T>& l, std::valarray<T>& intercept, bool& parallel) const
  // in 2d 'parallel' refers to whether l lies along any of the edges of the polygon
  // in 3d 'parallel' refers to whether l lies in all of the pieces of the polygon
{
  intercept.resize(0);
  parallel = false;
  if ( this->dim() != l.dim() )
    {
      throw; return T(0);
    }
  std::vector< std::valarray<T> > itemp(0);

  switch (this->dim())
    {
    case 2: // in 2d we only need to check the edges (much quicker)
      for (unsigned short i=0; i<this->size(); i++)
	{
	  bool this_edge = false;
	  if ( !is_zero((this->face(i))->intersect(l, intercept, this_edge)) )
	    {
	      if ( this_edge )
		{
		  parallel = true;
		}
	      for (size_t j=0; j<intercept.size()/this->dim(); j++)
		{
		  itemp.push_back(intercept[std::slice(j * this->dim(), this->dim(), 1)]);
		}
	    }
	}
      break;
    default: // in 3d we check each "piece" of the polygon
      unsigned short oec = 0;
      for (unsigned short i=0; i<(this->panel_).size(); i++)
	{
	  bool this_edge = false;
	  if ( !is_zero(this->panel_[i].intersect(l, intercept, this_edge)) )
	    {
	      if ( this_edge )
		{
		  oec++;
		}
	      for (unsigned short j=0; j<intercept.size()/this->dim(); j++)
		{
		  itemp.push_back(intercept[std::slice(j * this->dim(), this->dim(), 1)]);
		}
	    }
	}
      if ( oec == this->size() )
	{
	  parallel = true;
	}
      break;
    }
  
  remove_duplicates(itemp);
  intercept.resize(itemp.size()*this->dim());
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
/* ------------------------------- */

/* ---------- face ---------- */
template<class T> inline const compact_shape<T>*
polygon<T>::face(const unsigned short& en) const
{
  switch(this->size())
    {
    case 3: // triangle edge numbering is different to polygon edge numbering
      return (this->panel_[0]).face((en+2)%3);
      break;
    }
  return (this->panel_[en]).face(0);
}
/* -------------------------- */

/* ---------- centroid_p ---------- */
template<class T> inline const std::valarray<T>*
polygon<T>::centroid_p() const
{
  switch (this->dim())
    {
    case 2: // we don't use centroid_ in 2d
      return 0;
      break;
    }
  return &(this->centroid_);
}
/* -------------------------------- */


#endif /* _POLYGON_ */
