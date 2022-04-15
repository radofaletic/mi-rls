/*
 rotation
 
 A functional for rotating vectors
 
 Rado Faletic
 Department of Physics
 Faculty of Science
 Australian National University  ACT  0200
 Australia
 
 Rado.Faletic@anu.edu.au
 13th October 2004
 */


#ifndef _ROTATION_
#define _ROTATION_


/* ---------- header files ---------- */
#include <cmath>
#include <sstream>
#include <string>
#include <valarray>
/* -------- user header files ------- */
#include "angles.h"
#include "extra_math.h"
#include "matrix.h"
/* ---------------------------------- */


/* -------------------------------------------------- */
/* ---------- class & function declaration ---------- */
/* -------------------------------------------------- */


template<class T> class Rotation
{
private:
	Matrix<T> r_m_;
	bool set_r_m_;
	std::valarray<T> origin_;
	
	void set_dim(const unsigned short&);
	
public:
	
	Rotation(const unsigned short&);
	
	~Rotation();
	
	void reset();
	
	void clear();
	
	unsigned short dim() const;
	
	std::string print() const;
	
	std::string origin_print() const;
	
	void set_origin(const std::valarray<T>&);
	
	void clear_origin();
	void reset_origin();
	
	std::valarray<T> operator()(const std::valarray<T>&) const;
	std::valarray<T> O(const std::valarray<T>&) const;
	
	void set(const T&, const Angle::axes& = Angle::XY);
	void reset(const T&, const Angle::axes& = Angle::XY);
};


/* ------------------------------------------ */
/* ---------- function definitions ---------- */
/* ------------------------------------------ */


/* ---------- set_dim ---------- */
/* create a rotation matrix
 that is the dimension*
 dimesion identity matrix  */
template<class T> inline void
Rotation<T>::set_dim(const unsigned short& dimension)
{
	this->set_r_m_ = false;
	this->r_m_.init(dimension);
}
/* ------------------------------- */

/* ---------- Rotation ---------- */
/* create a rotation matrix of
 size dimension*dimension, ie
 the identity matrix            */
template<class T> inline
Rotation<T>::Rotation(const unsigned short& dimension)
{
	this->set_dim(dimension);
}
/* ------------------------------ */

/* ---------- ~Rotation ---------- */
template<class T> inline
Rotation<T>::~Rotation()
{
	this->clear();
}
/* ------------------------------- */

/* ---------- reset ---------- */
/* reset the rotation matrix
 to the identity matrix      */
template<class T> inline void
Rotation<T>::reset()
{
	this->set_r_m_ = false;
	this->r_m_.init(this->dim());
}
/* --------------------------- */

/* ---------- clear ---------- */
/* delete the rotation matrix
 and the origin              */
template<class T> inline void
Rotation<T>::clear()
{
	this->set_r_m_ = false;
	this->r_m_.clear();
	this->clear_origin();
}
/* --------------------------- */

/* ---------- dim ---------- */
/* return the dimension of
 the rotation matrix       */
template<class T> inline unsigned short
Rotation<T>::dim() const
{
	return this->r_m_.col_dim();
}
/* --------------------------- */

/* ---------- print ---------- */
/* print the rotation matrix   */
template<class T> std::string
Rotation<T>::print() const
{
	return this->r_m_.print();
}
/* --------------------------- */

/* ---------- origin_print ---------- */
/* print the origin of the rotation   */
template<class T> std::string
Rotation<T>::origin_print() const
{
	std::ostringstream tmp;
	tmp << '(';
	for (unsigned short i=0; i<(this->origin_).size(); i++)
	{
		if ( i )
		{
			tmp << ',';
		}
		tmp << this->origin_[i];
	}
	tmp << ')';
	return tmp.str();
}
/* ---------------------------------- */

/* ---------- set_origin ---------- */
/* set the origin of the rotation   */
template<class T> void
Rotation<T>::set_origin(const std::valarray<T>& new_origin)
{
	if ( is_zero(new_origin,T(1)) )
	{
		this->clear_origin();
		return;
	}
	this->origin_.resize(new_origin.size());
	this->origin_ = new_origin;
}
/* -------------------------------- */

/* ---------- clear_origin ---------- */
/* delete the rotation origin         */
template<class T> inline void
Rotation<T>::clear_origin()
{
	this->origin_.resize(0);
}
/* ---------------------------------- */

/* ---------- reset_origin ---------- */
/* reset the rotation origin to zero  */
template<class T> inline void
Rotation<T>::reset_origin()
{
	this->clear_origin();
}
/* ---------------------------------- */

/* ---------- operator() ---------- */
/* rotate the vector `vec' about the
 stored origin                    */
template<class T> std::valarray<T>
Rotation<T>::operator()(const std::valarray<T>& vec) const
{
	std::valarray<T> origvec = vec;
	
	// determine if the origin is to be shifted
	bool ors = ( this->origin_.size() );
	
	// shift the vector into the origin space
	if ( ors )
	{
		origvec -= this->origin_;
	}
	
	// there's little point doing calculations if they all return zero
	if ( is_zero(origvec,T(1)) )
	{
		return vec;
	}
	
	// do the rotation
	origvec = this->r_m_ * origvec;
	
	// shift the vector back into real space
	if ( ors )
	{
		origvec += this->origin_;
	}
	
	return origvec;
}
/* -------------------------------- */

/* ---------- O ---------- */
/* rotate the vector `vec'
 about the zero vector   */
template<class T> std::valarray<T>
Rotation<T>::O(const std::valarray<T>& vec) const
{
	return this->r_m_ * vec;
}
/* ----------------------- */

/* ---------- set ---------- */
/* assign the rotation matrix
 using the given angle and
 axis                      */
template<class T> void
Rotation<T>::set(const T& angle, const Angle::axes& input_axes)
{
	const unsigned short ldim = this->dim();
	Matrix<T> n_m(this->dim());
	unsigned short axis = 2;
	switch (input_axes)
	{
		case Angle::XY: case Angle::Z: // rotate the XY plane about the Z axis
			axis = 2;
			break;
		case Angle::YZ: case Angle::X: // rotate the YZ plane about the X axis
			axis = 0;
			break;
		case Angle::ZX: case Angle::Y: // rotate the ZX plane about the Y axis
			axis = 1;
			break;
	}
	switch(this->dim())
	{
		case 2: // 2-dimensional
			if ( is_zero(angle))
			{
				n_m(0,0) = 1;
				n_m(0,1) = 0;
				n_m(1,0) = 0;
				n_m(1,1) = 1;
			}
			else if ( eq(angle, T(90)) || eq(angle, T(-270)) )
			{
				n_m(0,0) = 0;
				n_m(0,1) = -1;
				n_m(1,0) = 1;
				n_m(1,1) = 0;
			}
			else if ( eq(angle, T(180)) || eq(angle, T(-180)) )
			{
				n_m(0,0) = -1;
				n_m(0,1) = 0;
				n_m(1,0) = 0;
				n_m(1,1) = -1;
			}
			else if ( eq(angle, T(270)) || eq(angle, T(-90)) )
			{
				n_m(0,0) = 0;
				n_m(0,1) = 1;
				n_m(1,0) = -1;
				n_m(1,1) = 0;
			}
			else
			{
				n_m(0,0) = cos(Angle::d2r(angle));
				n_m(0,1) = -sin(Angle::d2r(angle));
				n_m(1,0) = sin(Angle::d2r(angle));
				n_m(1,1) = cos(Angle::d2r(angle));
			}
			break;
		default: // 3-dimensional
			// use modulo arithmetic to determine the indexing
			// for the axis of rotation
			if ( is_zero(angle) )
			{
				n_m((axis+1)%ldim,(axis+1)%ldim) = 1;
				n_m((axis+1)%ldim,(axis+2)%ldim) = 0;
				n_m((axis+2)%ldim,(axis+1)%ldim) = 0;
				n_m((axis+2)%ldim,(axis+2)%ldim) = 1;
			}
			else if ( eq(angle, T(90)) || eq(angle, T(-270)) )
			{
				n_m((axis+1)%ldim,(axis+1)%ldim) = 0;
				n_m((axis+1)%ldim,(axis+2)%ldim) = -1;
				n_m((axis+2)%ldim,(axis+1)%ldim) = 1;
				n_m((axis+2)%ldim,(axis+2)%ldim) = 0;
			}
			else if ( eq(angle, T(180)) || eq(angle, T(-180)) )
			{
				n_m((axis+1)%ldim,(axis+1)%ldim) = -1;
				n_m((axis+1)%ldim,(axis+2)%ldim) = 0;
				n_m((axis+2)%ldim,(axis+1)%ldim) = 0;
				n_m((axis+2)%ldim,(axis+2)%ldim) = -1;
			}
			else if ( eq(angle, T(270)) || eq(angle, T(-90)) )
			{
				n_m((axis+1)%ldim,(axis+1)%ldim) = 0;
				n_m((axis+1)%ldim,(axis+2)%ldim) = 1;
				n_m((axis+2)%ldim,(axis+1)%ldim) = -1;
				n_m((axis+2)%ldim,(axis+2)%ldim) = 0;
			}
			else
			{
				n_m((axis+1)%ldim,(axis+1)%ldim) = cos(Angle::d2r(angle));
				n_m((axis+1)%ldim,(axis+2)%ldim) = -sin(Angle::d2r(angle));
				n_m((axis+2)%ldim,(axis+1)%ldim) = sin(Angle::d2r(angle));
				n_m((axis+2)%ldim,(axis+2)%ldim) = cos(Angle::d2r(angle));
			}
			break;
	}
	if ( this->set_r_m_ )
	{
		this->r_m_ = n_m * this->r_m_;
	}
	else
	{
		this->r_m_ = n_m;
		this->set_r_m_ = true;
	}
}
/* ------------------------- */

/* ---------- reset ---------- */
/* reset the rotation matrix,
 then reassign values to it  */
template<class T> inline void
Rotation<T>::reset(const T& angle, const Angle::axes& axis)
{
	this->reset();
	this->set(angle, axis);
}
/* --------------------------- */


#endif /* _ROTATION_ */
