/**
 line
 
 A collection of C++ routines for dealing with straight lines.
 
 Rado Faletic
 Department of Physics
 Faculty of Science
 Australian National University  ACT  0200
 Australia
 
 Rado.Faletic@anu.edu.au
 17th May 2004
 4th May 2022 
 */





#ifndef _LINE_
#define _LINE_





/* ---------- standard header files ---------- */
#include <sstream>
#include <string>
#include <valarray>





/* ------- user header files -------- */
#include "extra_math.h"
#include "shape.h"





/* ---------- class & function declaration ---------- */

// 'line' class
template<class T> class line : public shape<T>
{
private:
	std::valarray<T> zero_;
	std::valarray<T> gradient_;
public:
	void init(const std::valarray<T>&, const std::valarray<T>&);
	void init(const std::valarray<T>*, const std::valarray<T>*);
	line(const std::valarray<T>&, const std::valarray<T>&);
	line(const std::valarray<T>*, const std::valarray<T>*);
	line() { };
    ~line() { };
	unsigned short dim() const;
	std::string print() const;
	const std::valarray<T>* gradient_p() const;
	std::valarray<T> zero() const;
	std::valarray<T> evaluate(const T& = T(0)) const;
	T distance_p(const std::valarray<T>*) const;
	std::valarray<T> closest(const std::valarray<T>&) const;
	bool contains_p(const std::valarray<T>*) const;
	T intersect(const line<T>&, std::valarray<T>&, bool&) const;
};





/* ---------- function definitions ---------- */





/* ---------- init ---------- */
template<class T> void line<T>::init(const std::valarray<T>& gr, const std::valarray<T>& ze)
{
	if ( gr.size() != ze.size() || !gr.size() )
	{
		throw;
        return;
	}
	
	// allocate line properties
	this->gradient_.resize(gr.size());
	this->gradient_ = gr;
	this->zero_.resize(ze.size());
	this->zero_ = ze;
}





/* ---------- init ---------- */
template<class T> void line<T>::init(const std::valarray<T>* pt1, const std::valarray<T>* pt2)
{
	this->init(*pt2 - *pt1, *pt1);
}





/* ---------- line ---------- */
template<class T> inline line<T>::line(const std::valarray<T>& gr, const std::valarray<T>& ze)
{
	this->init(gr, ze);
}





/* ---------- line ---------- */
template<class T> inline line<T>::line(const std::valarray<T>* pt1, const std::valarray<T>* pt2)
{
	this->init(pt1, pt2);
}





/* ---------- dim ---------- */
template<class T> inline unsigned short line<T>::dim() const
{
	return this->zero_.size();
}





/* ---------- print ---------- */
template<class T> std::string line<T>::print() const
{
	std::ostringstream tmp;
	tmp << '(';
	for (unsigned short i=0; i<(this->gradient_).size(); i++)
	{
		if ( i )
		{
			tmp << ',';
		}
		tmp << (this->gradient_)[i];
	}
	tmp << ") * t + (";
	for (unsigned short i=0; i<(this->zero_).size(); i++)
	{
		if ( i )
		{
			tmp << ',';
		}
		tmp << (this->zero_)[i];
	}
	tmp << ')';
	return tmp.str();
}





/* ---------- gradient_p ---------- */
template<class T> inline const std::valarray<T>* line<T>::gradient_p() const
{
	return &(this->gradient_);
}





/* ---------- zero ---------- */
template<class T> inline std::valarray<T> line<T>::zero() const
{
	return this->zero_;
}





/* ---------- evaluate ---------- */
template<class T> inline std::valarray<T> line<T>::evaluate(const T& t_value) const
{
	return this->zero_ + ( this->gradient_ * t_value );
}





/* ---------- distance_p ---------- */
template<class T> inline T line<T>::distance_p(const std::valarray<T>* pt) const
{
	// use projected length along the line, then Pythagoras
	std::valarray<T> tmp = *pt - this->zero_;
	T d = std::abs( norm2(tmp) - std::abs( dot(tmp, this->gradient_) ) );
	if ( eq(d,T(1)) )
	{
		d = T(1);
	}
	return ( is_zero(d) ) ? T(0) : std::sqrt(d);
	//return norm(std::valarray<T>((*pt) - this->closest(*pt))); // OLD METHOD
}





/* ---------- closest ---------- */
template<class T> inline std::valarray<T> line<T>::closest(const std::valarray<T>& pt) const
{
	return this->evaluate(dot(std::valarray<T>(pt - this->zero_), this->gradient_));
}





/* ---------- contains_p ---------- */
template<class T> inline bool line<T>::contains_p(const std::valarray<T>* pt) const
{
	return is_zero(this->distance_p(pt));
}





/* ---------- intersect ---------- */
template<class T> T line<T>::intersect(const line<T>& line2, std::valarray<T>& intercept, bool& parallel) const
{
	intercept.resize(0);
	parallel = false;
	if ( this->dim() != line2.dim() )
	{
		throw;
        return T(0);
	}
	else if ( this->dim() == 1 )
	{
		parallel = true;
		return T(1);
	}
	else if ( v_pm_eq(&(this->gradient_), line2.gradient_p()) ) // parallel lines do not intersect
	{
		if ( this->contains(line2.zero()) )
		{
			parallel = true;
			return T(1);
		}
		return T(0);
	}
	std::valarray<T> zero2 = line2.zero();
	std::valarray<T> gradient2 = *(line2.gradient_p());
	if ( this->dim() == 3 ) // use Example 14, Swokowksi pg 727 to work out if the lines intersect
	{
		std::valarray<T> cp = cross(this->gradient_, gradient2);
		T d = std::abs( dot(cp, std::valarray<T>(zero2 - this->zero_)) ) / norm(cp);
		if ( !is_zero(d) )
		{
			return T(0);
		}
	}
	
	// solve for the free variables t_1 and t_2 by projection onto 2-space
	// NOTE: these do not necessarily solve for the interestion, but will be correct in the event that the lines actually do intersect in <dim>-space, otherwise intersect return falsly (above). Also, we take the average of the two solutions, to try and smear out numerical errors.
	T t_1 = T(0);
	T t_2 = T(0);
	for (unsigned short i=0; i<this->dim()-1; i++)
	{
		for (unsigned short j=i+1; j<this->dim(); j++)
		{
			T det = this->gradient_[i] * (-gradient2[j]) - (-gradient2[i]) * this->gradient_[j];
			if ( is_zero(det) )
			{
				continue;
			}
			t_1 = ( (-gradient2[j]) * ( zero2[i] - this->zero_[i] ) + gradient2[i] * ( zero2[j] - this->zero_[j] ) ) / det;
			t_2 = ( (-this->gradient_[j]) * ( zero2[i] - this->zero_[i] ) + this->gradient_[i] * ( zero2[j] - this->zero_[j] ) ) / det;
			break;
		}
	}
	
	// now create the actual intercept
	intercept.resize(this->dim());
	intercept = this->evaluate(t_1);
	if ( this->dim() == 3 )
	{
		intercept = line2.evaluate(t_2); // t_2 doesn't always seem to work in 2d... hmmm... funny numerics
		intercept /= T(2);
	}
	
	return T(1);
}





#endif /* _LINE_ */
