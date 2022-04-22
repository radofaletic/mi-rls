/**
 shape
 
 A collection of C++ routines for dealing with generic shapes.
 
 Rado Faletic
 Department of Physics
 Faculty of Science
 Australian National University  ACT  0200
 Australia
 
 Rado.Faletic@anu.edu.au
 21st June 2004
 22nd April 2022
 */





#ifndef _SHAPE_
#define _SHAPE_





/* ---------- standard header files ---------- */
#include <cstdarg>
#include <sstream>
#include <string>
#include <valarray>





/* ---------- class declaration ---------- */
template<class T> class line;

// `shape' base class
template<class T> class shape
{
public:
	shape() { };
	~shape() { };
	virtual unsigned short dim() const = 0;
	virtual std::string print() const = 0;
	virtual T distance_p(const std::valarray<T>*) const = 0;
	T distance(const std::valarray<T>& p) const { return this->distance_p(&p); };
	virtual bool contains_p(const std::valarray<T>*) const = 0;
	bool contains(const std::valarray<T>& p) const { return this->contains_p(&p); };
	virtual T intersect(const line<T>&, std::valarray<T>&, bool&) const = 0;
	bool intersect(const line<T>& l, bool& e) const { std::valarray<T> t; return this->intersect(l,t,e); };
};

// `compact_shape' base class
template<class T> class compact_shape : public shape<T>
{
private:
	std::valarray< const std::valarray<T>* > node_;
protected:
	void assign_nodes(const unsigned int ...);
	std::valarray<T> operator[] (const unsigned short&) const;
	//T measure_;
public:
	compact_shape() { };
	const std::valarray<T>* p(const unsigned short&) const;
	unsigned short dim() const;
	unsigned short size() const;
	std::string print() const;
	T measure() const;
	virtual T make_measure() const = 0;
	virtual T scale() const = 0;
	virtual const compact_shape<T>* face(const unsigned short&) const = 0;
	std::size_t adjacent(const compact_shape<T>&, const bool& = true) const;
};





/* ---------- function definitions ---------- */





/* ---------- assign_nodes ---------- */
template<class T> void compact_shape<T>::assign_nodes(const unsigned int nn ...)
{
	// get the nodes
	this->node_.resize(nn);
	std::va_list ap;
	va_start(ap, nn);
	for (unsigned short i=0; i<nn; i++)
	{
		this->node_[i] = va_arg(ap, const std::valarray<T>*);
	}
	va_end(ap);
	
	// check for consistency
	for (unsigned short i=0; i<nn; i++)
	{
		if ( (this->node_[i])->size() != (this->node_[(i+1)%nn])->size() ||
			v_eq(this->node_[i], this->node_[(i+1)%nn]) )
		{
			throw;
            return;
		}
	}
}





/* ---------- operator[] ---------- */
template<class T> inline std::valarray<T> compact_shape<T>::operator[](const unsigned short& n) const
{
	return *(this->node_[n]);
}





/* ---------- p ---------- */
template<class T> inline const std::valarray<T>* compact_shape<T>::p(const unsigned short& n) const
{
	return this->node_[n];
}





/* ---------- dim ---------- */
template<class T> inline unsigned short compact_shape<T>::dim() const
{
	return (this->node_[0])->size();
}





/* ---------- size ---------- */
template<class T> inline unsigned short compact_shape<T>::size() const
{
	return this->node_.size();
}





/* ---------- print ---------- */
template<class T> std::string compact_shape<T>::print() const
{
	std::ostringstream tmp;
	tmp << "{ ";
	for (unsigned short i=0; i<this->node_.size(); i++)
	{
		if ( i )
		{
			tmp << ", ";
		}
		tmp << "(";
		for (unsigned short j=0; j<(this->node_[i])->size(); j++)
		{
			if ( j )
			{
				tmp << ',';
			}
			tmp << (*(this->node_[i]))[j];
		}
		tmp << ")";
	}
	tmp << " }";
	return tmp.str();
}





/* ---------- measure ---------- */
template<class T> inline T compact_shape<T>::measure() const
{
	//return this->measure_;
	return this->make_measure();
}





/* ---------- adjacent ---------- */
template<class T> std::size_t compact_shape<T>::adjacent(const compact_shape<T>& s, const bool& quick) const
{
    std::size_t counter = 0;
	for (std::size_t i=0; i<this->node_.size(); i++)
	{
		for (std::size_t j=0; j<s.size(); j++)
		{
			if ( this->node_[i] == s.p(j) )
			{
				counter++;
				if ( quick )
				{
					return counter;
				}
			}
		}
	}
	return counter;
}





#endif /* _SHAPE_ */
