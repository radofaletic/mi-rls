/**
 extra_math
 
 A collection of C++ routines for extra arithmetic functionionality.
 The vector functions available are:
 - dot product
 - cross product
 - norm
 - normalise
 There are also some useful template classes defining some scalar arithmetic functionals.
 
 Rado Faletic
 Department of Physics
 Faculty of Science
 Australian National University  ACT  0200
 Australia
 
 Rado.Faletic@anu.edu.au
 21st September 2004
 22nd April 2022, updated to C++20 
 */





#ifndef _EXTRA_MATH_
#define _EXTRA_MATH_





/* ---------- standard header files ---------- */
#include <algorithm>
#include <cmath>
#include <limits>
#include <valarray>
#include <vector>





/* ------------ user header files ------------ */
#include "front-end.h"





/* ---------- function declarations ---------- */

// float dot_product = dot(vec1,vec2);
template<class T> T dot(const std::valarray<T>&, const std::valarray<T>&);

// float vec_norm_squared = norm2(vec);
template<class T> T norm2(const std::valarray<T>&);

// float vec_norm = norm(vec);
template<class T> T norm(const std::valarray<T>&);

// float norm_of_vec2-vec1 = norm(vec1, vec2);
template<class T> T norm(const std::valarray<T>*, const std::valarray<T>*);

// float vec_norm = normalise(vec);
template<class T> T normalise(std::valarray<T>&);

// std::valarray<float> vec_cross = cross(vec1,vec2);
template<class T> std::valarray<T> cross(const std::valarray<T>&, const std::valarray<T>&);

// std::valarray<float> vec_cross_of_differences = cross(vec1,vec2,vec3);
template<class T> std::valarray<T> cross(const std::valarray<T>*, const std::valarray<T>*, const std::valarray<T>*);

// bool are_equal_within_limits = eq(a, b);
template<class T> bool eq(const T&, const T& = T(0), const T& = T(1));

// bool are_same = v_eq(v1, v2);
template<class T> bool v_eq(const std::valarray<T>&, const std::valarray<T>&, const T& = T(1));

// bool are_same = v_eq(v1, v2);
template<class T> bool v_eq(const std::valarray<T>*, const std::valarray<T>*, const T& = T(1));

// bool is_num_identity_vector = v_eq(v1, num);
template<class T> bool v_eq(const std::valarray<T>&, const T& = T(0), const T& = T(1));

// bool are_same_or_negatives = v_pm_eq(v1, v2);
template<class T> bool v_pm_eq(const std::valarray<T>*, const std::valarray<T>*, const T& = T(1));

// bool is_zero = is_zero(a);
template<class T> bool is_zero(const T&, const T& = T(1));

// bool is_zero = is_zero(v);
template<class T> bool is_zero(const std::valarray<T>&, const T& = T(1));

// void vertical_flip(myvec, size_of_row);
template<class T> void vertical_flip(std::valarray<T>&, const std::size_t& = 1);

// void horizontal_flip(myvec, size_of_col);
template<class T> void horizontal_flip(std::valarray<T>&, const std::size_t& = 1);

// remove_duplicates(vec);
template<class T> void remove_duplicates(std::vector<T>&);

// remove_duplicates(vec_of_vecs);
template<class T> void remove_duplicates(std::vector< std::valarray<T> >&);

// float distance = maximise_distance(vec_of_vecs);
template<class T> T maximise_distance(const std::vector< std::valarray<T> >&);

// rotate_matrix(Nrows, Ncols, image);
template<class T> void rotate_matrix(std::size_t&, std::size_t&, std::valarray<T>&, const bool& = true);

// int power = ipow(init_int, ipower);
template<class T> T ipow(const T&, const std::size_t&);





/* ---------- functionals (classes) ---------- */

/* ---------- decreaser ---------- */
template<class T> class decreaser
{
private:
    T pivot;
public:
    decreaser(const T& my_value = T(0)) : pivot(my_value) { };
    void operator()(T& value) const { if ( value >=pivot ) value--; };
};

/* ---------- increaser ---------- */
template<class T> class increaser
{
private:
    T pivot;
public:
    increaser(const T& my_value = T(0)) : pivot(my_value) { };
    void operator()(T& value) const { if ( value >= pivot ) value++; };
};

/* ---------- AddVector ---------- */
/* this is a functional to add one
 vector to another               */
template<class T> class AddVector
{
private:
    std::valarray<T> add_vec;
public:
    AddVector(const std::valarray<T>& new_vec) : add_vec(new_vec) { };
    void operator()(std::valarray<T>& my_vec) const { my_vec += add_vec; };
};

/* ---------- SubtractVector ---------- */
/* this is a functional to subtract one
 vector from another                  */
template<class T> class SubtractVector
{
private:
    std::valarray<T> sub_vec;
public:
    SubtractVector(const std::valarray<T>& new_vec) : sub_vec(new_vec) { };
    void operator()(std::valarray<T>& my_vec) const { my_vec -= sub_vec; };
};

/* ---------- Abs ---------- */
template<class T> void Abs(std::valarray<T>& my_vec) { my_vec = std::abs(my_vec); };





/* ---------- classes ---------- */

template<class T> class two_numbers
{
public:
    std::size_t n_integer;
    T n_real;
    two_numbers() : n_integer(0), n_real(T(0)) { };
    two_numbers(const std::size_t& nint, const T& nre) : n_integer(nint), n_real(nre) { };
    ~two_numbers() { };
};
template<class T> bool operator==(const two_numbers<T>& a, const two_numbers<T>& b) { return a.n_integer == b.n_integer; };
template<class T> bool operator!=(const two_numbers<T>& a, const two_numbers<T>& b) { return a.n_integer != b.n_integer; };
template<class T> bool operator<(const two_numbers<T>& a, const two_numbers<T>& b) { return a.n_integer < b.n_integer; };
template<class T> bool operator<=(const two_numbers<T>& a, const two_numbers<T>& b) { return a.n_integer <= b.n_integer; };
template<class T> bool operator>(const two_numbers<T>& a, const two_numbers<T>& b) { return a.n_integer > b.n_integer; };
template<class T> bool operator>=(const two_numbers<T>& a, const two_numbers<T>& b) { return a.n_integer >= b.n_integer; };





/* ---------- function definitions ---------- */





/* ---------- dot ---------- */
/// returns the dot product ( u . v )
template<class T> inline T dot(const std::valarray<T>& u, const std::valarray<T>& v)
{
    T dp = (u*v).sum();
    return ( is_zero(dp) ) ? T(0) : dp;
}





/* ---------- norm2 ---------- */
/// returns the norm squared |v|^2
template<class T> inline T norm2(const std::valarray<T>& v)
{
    return dot(v,v);
}





/* ---------- norm ---------- */
/// returns the norm ( |v| )
template<class T> inline T norm(const std::valarray<T>& v)
{
    return std::sqrt(norm2(v));
}





/* ---------- norm ---------- */
/// returns the norm ( |v2-v1| )
template<class T> inline T norm(const std::valarray<T>* v1, const std::valarray<T>* v2)
{
    return std::sqrt(norm2(std::valarray<T>(*v2 - *v1)));
}





/* ---------- normalise ---------- */
/// returns the unit vector ( v / |v| )
template<class T> inline T normalise(std::valarray<T>& v)
{
    // divide the vector by its norm
    T ip = norm2(v);
    if ( !ip )
    {
        return T(0);
    }
    else if ( eq(ip, T(1)) )
    {
        return T(1);
    }
    ip = std::sqrt(ip);
    v /= ip;
    return ip;
}





/* ---------- cross ---------- */
/// returns the cross product ( v1 ⨯ v2 )
template<class T> inline std::valarray<T> cross(const std::valarray<T>& v1, const std::valarray<T>& v2)
{
    std::valarray<T> x(T(0), v1.size());
    if ( v1.size() > 3 || v1.size() != v2.size() )
    {
        debug("cross","can't compute cross product of vectors");
        x.resize(0);
        throw; return x;
    }
    if ( x.size() == 3 )
    {
        x[0] = v1[1] * v2[2] - v1[2] * v2[1];
        x[1] = v1[2] * v2[0] - v1[0] * v2[2];
        x[2] = v1[0] * v2[1] - v1[1] * v2[0];
    }
    return x;
}





/* ---------- cross ---------- */
/// returns the cross product ( (v2-v1)x(v3-v1) )
template<class T> inline std::valarray<T> cross(const std::valarray<T>* v1, const std::valarray<T>* v2, const std::valarray<T>* v3)
{
    std::valarray<T> x(T(0), v1->size());
    if ( v1->size() > 3 || v1->size() != v2->size() || v2->size() != v3->size() )
    {
        debug("cross","can't compute cross product of p_vectors");
        x.resize(0);
        throw; return x;
    }
    if ( x.size() == 3 )
    {
        x[0] = ((*v2)[1]-(*v1)[1]) * ((*v3)[2]-(*v1)[2]) - ((*v2)[2]-(*v1)[2]) * ((*v3)[1]-(*v1)[1]);
        x[1] = ((*v2)[2]-(*v1)[2]) * ((*v3)[0]-(*v1)[0]) - ((*v2)[0]-(*v1)[0]) * ((*v3)[2]-(*v1)[2]);
        x[2] = ((*v2)[0]-(*v1)[0]) * ((*v3)[1]-(*v1)[1]) - ((*v2)[1]-(*v1)[1]) * ((*v3)[0]-(*v1)[0]);
    }
    return x;
}





/* ---------- eq ---------- */
/// determins if a and b are equal, allowing for numerical innacuracies
template<class T> inline bool eq(const T& a, const T& b, const T& e)
{
    return ( a == b || std::abs(b-a) <= std::abs(e) * std::numeric_limits<T>::epsilon() ) ? true : false;
}





/* ---------- v_eq ---------- */
template<class T> inline bool v_eq(const std::valarray<T>& u, const std::valarray<T>& v, const T& e)
{
    if ( !v.size() || u.size() != v.size() )
    {
        return false;
    }
    for (std::size_t i=0; i<u.size(); i++)
    {
        if ( !eq(u[i], v[i], e) )
        {
            return false;
        }
    }
    return true;
}





/* ---------- v_eq ---------- */
template<class T> inline bool v_eq(const std::valarray<T>* u, const std::valarray<T>* v, const T& e)
{
    if ( !v->size() || u->size() != v->size() )
    {
        return false;
    }
    for (std::size_t i=0; i<u->size(); i++)
    {
        if ( !eq((*u)[i], (*v)[i], e) )
        {
            return false;
        }
    }
    return true;
}





/* ---------- v_eq ---------- */
template<class T> inline bool v_eq(const std::valarray<T>& v, const T& n, const T& e)
{
    if ( !v.size() )
    {
        return false;
    }
    for (std::size_t i=0; i<v.size(); i++)
    {
        if ( !eq(v[i], n, e) )
        {
            return false;
        }
    }
    return true;
}





/* ---------- v_pm_eq ---------- */
template<class T> inline bool v_pm_eq(const std::valarray<T>* u, const std::valarray<T>* v, const T& e)
{
    if ( !v->size() || u->size() != v->size() )
    {
        return false;
    }
    bool p = true;
    bool m = true;
    for (std::size_t i=0; i<u->size(); i++)
    {
        if ( !eq((*u)[i], (*v)[i], e) )
        {
            p = false;
        }
        if ( !eq((*u)[i], -(*v)[i], e) )
        {
            m = false;
        }
    }
    return ( p || m );
}





/// determines if a and b and equal, allowing for numerical innacuracies
template<class T> inline bool is_zero(const T& a, const T& e)
{
    return eq(a, T(0), e);
}





/* ---------- is_zero ---------- */
template<class T> inline bool is_zero(const std::valarray<T>& v, const T& e)
{
    return v_eq(v, T(0), e);
}





/* ---------- vertical_flip ---------- */
template<class T> inline void vertical_flip(std::valarray<T>& v, const std::size_t& row_size)
{
    if ( v.size()%row_size )
    {
        throw; return;
    }
    for (std::size_t i=0; i<(v.size()/row_size)/2; i++)
    {
        for (std::size_t j=0; j<row_size; j++)
        {
            std::swap(v[i*row_size+j], v[v.size()-((i+1)*row_size)+j]);
        }
    }
}





/* ---------- horizontal_flip ---------- */
template<class T> inline void horizontal_flip(std::valarray<T>& v, const std::size_t& col_size)
{
    if ( v.size()%col_size )
    {
        throw; return;
    }
    for (std::size_t j=0; j<v.size()/col_size; j++)
    {
        std::reverse(&(v[j*col_size]), &(v[j*col_size+col_size]));
    }
}





/* ---------- remove_duplicates ---------- */
template<class T> void remove_duplicates(std::vector<T>& v)
{
    std::sort(v.begin(), v.end());
    typename std::vector<T>::iterator v_p = std::unique(v.begin(), v.end());
    v.erase(v_p, v.end());
}





/* ---------- remove_duplicates ---------- */
template<class T> void remove_duplicates(std::vector< std::valarray<T> >& v)
{
    bool erc = false;
    while ( !erc && v.size() > 1 )
    {
        for (std::size_t i=0; i<v.size(); i++)
        {
            if ( !v[i].size() )
            {
                v.erase(v.begin()+i);
                break;
            }
            for (std::size_t j=i+1; j<v.size(); j++)
            {
                if ( v_eq(v[i], v[j], T(100)) )
                    // an 'e' of 10^2 * std::epsilon<T> seems to work here
                {
                    v.erase(v.begin()+j);
                    erc = true;
                    break;
                }
            }
            if ( erc )
            {
                erc = false;
                break;
            }
            if ( i == v.size() - 1 )
            {
                erc = true;
            }
        }
    }
}





/* ---------- maximise_distance ---------- */
template<class T> T maximise_distance(const std::vector< std::valarray<T> >& v)
{
    T distance = T(0);
    for (std::size_t i=0; i<v.size(); i++)
    {
        for (std::size_t j=i+1; j<v.size(); j++)
        {
            distance = std::max(distance, norm(std::valarray<T>(v[j] - v[i])));
        }
    }
    return distance;
}





/* ---------- rotate_matrix ---------- */
/// rotate an image 90 degrees in the clockwise direction, or anti-clockwise if dir is false
template<class T> void rotate_matrix(std::size_t& Nrows, std::size_t& Ncols, std::valarray<T>& m, const bool& dir)
{
    std::valarray<T> orig = m;
    for (std::size_t i=0; i<Nrows; i++)
    {
        if ( !dir )
        {
            m[std::slice(Nrows-1-i,Ncols,Nrows)] = orig[std::slice(i*Ncols,Ncols,1)];
        }
        else
        {
            m[std::slice(i*Nrows,Nrows,1)] = orig[std::slice(Ncols-1-i,Nrows,Ncols)];
        }
    }
    std::swap(Nrows, Ncols);
}





/* ---------- ipow ---------- */
template<class T> T ipow(const T& base, const std::size_t& power)
{
    T ret = T(1);
    for (std::size_t i=0; i<power; i++)
    {
        ret *= base;
    }
    return ret;
}





#endif /* _EXTRA_MATH_ */
