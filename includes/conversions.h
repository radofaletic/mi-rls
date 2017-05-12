/*
  conversions

  convert various formats strings

  Rado Faletic
  Department of Physics
  Faculty of Science
  Australian National University  ACT  0200
  Australia

  Rado.Faletic@anu.edu.au
  8th July 2004
*/


#ifndef _CONVERSIONS_
#define _CONVERSIONS_


/* ---------- standard header files ---------- */
#include <algorithm>
#include <iomanip>
#include <sstream>
#include <string>
#include <valarray>
/* ------------------------------------------- */


/* ------------------------------------------- */
/* ---------- function declarations ---------- */
/* ------------------------------------------- */


template<class T> std::string ctos(const T&);

template<class T> std::string ntos(const T&);
long int stoi(const std::string&);

template<class T> std::string ntoIEEE(const T&);
double IEEEton(const std::string&);

template<class T> std::string vtoIEEE(const std::valarray<T>&);
template<class T> std::valarray<T> IEEEtov(const std::string&);

template<class T> std::string vtos(const T&);

template<class T> std::string vtos(const T&, const T&);

template<class T> std::string vtos(const T&, const T&, const T&);

bool is_nbo();

template<class T> void byte_swap(T&, const unsigned short& = sizeof(T));

std::string indent(const std::string&, const unsigned short& = 1);


/* ------------------------------------------ */
/* ---------- function definitions ---------- */
/* ------------------------------------------ */


/* ---------- ctos ---------- */
/* returns a string from a
   char                       */
template<class T> inline std::string
ctos(const T& character_string)
{
  return std::string(character_string);
}
/* -------------------------- */

/* ---------- ntos ---------- */
/* returns a string from a
   number                     */
template<class T> std::string
ntos(const T& number)
{
  std::ostringstream word;
  word << number;
  return word.str();
}
/* -------------------------- */

/* ---------- stoi ---------- */
/* returns an integer from a
   string                     */
long int
stoi(const std::string& number)
{
  std::istringstream word(number);
  long int ret;
  word >> ret;
  return ret;
}
/* -------------------------- */

/* ---------- ntoIEEE ---------- */
/* returns an IEEE string from a
   number                        */
template<class T> std::string
ntoIEEE(const T& number)
{
  short ntoi_precision = 8;
  std::ostringstream word;
  bool is_i = false;
  std::string bw = "";
  std::string ew = "";
  if ( std::floor(number) != std::ceil(number) )
    {
      word.setf(std::ios_base::scientific, std::ios_base::floatfield);
      word << std::setprecision(ntoi_precision) << number;
      bw = word.str();
      std::string::size_type get_e = bw.find("e");
      if ( get_e != bw.npos )
	{
	  bw = bw.substr(0,get_e);
	}
      for (size_t i=bw.size(); 0<i; i--)
	{
	  if ( bw[i-1] == '0' )
	    {
	      bw = bw.substr(0,bw.size()-1);
	    }
	  else if ( bw[i-1] == '.' )
	    {
	      bw = bw.substr(0,bw.size()-1);
	      break;
	    }
	  else
	    {
	      break;
	    }
	}
      ew = word.str();
      get_e = ew.find("e");
      if ( get_e != ew.npos )
	{
	  ew = ew.substr(get_e, ew.size()-get_e);
	  std::string::size_type get_pm = ew.find("+");
	  if ( get_pm == ew.npos )
	    {
	      get_pm = ew.find("-");
	    }
	  if ( get_pm != ew.npos && get_pm < ew.size()-1 )
	    {
	      get_pm++;
	      if ( ew.substr(get_pm,ew.size()-get_pm) == std::string("0") 
		   || ew.substr(get_pm,ew.size()-get_pm) == std::string("00") 
		   || ew.substr(get_pm,ew.size()-get_pm) == std::string("000") )
		{
		  ew = "";
		}
	    }
	}
      else
	{
	  ew = "";
	}
    }
  else
    {
      is_i = true;
      word << (long int)(number);
    }
  return ( is_i ) ? word.str() : bw+ew;
}
/* ----------------------------- */

/* ---------- IEEEton ---------- */
/* returns an number from an IEEE
   string                        */
double
IEEEton(const std::string& number)
{
  std::istringstream word(number);
  double ret;
  word >> ret;
  return ret;
}
/* ----------------------------- */

/* ---------- vtoIEEE ---------- */
/* returns an IEEE string from a
   vector of number              */
template<class T> std::string
vtoIEEE(const std::valarray<T>& v)
{
  std::ostringstream word;
  for (size_t i=0; i<v.size(); i++)
    {
      if ( i )
	{
	  word << ',';
	}
      word << ntoIEEE(v[i]);
    }
  return word.str();
}
/* ----------------------------- */

/* ---------- IEEEtov ---------- */
/* returns a vector from an IEEE
   string                        */
template<class T> std::valarray<T>
IEEEtov(const std::string& number)
{
  size_t nsz = std::count(number.begin(), number.end(), ',') + 1;
  std::valarray<T> ret(nsz);
  size_t fp = 0;
  size_t sp = 1;
  for (size_t i=0; i<ret.size(); i++)
    {
      for (size_t j=sp; j<number.size(); j++)
	{
	  if ( number[j] == ',' || j == number.size() - 1 )
	    {
	      sp = ( j == number.size() - 1 ) ? j+1 : j;
	      break;
	    }
	}
      ret[i] = IEEEton(number.substr(fp,sp-fp));
      fp = sp + 1;
      sp = fp + 1;
    }
  return ret;
}
/* ----------------------------- */

/* ---------- vtos ---------- */
/* returns a string from a
   vector/valarray            */
template<class T> std::string
vtos(const T& v)
{
  std::ostringstream word;
  word << '(';
  for (size_t i=0; i<v.size(); i++)
    {
      if (i != 0 )
	{
	  word << ',';
	}
      word << v[i];
    }
  word << ')';
  return word.str();
}
/* -------------------------- */

/* ---------- vtos ---------- */
/* returns a strings from two
   vectors/valarrays          */
template<class T> std::string
vtos(const T& u, const T& v)
{
  std::ostringstream word;
  word << '(';
  for (size_t i=0; i<u.size(); i++)
    {
      if (i != 0 )
	{
	  word << ',';
	}
      word << '{' << u[i] << ' ' << v[i] << '}';
    }
  word << ')';
  return word.str();
}
/* -------------------------- */

/* ---------- vtos ---------- */
/* returns a string from three
   vectors/valarrays          */
template<class T> std::string
vtos(const T& u, const T& v, const T& w)
{
  std::ostringstream word;
  word << '(';
  for (size_t i=0; i<u.size(); i++)
    {
      if (i != 0 )
	{
	  word << ',';
	}
      word << '{' << u[i] << ' ' << v[i] << ' ' << w[i] << '}';
    }
  word << ')';
  return word.str();
}
/* -------------------------- */

/* ---------- is_nbo ---------- */
/* true if system is in Network
   Byte Order (Big Endian), false
   otherwise                    */
bool
is_nbo()
{
  int mynumber = 67305985;
  char* c = (char*) &mynumber;
  int c1 = int(*(c+0));
  int c2 = int(*(c+1));
  int c3 = int(*(c+2));
  int c4 = int(*(c+3));
  int newnum = (c4*1) + (c3*256) + (c2*65536) + (c1*16777216);
  if ( newnum == mynumber )
    {
      return true;
    }
  return false;
}
/* ---------------------------- */

/* ---------- byte_swap ---------- */
/* given any value reverse the
   byte order, ie change the
   endiannes                       */
template<class T> void
byte_swap(T& number, const unsigned short& size_of)
{
  T* number_ref = &number;
  unsigned char* swapper = (unsigned char*) number_ref;
  for (unsigned short i=0; i<size_of/2; i++)
    {
      std::swap(swapper[i], swapper[size_of-1-i]);
    }
  return;
}

/* ------------------------------- */

/* ---------- indent ---------- */
/* indent the given string by
   `n' tabulatures              */
std::string
indent(const std::string& str, const unsigned short& n)
{
  std::string tabs("");
  for (unsigned short i=0; i<n; i++)
    {
      tabs += "\t";
    }
  std::string s = tabs + str;
  std::string::size_type pb = 0;
  while ( pb != s.size() )
    {
      if ( s[pb] == '\n' )
	{
	  s.replace(pb,1,"\n"+tabs);
	}
      pb++;
    }
  return s;
}
/* ---------------------------- */


#endif /* _CONVERSIONS_ */
