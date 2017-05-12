
#include <algorithm>
#include <functional>
#include <iostream>
#include <limits>
#include <numeric>
#include <typeinfo>
#include <valarray>
#include <vector>

template<class C> class Add {
  C val;
public:
  Add(C c) { val = c; }
  void operator() (C& c) const { c += val; }
};

template<class T> bool equal(const T& a, const T& b)
{
  T e = std::numeric_limits<float>::epsilon();
  std::cout << "  e = " << e << std::endl;
  if ( std::abs(b-a) <= e )
    {
      return true;
    }
  return false;
}

int main(void)
{
  std::cout << "\nABOUT EXTREMES\n--------------" << std::endl;
  std::cout << "    int infinity = " << std::numeric_limits<int>::infinity() << std::endl;
  float i = std::numeric_limits<float>::infinity();
  std::cout << "  float infinity = " << i << std::endl;

  float a = 345.23;
  std::cout << a << "/infinity = " << a/i << std::endl;

  std::cout << "  float epsilon = " << std::numeric_limits<float>::epsilon() << std::endl
	    << " double epsilon = " << std::numeric_limits<double>::epsilon() << std::endl
	    << "ldouble epsilon = " << std::numeric_limits<long double>::epsilon() << std::endl;

  std::cout << "\nABOUT VECTORS\n--------------" << std::endl;
  std::vector<float> v1(5), v2(5);

  std::cout << "   v1 =";
  v1[0] = 2; v1[1] = 1; v1[2] = 4; v1[3] = 2; v1[4] = 9;\
  for (int j=0; j<v1.size(); j++)
    {
      std::cout << ' ' << v1[j];
    }
  std::cout << std::endl
	    << "   v2 =";
  v2[0] = 2; v2[1] = 2; v2[2] = 8; v2[3] = 8; v2[4] = 6;
  for (int j=0; j<v2.size(); j++)
    {
      std::cout << ' ' << v2[j];
    }
  std::cout << std::endl;

  std::cout << "v1.v2 = ";
  float b = std::inner_product(v1.begin(), v1.end(), v2.begin(), float(0));
  std::cout << b << std::endl;

  std::for_each(v1.begin(), v1.end(), Add<float>(5));
  std::cout << " v1+5 =";
  for (int j=0; j<v1.size(); j++)
    {
      std::cout << ' ' << v1[j];
    }
  std::cout << std::endl;

  std::valarray<double> mv1(3); mv1[0] = 5.1; mv1[1] = 5.3; mv1[2] = 7.8;
  std::valarray<double> mv2 = mv1; mv2[2] = 8.0;
  std::cout << "\n  mv1   = " << typeid(mv1).name();
  for (unsigned short q=0; q<mv1.size(); q++) std::cout << ", " << mv1[q];
  std::cout << "\n  mv2   = " << typeid(mv2).name();
  for (unsigned short q=0; q<mv2.size(); q++) std::cout << ", " << mv2[q];
  std::valarray<bool> mv21 = ( mv2 == mv1);
  std::cout << "\nmv2-mv1 = " << typeid(mv21).name();
  for (unsigned short q=0; q<(mv21).size(); q++) std::cout << ", " << (mv21)[q];

  std::cout << "\nABOUT TYPES\n-----------" << std::endl;
  std::cout << "sizeof(int) = " << sizeof(int) << std::endl
	    << "typeid(int) = " << typeid(int).name() << std::endl
	    << "int::digits = " << std::numeric_limits<int>::digits << std::endl;

  std::cout << "\nABOUT ERRORS\n------------" << std::endl;
  double aa = 2.000005;
  double bb = 1.999998;
  std::cout << " aa = " << aa << std::endl
	    << " bb = " << bb << std::endl;
  if ( equal(aa,bb) )
    {
      std::cout << "aa == bb" << std::endl;
    }
  else
    {
      std::cout << "aa != bb" << std::endl;
    }
}
