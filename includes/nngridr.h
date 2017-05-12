/*
  NNGRIDR - Natural Neighbour GRIDding algoRithm(c) 1991 D F Watson,
                                            2003,2004,2005 R Faletic
*/

#ifndef _NNGRIDR_
#define _NNGRIDR_


/* standard library includes */
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <limits>
#include <string>
#include <valarray>
#include <vector>


/* local includes */

#include "angles.h"
#include "conversions.h"
#include "front-end.h"


/* "NNG_datum" class */
template<class T> class NNG_datum                           
{
public:
  std::valarray<T> values;
  NNG_datum<T>*    nextdat;
  NNG_datum() : values(3), nextdat(0) { };
};

/* "NNG_simp" class */
template<class T> class NNG_simp
{
public:
  std::valarray<int> vert;
  std::valarray<T>   cent;
  NNG_simp<T>*       nextsimp;
  NNG_simp() : vert(3), cent(3), nextsimp(0) { };
};

/* "NNG_temp" class */
class NNG_temp
{
public:
  std::valarray<int> end;
  NNG_temp*          nexttemp;
  NNG_temp() : end(2), nexttemp(0) { };
};

/* "NNG_neig" class */
template<class T> class NNG_neig
{
public:
  int           neinum;
  T             narea;
  T             coord;
  NNG_neig<T>*  nextneig;
  NNG_neig() : nextneig(0) { };
};

/* "NNG_Matrix" class */
template<class T> class NNG_Matrix
{
private:
  size_t            Mrows_;
  size_t            Mcols_;
  std::valarray<T>* Mdata_;
public:
  NNG_Matrix() { Mrows_ = 0; Mcols_ = 0; Mdata_ = 0; };
  NNG_Matrix(const size_t& nrows, const size_t& ncols, const T& value = T(0)) { init(nrows, ncols, value); };
  ~NNG_Matrix() { free(); };
  void init(const size_t& nrows, const size_t& ncols, const T& value = T(0))
  { 
    Mrows_ = ( nrows < 2 ) ? 2 : nrows;
    Mcols_ = ( ncols < 2 ) ? 2 : ncols;
    Mdata_ = new std::valarray<T>(value, Mrows_*Mcols_);
  };
  void free() { Mrows_ = 0; Mcols_ = 0; delete Mdata_; };
  T operator()(const size_t& row, const size_t& col) const { return (*Mdata_)[row*Mcols_+col]; };
  T& operator()(const size_t& row, const size_t& col) { return (*Mdata_)[row*Mcols_+col]; };
};

/* "nngridr_ini" class */
template<class T> class nngridr_ini
{
public:
  // variables from "nngridr.ini"
  bool        igrad;
  bool        ichoro;
  int         x_nodes;
  int         y_nodes;
  bool        non_neg;
  short       arriba;
  T           bI;
  T           bJ;
  T           xstart;
  T           ystart;
  T           xterm;
  T           yterm;
  T           nuldat;
  bool        extrap;
  // internal variables
  bool        updir;
  int         magcnt;
  bool        mxmn;
  T           xspace;
  T           yspace;
  // functions
  void        ShowState(const bool&, const bool& = false);
  // initialiser
  nngridr_ini()
  {
    igrad   = false; // estimate gradients?
    ichoro  = false; // calculate Choropleth data?
    x_nodes = 10; // number of points in the x direction
    y_nodes = 10; // number of points in the y direction
    non_neg = false; // if non_neg=true then don't allow negative interpolation
    arriba  = 1; // +1=L_to_R, -1=R_to_Leog ou
    bI      = 1.5;
    bJ      = 7.0;
    xstart  = T(0); // smallest x
    ystart  = T(0); // smallest y
    xterm   = T(1); // largest x
    yterm   = T(1); // largest y
    nuldat  = T(-999); // what number represents a NULL data value
    extrap  = true; // extrapolation outside of the convex hull
    mxmn    = false; // write out max/min
    updir   = true; // if updir=true write S_to_N
    xspace  = (xterm - xstart) / (x_nodes - 1); // spacing in the x direction
    yspace  = (yterm - ystart) / (y_nodes - 1); // spacing in the y direction
  };
};

/* "nngrid" class */
template<class T> class nngrid
{
public:
  // variables
  NNG_Matrix<T>       points;
  NNG_Matrix<T>       joints;
  std::valarray<int>  jndx;
  NNG_datum<T>*       rootdat;
  NNG_datum<T>*       curdat;
  NNG_simp<T>*        rootsimp;
  NNG_simp<T>*        cursimp;
  NNG_simp<T>*        lastsimp;
  NNG_simp<T>*        prevsimp;
  NNG_temp*           roottemp;
  NNG_temp*           curtemp;
  NNG_temp*           lasttemp;
  NNG_temp*           prevtemp;
  NNG_neig<T>*        rootneig;
  NNG_neig<T>*        curneig;
  NNG_neig<T>*        lastneig;
  std::valarray<T>    bigtri;
  int                 datcnt;
  int                 datcnt3;
  int                 numtri;
  int                 numnei;
  int                 neicnt;
  bool                goodflag;
  bool                ext;
  T                   wbit;
  T                   asum;
  T                   xx;
  T                   maxhoriz;
  std::valarray<T>    maxxy;
  std::valarray<T>    work3;
  T                   sumx;
  T                   sumy;
  T                   sumz;
  T                   sumx2;
  T                   sumy2;
  T                   sumxy;
  T                   sumyz;
  T                   sumzx;
  T                   bbb;
  T                   ccc;
  T                   bignum;
  T                   range;
  // initialiser
  nngrid() : rootdat(0),
	     rootsimp(0),
	     roottemp(0),
	     rootneig(0),
	     bigtri(9),
	     datcnt(0),
	     maxxy(6),
	     work3(9),
	     range(T(10))
  {
    bigtri[0] = T(-1); bigtri[1] = T(-1); bigtri[2] = T(0);
    bigtri[3] = T(5);  bigtri[4] = T(-1); bigtri[5] = T(0);
    bigtri[6] = T(-1); bigtri[7] = T(5);  bigtri[8] = T(0);
    bignum    = std::pow(T(10), std::numeric_limits<T>::max_exponent10 - 1);
  };
};


/* function declarations from nnchead.h */
template<class T> bool ReadData(nngrid<T>&, nngridr_ini<T>&, std::vector< std::valarray<T> >&, std::valarray<T>&, const bool& = false);
template<class T> void ChoroPleth(nngrid<T>&, nngridr_ini<T>&);
template<class T> void FindNeigh(nngrid<T>&, const int&);
template<class T> void TriCentr(nngrid<T>&);
template<class T> void TriNeigh(nngrid<T>&);
template<class T> void Gradient(nngrid<T>&);
template<class T> bool MakeGrid(nngrid<T>&, const nngridr_ini<T>&, std::vector< std::valarray<T> >&, std::valarray<T>&);
template<class T> void FindProp(nngrid<T>&, const T&, const T&);
template<class T> T Surface(nngrid<T>&);
template<class T> T Meld(nngrid<T>&, const nngridr_ini<T>&, T, const T&, const T&);


/* local function declarations */
template<class T> void nngridrc(const std::vector< std::valarray<T> >&, std::valarray<T>&, const size_t&, const size_t&, const bool& = false);
template<class T> void nngridr(std::vector< std::valarray<T> >&, std::valarray<T>&, const size_t&, const size_t&, const bool& = false);


template<class T> void
nngridrc(const std::vector< std::valarray<T> >& ucpoints,
	 std::valarray<T>& udata,
	 const size_t& Nrows,
	 const size_t& Ncols,
	 const bool& verbose)
{
  std::vector< std::valarray<T> > upoints = ucpoints;
  nngridr(upoints, udata, Nrows, Ncols, verbose);
}

template<class T> void
nngridr(std::vector< std::valarray<T> >& upoints,
	std::valarray<T>& udata,
	const size_t& Nrows,
	const size_t& Ncols,
	const bool& verbose)
{
  if ( upoints.size() != udata.size() )
    {
      debug("nngridr","number of points = "+ntos(upoints.size())+", but number of datum = "+ntos(udata.size()));
      throw; return;
    }
  nngridr_ini<T> ini;
  nngrid<T>      nng;

  ini.x_nodes = Ncols;
  ini.y_nodes = Nrows;
  ini.xstart = ini.xterm = upoints[0][0];
  ini.ystart = ini.yterm = upoints[0][1];
  for (size_t i=0; i<upoints.size(); i++)
    {
      ini.xstart = std::min(upoints[i][0], ini.xstart);
      ini.xterm = std::max(upoints[i][0], ini.xterm);
      ini.ystart = std::min(upoints[i][1], ini.ystart);
      ini.yterm = std::max(upoints[i][1], ini.yterm);
    }
  T ini_scale = ( ini.xterm - ini.xstart + ini.yterm - ini.ystart ) / 2;
  if ( ini_scale < 2 )
    {
      ini_scale = 2;
    }
  //ini.xstart -= ini_scale * std::numeric_limits<T>::epsilon();
  //ini.xterm += ini_scale * std::numeric_limits<T>::epsilon();
  //ini.ystart -= ini_scale * std::numeric_limits<T>::epsilon();
  //ini.yterm += ini_scale * std::numeric_limits<T>::epsilon();
  //ini.xstart -= 10E-5;
  //ini.xterm += 10E-5;
  //ini.ystart -= 10E-5;
  //ini.yterm += 10E-5;
  ini.xspace = ( ini.xterm - ini.xstart ) / ( ini.x_nodes - 1 );
  ini.yspace = ( ini.yterm - ini.ystart ) / ( ini.y_nodes - 1 );
  if ( ini.yspace < ini.xspace )
    {
      T extra = ( ( ini.xterm - ini.xstart ) - ( ini.yterm - ini.ystart ) ) / 2;
      ini.ystart -= extra;
      ini.yterm += extra;
      ini.yspace = ( ini.yterm - ini.ystart ) / ( ini.y_nodes - 1 );
    }
  else if ( ini.xspace < ini.yspace )
    {
      T extra = ( ( ini.yterm - ini.ystart ) - ( ini.xterm - ini.xstart ) ) / 2;
      ini.xstart -= extra;
      ini.xterm += extra;
      ini.xspace = ( ini.xterm - ini.xstart ) / ( ini.x_nodes - 1 );
    }

  ini.ShowState(false, verbose);
  if ( verbose )
    {
      debug("nngridr","..................... Execution Proceeding ......................");
    }
  if ( ReadData(nng, ini, upoints, udata, verbose) )
    {
      if ( verbose )
	{
	  debug(ntos(nng.datcnt)+" data were input");
	}
      if ( ini.ichoro )
	{
	  ChoroPleth(nng, ini);
	  if ( verbose )
	    {
	      debug("Choropleth calculations completed");
	    }
	}
      if ( ini.igrad )
	{
	  Gradient(nng);
	  if ( verbose )
	    {
	      debug("Gradient estimation completed");
	    }
	}
      MakeGrid(nng, ini, upoints, udata);
      if ( verbose )
	{
	  debug("nngridr","Grid completed");
	}
    }
}


template<class T> void nngridr_ini<T>::ShowState(const bool& atime, const bool& verbose)
{
  if ( verbose )
    {
      debug("nngridr::ShowState","NNGRIDR - Natural Neighbour GRIDding algoRithm (c) 1991 D F Watson");
      debug("                                                   2003 R Faletic");
    }
  magcnt = 0;
  if ( mxmn )
    {
      if ( verbose )
	{
	  debug("Print data maximums and minimums");
	}
      magcnt++;
    }
  if ( verbose )
    {
      if ( atime )
	{
	  debug("Print execution times");
	}
      debug("..................... Interpolation Controls .....................");
      debug(((ichoro)?std::string("Choropleth"):std::string("Functional"))+" interpolation");
      if ( igrad )
	{
	  debug("Using estimated gradients");
	  debug("Tautness parameter #1 "+ntos(bI));
	  debug("Tautness parameter #2 "+ntos(bJ));
	  if ( non_neg )
	    {
	      debug("Negative interpolation not allowed");
	    }
	  else
	    {
	      debug("Negative interpolation is allowed");
	    }
	}
      else
	{
	  debug("Don't use gradients");
	}
      if ( extrap )
	{
	  debug("Extrapolation allowed");
	}
      else 
	{
	  debug("No extrapolation");
	  debug("Null value is "+ntos(nuldat));
	}
      debug("................... Gridded Region Parameters ....................");
      debug("Nodes along x: "+ntos(x_nodes));
      debug("Nodes along y: "+ntos(y_nodes));
      
      debug("Least x: "+ntos(xstart));
      debug("Least y: "+ntos(ystart));
      
      debug("x spacing: "+ntos(xspace));
      debug("y spacing: "+ntos(yspace));
      
      debug("Greatest x: "+ntos(xterm));
      debug("Greatest y: "+ntos(yterm));
      debug("....................... Output Parameters ........................");
      if ( updir )
	{
	  debug("Output order from south to north");
	}
      else
	{
	  debug("Output order from north to south");
	}
    }
}


template<class T> bool ReadData(nngrid<T>& nng,
				nngridr_ini<T>& ini,
				std::vector< std::valarray<T> >& upoints,
				std::valarray<T>& udata,
				const bool& verbose)
{
  std::valarray<T> temp(3); // double temp[3];

  if ( ! nng.rootdat ) 
    {
      nng.rootdat = new NNG_datum<T>; // IMakeDatum();
      nng.rootsimp = new NNG_simp<T>; // IMakeSimp();
      nng.roottemp = new NNG_temp; // IMakeTemp();
      nng.rootneig = new NNG_neig<T>; // IMakeNeig();
      (nng.rootdat)->values[0] = (nng.rootdat)->values[1]
	= (nng.rootdat)->values[2] = T(0);
    }
  else 
    {
      nng.jndx.resize(0);
      nng.points.free();
      nng.joints.free();
    }
  nng.curdat = nng.rootdat;
  nng.datcnt = 0;
  if ( ini.mxmn )
    {
      nng.maxxy[0] = nng.maxxy[1] = nng.maxxy[2] = -(nng.maxxy[3] = nng.maxxy[4] = nng.maxxy[5] = nng.bignum);
    }
  T minx = ini.xstart;
  T maxx = ini.xterm;
  T miny = ini.ystart;
  T maxy = ini.yterm;
  for (size_t i=0; i<upoints.size(); i++)
    {
      temp[0] = upoints[i][0];
      temp[1] = upoints[i][1];
      temp[2] = udata[i];
      if ( ini.mxmn )
	{
	  for (int i1=0; i1<3; i1++)
	    {
	      if ( nng.maxxy[i1] < temp[i1] ) 
		{
		  nng.maxxy[i1] = temp[i1]; 
		}
	      if ( nng.maxxy[3+i1] > temp[i1] ) 
		{
		  nng.maxxy[3+i1] = temp[i1]; 
		}
	    }
	}
      if ( minx < temp[0] && temp[0] < maxx && miny < temp[1] && temp[1] < maxy )
	{
	  if ( ! (nng.curdat)->nextdat ) 
	    {
	      (nng.curdat)->nextdat = new NNG_datum<T>; // IMakeDatum();
	    }
	  nng.curdat = (nng.curdat)->nextdat;
	  nng.datcnt++;
	  (nng.curdat)->values = temp;
	}
    }
  if ( verbose && ini.mxmn )
    {
      debug("Data maximums "+ntoIEEE(nng.maxxy[0])+"  "+ntoIEEE(nng.maxxy[1])+"  "+ntoIEEE(nng.maxxy[2]));
      debug("Data minimums "+ntoIEEE(nng.maxxy[3])+"  "+ntoIEEE(nng.maxxy[4])+"  "+ntoIEEE(nng.maxxy[5]));
      debug("Data differ's "+ntoIEEE(nng.maxxy[0] - nng.maxxy[3])+"  "+ntoIEEE(nng.maxxy[1] - nng.maxxy[4])+"  "+ntoIEEE(nng.maxxy[2] - nng.maxxy[5]));
    }
  if ( nng.datcnt > 3 )
    {
      nng.datcnt3 = nng.datcnt + 3;
      nng.jndx.resize(nng.datcnt3); // IntVect(datcnt3);
      nng.sumx = nng.sumy = nng.sumz
	= nng.sumx2 = nng.sumy2
	= nng.sumxy = nng.sumzx = nng.sumyz
	= T(0);
      nng.maxxy[0] = nng.maxxy[1] = nng.maxxy[2]
	= -(nng.maxxy[3] = nng.maxxy[4] = nng.maxxy[5] = nng.bignum);
      if ( ini.igrad )
	{
	  nng.points.init(nng.datcnt+4, 6);
	}
      else
	{
	  nng.points.init(nng.datcnt+4, 3);
	}
      nng.joints.init(nng.datcnt3, 2); 
      nng.curdat = (nng.rootdat)->nextdat;
      (nng.rootdat)->nextdat = 0;
      for (int i0=0; i0<nng.datcnt; i0++)
	{
	  nng.sumx += nng.points(i0,0) = (nng.curdat)->values[0];
	  nng.sumx2 += nng.points(i0,0) * nng.points(i0,0);
	  if  (nng.maxxy[0] < nng.points(i0,0) ) 
	    {
	      nng.maxxy[0] = nng.points(i0,0);  
	    }
	  if ( nng.maxxy[3] > nng.points(i0,0) ) 
            {
	      nng.maxxy[3] = nng.points(i0,0);  
	    }
	  nng.sumy += nng.points(i0,1) = (nng.curdat)->values[1];
	  nng.sumy2 += nng.points(i0,1) * nng.points(i0,1);
	  nng.sumxy += nng.points(i0,0) * nng.points(i0,1);
	  if ( nng.maxxy[1] < nng.points(i0,1) ) 
	    {
	      nng.maxxy[1] = nng.points(i0,1);  
	    }
	  if ( nng.maxxy[4] > nng.points(i0,1) ) 
            {
	      nng.maxxy[4] = nng.points(i0,1);  
	    }
	  nng.sumz += nng.points(i0,2) = (nng.curdat)->values[2];
	  nng.sumzx += nng.points(i0,2) * nng.points(i0,0);
	  nng.sumyz += nng.points(i0,1) * nng.points(i0,2);
	  if ( nng.maxxy[2] < nng.points(i0,2) ) 
	    {
	      nng.maxxy[2] = nng.points(i0,2); 
	    }
	  if ( nng.maxxy[5] > nng.points(i0,2) ) 
	    {
	      nng.maxxy[5] = nng.points(i0,2); 
	    }
	  NNG_datum<T>* holddat = nng.curdat;
	  nng.curdat = (nng.curdat)->nextdat;
	  delete holddat;
	}
      T det = (nng.datcnt * (nng.sumx2 * nng.sumy2 - nng.sumxy * nng.sumxy)) -
	(nng.sumx * (nng.sumx * nng.sumy2 - nng.sumy * nng.sumxy)) +
	(nng.sumy * (nng.sumx * nng.sumxy - nng.sumy * nng.sumx2));
      T aaa = ((nng.sumz * (nng.sumx2 * nng.sumy2 - nng.sumxy * nng.sumxy)) -
	       (nng.sumzx * (nng.sumx * nng.sumy2 - nng.sumy * nng.sumxy)) +
	       (nng.sumyz * (nng.sumx * nng.sumxy - nng.sumy * nng.sumx2))) / det;
      nng.bbb = ((nng.datcnt * (nng.sumzx * nng.sumy2 - nng.sumyz * nng.sumxy)) -
		 (nng.sumz * (nng.sumx * nng.sumy2 - nng.sumy * nng.sumxy)) +
		 (nng.sumy * (nng.sumx * nng.sumyz - nng.sumy * nng.sumzx))) / det;
      nng.ccc = ((nng.datcnt * (nng.sumx2 * nng.sumyz - nng.sumxy * nng.sumzx)) -
		 (nng.sumx * (nng.sumx * nng.sumyz - nng.sumy * nng.sumzx)) +
		 (nng.sumz * (nng.sumx * nng.sumxy - nng.sumy * nng.sumx2))) / det;
      if ( verbose && ini.mxmn )
	{
	  debug("Grid region maximums "+ntoIEEE(nng.maxxy[0])+" "+ntoIEEE(nng.maxxy[1])+" "+ntoIEEE(nng.maxxy[2]));
	  debug("Grid region minimums "+ntoIEEE(nng.maxxy[3])+" "+ntoIEEE(nng.maxxy[4])+" "+ntoIEEE(nng.maxxy[5]));
	  debug("Grid region differ's "+ntoIEEE(nng.maxxy[0] - nng.maxxy[3])+" "+ntoIEEE(nng.maxxy[1] - nng.maxxy[4])+" "+ntoIEEE(nng.maxxy[2] - nng.maxxy[5]));
	}
      if ( nng.maxxy[0] < maxx ) 
	{
	  nng.maxxy[0] = maxx; 
	}
      if ( nng.maxxy[3] > minx ) 
	{
	  nng.maxxy[3] = minx; 
	}
      if ( nng.maxxy[1] < maxy ) 
	{
	  nng.maxxy[1] = maxy; 
	}
      if ( nng.maxxy[4] > miny ) 
	{
	  nng.maxxy[4] = miny; 
	}
      for (int i0=0; i0<3; i0++) 
	{
	  nng.maxxy[i0] -= nng.maxxy[3+i0];
	}
      nng.maxhoriz = nng.maxxy[0]; 
      if ( nng.maxhoriz < nng.maxxy[1] ) 
	{
	  nng.maxhoriz = nng.maxxy[1];
	}
      //nng.wbit = nng.maxhoriz * std::numeric_limits<T>::epsilon();
      nng.wbit = nng.maxhoriz * 10E-5;
      if ( nng.maxxy[0] / nng.maxxy[1] > T(2) || nng.maxxy[1] / nng.maxxy[0] > T(2) )
	{
	  ini.igrad = false;
	  debug("nngridr::ReadData","W A R N I N G  The ratio of width to length of this gridded region may be too extreme for good interpolation. Changing the block proportions, or rescaling the x or y coordinate may be indicated.");
	}
      if ( !ini.ichoro && ini.igrad )
	{
	  if ( nng.maxxy[2] / nng.maxxy[0] > T(60) || nng.maxxy[2] / nng.maxxy[1] > T(60) ) 
	    {
	      ini.igrad = false;
	      debug("nngridr::ReadData","N.B. The ratio of vertical to horizontal scales is too large for gradient estimation. Rescale the data if gradients are required.");
	    }
	  if ( nng.maxxy[2] / nng.maxxy[0] < T(.017) || nng.maxxy[2] / nng.maxxy[1] < T(.017) ) 
	    {
	      ini.igrad = false;
	      debug("nngridr::ReadData","N.B. The ratio of vertical to horizontal scales is too small for gradient estimation. Rescale the data if gradients are required.");
	    }
	}
      for (int i0=0; i0<3; i0++)
	{
	  nng.points(nng.datcnt+i0,0) = nng.maxxy[3] + nng.bigtri[3*i0] * nng.maxxy[0] * nng.range;
	  nng.points(nng.datcnt+i0,1) = nng.maxxy[4] + nng.bigtri[3*i0+1] * nng.maxxy[1] * nng.range;
	  nng.points(nng.datcnt+i0,2) = aaa + nng.bbb * nng.points(nng.datcnt+i0,0) + nng.ccc * nng.points(nng.datcnt+i0,1);
	}
      nng.rootdat = 0;
    }
  else
    {
      debug("nngridr::ReadData","Insufficient data in gridded region to triangulate... increase the size of the gridded region");
      throw; return false;
    }
  std::srand(367);     
  for (int i0=0; i0<nng.datcnt; i0++) 
    {
      for (int i1=0; i1<2; i1++)
	{
	  nng.points(i0,i1) += nng.wbit * (0.5 - (T)std::rand() / RAND_MAX);
	}
    }
  return true;
}


template<class T> void ChoroPleth(nngrid<T>& nng,
				  nngridr_ini<T>& ini)
{
  std::valarray<T> cmax(4); // double cmax[2][2];
  nng.maxxy[2] = -nng.bignum;
  nng.maxxy[5] = nng.bignum;
  nng.sumz = nng.sumzx = nng.sumyz = T(0);
  for (int i0=0; i0<nng.datcnt; i0++)
    {
      FindNeigh(nng, i0);
      cmax[0] = cmax[1] = -nng.bignum; // cmax[0][0] = cmax[0][1] = -nng.bignum;
      cmax[2] = cmax[3] = nng.bignum; //cmax[1][0] = cmax[1][1] = nng.bignum;
      nng.cursimp = nng.rootsimp;
      for (int i1=0; i1<nng.numtri; i1++)
	{
	  nng.cursimp = (nng.cursimp)->nextsimp;
	  nng.joints(i1,0) = (nng.cursimp)->cent[0];
	  if ( cmax[0] < nng.joints(i1,0) ) // if ( cmax[0][0] < joints[i1][0] ) 
	    {
	      cmax[0] = nng.joints(i1,0); // cmax[0][0] = joints[i1][0];
	    }
	  if ( cmax[2] > nng.joints(i1,0) ) // if ( cmax[1][0] > joints[i1][0] ) 
            {
	      cmax[2] = nng.joints(i1,0); // cmax[1][0] = joints[i1][0];
	    }
	  nng.joints(i1,1) = (nng.cursimp)->cent[1];
	  if ( cmax[1] < nng.joints(i1,1) ) // if ( cmax[0][1] < joints[i1][1] ) 
	    {
	      cmax[1] = nng.joints(i1,1); // cmax[0][1] = joints[i1][1];
	    }
	  if ( cmax[3] > nng.joints(i1,1) ) // if ( cmax[1][1] > joints[i1][1] ) 
	    {
	      cmax[3] = nng.joints(i1,1); // cmax[1][1] = joints[i1][1];
	    }
	}
      cmax[0] -= cmax[2]; // cmax[0][0] -= cmax[1][0];
      cmax[1] -= cmax[3]; // cmax[0][1] -= cmax[1][1];
      for (int i1=0; i1<3; i1++)
	{
	  nng.joints(nng.datcnt+i1,0) = cmax[2] + nng.bigtri[3*i1] * cmax[0] * nng.range; // joints[datcnt+i1][0] = cmax[1][0] + bigtri[i1][0] * cmax[0][0] * nng.range;
	  nng.joints(nng.datcnt+i1,1) = cmax[3] + nng.bigtri[3*i1+1] * cmax[1] * nng.range; // joints[datcnt+i1][1] = cmax[1][1] + bigtri[i1][1] * cmax[0][1] * nng.range;
	}
      nng.neicnt = nng.numtri;
      TriCentr(nng);
      nng.cursimp = nng.rootsimp;
      nng.asum = T(0);
      for (int i1=0; i1<nng.numtri; i1++)
	{
	  nng.cursimp = (nng.cursimp)->nextsimp;
	  if ( (nng.cursimp)->vert[0] < nng.datcnt )
	    {
	      nng.asum += std::abs( (nng.joints((nng.cursimp)->vert[1],0) - 
				      nng.joints((nng.cursimp)->vert[0],0)) *
				     (nng.joints((nng.cursimp)->vert[2],1) - 
				      nng.joints((nng.cursimp)->vert[0],1)) -
				     (nng.joints((nng.cursimp)->vert[2],0) - 
				      nng.joints((nng.cursimp)->vert[0],0)) *
				     (nng.joints((nng.cursimp)->vert[1],1) - 
				      nng.joints((nng.cursimp)->vert[0],1)) ) / 2;
	    }
	}
      if ( nng.asum > 0 )
	{
	  nng.points(i0,2) /= nng.asum; // conditional added 4/4/95
	}
      else
	{
	  nng.points(i0,2) = 0;
	}
      nng.sumz += nng.points(i0,2);
      nng.sumzx += nng.points(i0,2) * nng.points(i0,0);
      nng.sumyz += nng.points(i0,1) * nng.points(i0,2);
      if ( nng.maxxy[2] < nng.points(i0,2) ) 
	{
	  nng.maxxy[2] = nng.points(i0,2); 
	}
      if ( nng.maxxy[5] > nng.points(i0,2) ) 
	{
	  nng.maxxy[5] = nng.points(i0,2); 
	}
    }
  nng.maxxy[2] -= nng.maxxy[5];
  if ( ini.igrad )
    {
      if ( nng.maxxy[2] / nng.maxhoriz > T(60) ||
	   nng.maxxy[2] / nng.maxxy[1] > T(60) ) 
	{
	  ini.igrad = false;
	  debug("nngridr::ChoroPleth","N.B. The ratio of vertical to horizontal scales is too large for gradient estimation. Rescale the data if gradients are required.");
	}
      if ( nng.maxxy[2] / nng.maxhoriz < T(.017) ||
	   nng.maxxy[2] / nng.maxxy[1] < T(.017) ) 
	{
	  ini.igrad = false;
	  debug("nngridr::ChoroPleth","N.B. The ratio of vertical to horizontal scales is too small for gradient estimation. Rescale the data if gradients are required.");
	}
    }
  T det = (nng.datcnt * (nng.sumx2 * nng.sumy2 - nng.sumxy * nng.sumxy)) -
    (nng.sumx * (nng.sumx * nng.sumy2 - nng.sumy * nng.sumxy)) +
    (nng.sumy * (nng.sumx * nng.sumxy - nng.sumy * nng.sumx2));
  T aaa = ((nng.sumz * (nng.sumx2 * nng.sumy2 - nng.sumxy * nng.sumxy)) -
	   (nng.sumzx * (nng.sumx * nng.sumy2 - nng.sumy * nng.sumxy)) +
	   (nng.sumyz * (nng.sumx * nng.sumxy - nng.sumy * nng.sumx2))) / det;
  nng.bbb = ((nng.datcnt * (nng.sumzx * nng.sumy2 - nng.sumyz * nng.sumxy)) -
	     (nng.sumz * (nng.sumx * nng.sumy2 - nng.sumy * nng.sumxy)) +
	     (nng.sumy * (nng.sumx * nng.sumyz - nng.sumy * nng.sumzx))) / det;
  nng.ccc = ((nng.datcnt * (nng.sumx2 * nng.sumyz - nng.sumxy * nng.sumzx)) -
	     (nng.sumx * (nng.sumx * nng.sumyz - nng.sumy * nng.sumzx)) +
	     (nng.sumz * (nng.sumx * nng.sumxy - nng.sumy * nng.sumx2))) / det;
  for (int i1=0; i1<3; i1++)
    {
      nng.points(nng.datcnt+i1,2) =
	aaa + nng.bbb * nng.points(nng.datcnt+i1,0) +
	nng.ccc * nng.points(nng.datcnt+i1,1);
      if ( ini.igrad )
	{
	  nng.points(nng.datcnt+i1,3) = -nng.bbb;
	  nng.points(nng.datcnt+i1,4) = -nng.ccc;
	  nng.points(nng.datcnt+i1,5) = T(1);
	}
    }
}


template<class T> void Gradient(nngrid<T>& nng)
{
  for (int i0=0; i0<nng.datcnt; i0++)
    {
      FindNeigh(nng, i0);
      if ( !nng.ext ) 
	{
	  TriNeigh(nng);
	  T wxd = nng.points(i0,0);
	  T wyd = nng.points(i0,1);
	  FindProp(nng, wxd, wyd);
	  T xc = Surface(nng);
	  T wxde = wxd + nng.wbit;
	  FindProp(nng, wxde, wyd);
	  T xe = Surface(nng);
	  T wydn = wyd + nng.wbit;
	  FindProp(nng, wxd, wydn);
	  T xn = Surface(nng);
	  nng.points(i0,3) = (xc - xe) / nng.wbit;
	  nng.points(i0,4) = (xc - xn) / nng.wbit;
	  nng.asum /= Angle::pi; 
	  nng.points(i0,5) =
	    1 - std::sqrt(nng.asum / (nng.asum + (nng.points(i0,2) - xc) *
				      (nng.points(i0,2) - xc) ) );
	}
      else     
	{
	  nng.points(i0,3) = nng.points(i0,4) = nng.points(i0,5) = nng.xx = T(0);
	  nng.cursimp = nng.rootsimp;
	  for (int i1=0; i1<nng.numtri; i1++)
	    {
	      nng.cursimp = (nng.cursimp)->nextsimp;
	      for (int i2=0; i2<2; i2++) 
		{
		  for (int i3=0; i3<3; i3++)
		    {
		      nng.work3[3*i2+i3] = nng.points((nng.cursimp)->vert[0],i3) - nng.points((nng.cursimp)->vert[i2+1],i3);
		    }
		}
	      nng.work3[6] = nng.work3[1] * nng.work3[5] - nng.work3[4] * nng.work3[2];
	      nng.work3[7] = nng.work3[2] * nng.work3[3] - nng.work3[5] * nng.work3[0];
	      nng.work3[8] = nng.work3[0] * nng.work3[4] - nng.work3[3] * nng.work3[1];
	      T u2 = 1;
	      if ( nng.work3[8] < T(0) )
		{
		  u2 = -1;
		}
	      nng.xx += std::sqrt(nng.work3[6]*nng.work3[6] + nng.work3[7]*nng.work3[7] + nng.work3[8]*nng.work3[8]);
	      for (int i2=0; i2<3; i2++)
		{
		  nng.points(i0,i2+3) += nng.work3[6+i2] * u2;
		}
	    }
	  nng.xx = 1 - std::sqrt(nng.points(i0,3)*nng.points(i0,3) + 
				 nng.points(i0,4)*nng.points(i0,4) + 
				 nng.points(i0,5)*nng.points(i0,5)) / nng.xx;
	  nng.points(i0,3) /= nng.points(i0,5);
	  nng.points(i0,4) /= nng.points(i0,5);
	  nng.points(i0,5) = nng.xx; 
	}
    }
  for (int i0=0; i0<3; i0++)
    {
      nng.points(nng.datcnt+i0,3) = -nng.bbb;
      nng.points(nng.datcnt+i0,4) = -nng.ccc;
      nng.points(nng.datcnt+i0,5) = T(1);
    }
}


template<class T> void FindNeigh(nngrid<T>& nng,
				 const int& ipt)
{
  if ( ! (nng.rootsimp)->nextsimp ) 
    {
      (nng.rootsimp)->nextsimp = new NNG_simp<T>; // IMakeSimp();
    }
  nng.cursimp = (nng.rootsimp)->nextsimp;
  (nng.cursimp)->vert[0] = nng.datcnt;
  (nng.cursimp)->vert[1] = nng.datcnt + 1;
  (nng.cursimp)->vert[2] = nng.datcnt + 2;
  (nng.cursimp)->cent[0] = (nng.cursimp)->cent[1] = 0.5;
  (nng.cursimp)->cent[2] = nng.bignum;
  nng.numtri = 1;
  nng.lasttemp = nng.roottemp;
  for (int i2=0; i2<3; i2++)
    {
      int j1 = 0;
      if ( j1 == i2 )
	{
	  j1++;
	}
      int j2 = j1 + 1;
      if ( j2 == i2 )
	{
	  j2++;
	}
      if ( ! (nng.lasttemp)->nexttemp ) 
	{
	  (nng.lasttemp)->nexttemp = new NNG_temp; // IMakeTemp();
	}
      nng.lasttemp = (nng.lasttemp)->nexttemp;
      (nng.lasttemp)->end[0] = (nng.cursimp)->vert[j1];
      (nng.lasttemp)->end[1] = (nng.cursimp)->vert[j2];
    }
  nng.curtemp = nng.roottemp;
  for (int i1=0; i1<3; i1++)
    {
      nng.curtemp = (nng.curtemp)->nexttemp;
      for (int i2=0; i2<2; i2++)
	{
	  nng.work3[3*i2] = nng.points((nng.curtemp)->end[i2],0) - nng.points(ipt,0);
	  nng.work3[3*i2+1] = nng.points((nng.curtemp)->end[i2],1) - nng.points(ipt,1);
	  nng.work3[3*i2+2] =
	    nng.work3[3*i2] * (nng.points((nng.curtemp)->end[i2],0) + nng.points(ipt,0)) / 2 +
	    nng.work3[3*i2+1] * (nng.points((nng.curtemp)->end[i2],1) + nng.points(ipt,1)) / 2;
	}
      nng.xx = nng.work3[0] * nng.work3[4] - nng.work3[3] * nng.work3[1];
      (nng.cursimp)->cent[0] = (nng.work3[2] * nng.work3[4] - nng.work3[5] * nng.work3[1]) / nng.xx;
      (nng.cursimp)->cent[1] = (nng.work3[0] * nng.work3[5] - nng.work3[3] * nng.work3[2]) / nng.xx;
      (nng.cursimp)->cent[2] =
	(nng.points(ipt,0) - (nng.cursimp)->cent[0]) *
	(nng.points(ipt,0) - (nng.cursimp)->cent[0]) +
	(nng.points(ipt,1) - (nng.cursimp)->cent[1]) *
	(nng.points(ipt,1) - (nng.cursimp)->cent[1]);
      (nng.cursimp)->vert[0] = (nng.curtemp)->end[0];
      (nng.cursimp)->vert[1] = (nng.curtemp)->end[1];
      (nng.cursimp)->vert[2] = ipt;
      nng.lastsimp = nng.cursimp;
      if ( ! (nng.cursimp)->nextsimp ) 
	{
	  (nng.cursimp)->nextsimp = new NNG_simp<T>; // IMakeSimp();
	}
      nng.cursimp = (nng.cursimp)->nextsimp; 
    }
  nng.numtri += 2;
  for (int i0=0; i0<nng.datcnt; i0++)
    {
      if ( i0 != ipt )
	{
	  int j4 = 0;
	  int j3 = -1;
	  nng.lasttemp = nng.roottemp;
	  nng.cursimp = nng.rootsimp;
	  for (int i1=0; i1<nng.numtri; i1++)
	    {
	      nng.prevsimp = nng.cursimp;
	      nng.cursimp = (nng.cursimp)->nextsimp;
	      nng.xx = (nng.cursimp)->cent[2] - 
		(nng.points(i0,0) - (nng.cursimp)->cent[0]) *
		(nng.points(i0,0) - (nng.cursimp)->cent[0]);
	      if ( nng.xx > T(0) )
		{
		  nng.xx -=
		    (nng.points(i0,1) - (nng.cursimp)->cent[1]) *
		    (nng.points(i0,1) - (nng.cursimp)->cent[1]);
		  if ( nng.xx > T(0) )
		    {
		      j4--;
		      for (int i2=0; i2<3; i2++)
			{
			  bool NextOne = false;
			  int j1 = 0; 
			  if ( j1 == i2 )
			    {
			      j1++; 
			    }
			  int j2 = j1 + 1;    
			  if ( j2 == i2 )
			    {
			      j2++;
			    }
			  if ( j3 > 1 )
			    {
			      int j5 = j3;
			      nng.curtemp = nng.roottemp;
			      for (int i3=0; i3<=j5; i3++)
				{
				  nng.prevtemp = nng.curtemp;
				  nng.curtemp = (nng.curtemp)->nexttemp;
				  if ( (nng.cursimp)->vert[j1] == (nng.curtemp)->end[0] )
				    {
				      if ( (nng.cursimp)->vert[j2] == (nng.curtemp)->end[1] )
					{
					  if ( nng.curtemp == nng.lasttemp ) 
					    {
					      nng.lasttemp = nng.prevtemp;
					    }
					  else
					    {
					      (nng.prevtemp)->nexttemp = (nng.curtemp)->nexttemp;
					      (nng.curtemp)->nexttemp = (nng.lasttemp)->nexttemp;
					      (nng.lasttemp)->nexttemp = nng.curtemp;
					    }
					  j3--;
					  NextOne = true;
					  break;
					}
				    }
				}
			      if ( NextOne )
				{
				  continue;
				}
			    }
			  if ( ! (nng.lasttemp)->nexttemp ) 
			    {
			      (nng.lasttemp)->nexttemp = new NNG_temp; // IMakeTemp();
			    }
			  nng.lasttemp = (nng.lasttemp)->nexttemp;
			  j3++;
			  (nng.lasttemp)->end[0] = (nng.cursimp)->vert[j1];
			  (nng.lasttemp)->end[1] = (nng.cursimp)->vert[j2];
			}
		      if ( nng.cursimp == nng.lastsimp ) 
			{
			  nng.lastsimp = nng.prevsimp;
			}
		      else
			{
			  (nng.prevsimp)->nextsimp = (nng.cursimp)->nextsimp;
			  (nng.cursimp)->nextsimp = (nng.lastsimp)->nextsimp;
			  (nng.lastsimp)->nextsimp = nng.cursimp;
			  nng.cursimp = nng.prevsimp;
			}
		    }
		}
	    }
	  if (j3 > -1)
	    {
	      nng.curtemp = nng.roottemp;
	      nng.cursimp = (nng.lastsimp)->nextsimp;
	      for (int i1=0; i1<=j3; i1++)
		{
		  nng.curtemp = (nng.curtemp)->nexttemp;
		  if ( (nng.curtemp)->end[0] == ipt || (nng.curtemp)->end[1] == ipt )
		    {
		      for (int i2=0; i2<2; i2++)
			{
			  nng.work3[3*i2] = nng.points((nng.curtemp)->end[i2],0) - nng.points(i0,0);
			  nng.work3[3*i2+1] = nng.points((nng.curtemp)->end[i2],1) - nng.points(i0,1);
			  nng.work3[3*i2+2] =
			    nng.work3[3*i2] * (nng.points((nng.curtemp)->end[i2],0) + nng.points(i0,0)) / 2 +
			    nng.work3[3*i2+1] * (nng.points((nng.curtemp)->end[i2],1) + nng.points(i0,1)) / 2;
			}
		      nng.xx = nng.work3[0] * nng.work3[4] - nng.work3[3] * nng.work3[1];
		      (nng.cursimp)->cent[0] = (nng.work3[2] * nng.work3[4] - nng.work3[5] * nng.work3[1]) / nng.xx;
		      (nng.cursimp)->cent[1] = (nng.work3[0] * nng.work3[5] - nng.work3[3] * nng.work3[2]) / nng.xx;
		      (nng.cursimp)->cent[2] = 
			(nng.points(i0,0) - (nng.cursimp)->cent[0]) * (nng.points(i0,0) - (nng.cursimp)->cent[0]) +
			(nng.points(i0,1) - (nng.cursimp)->cent[1]) * (nng.points(i0,1) - (nng.cursimp)->cent[1]);
		      (nng.cursimp)->vert[0] = (nng.curtemp)->end[0];
		      (nng.cursimp)->vert[1] = (nng.curtemp)->end[1];
		      (nng.cursimp)->vert[2] = i0;
		      nng.lastsimp = nng.cursimp;
		      if ( ! (nng.cursimp)->nextsimp ) 
			{
			  (nng.cursimp)->nextsimp = new NNG_simp<T>; // IMakeSimp();
			}
		      nng.cursimp = (nng.cursimp)->nextsimp; 
		      j4++;
		    }
		}
	      nng.numtri += j4;
	    }
	}
    }
  std::fill_n(&((nng.jndx)[0]),nng.datcnt,0);
  nng.cursimp = nng.rootsimp;
  nng.ext=false;
  for (int i1=0; i1<nng.numtri; i1++)
    {
      nng.cursimp = (nng.cursimp)->nextsimp;
      for (int i2=0; i2<3; i2++) 
	{
	  if ((nng.cursimp)->vert[i2] < nng.datcnt)
	    {
	      if ((nng.cursimp)->vert[i2] != ipt) 
		{
		  (nng.jndx)[(nng.cursimp)->vert[i2]] = 1;
		}
	    }
	  else
	    {
	      nng.ext = true; 
	    }
	}
    }
}


template<class T> void TriCentr(nngrid<T>& nng)
{
  if ( ! (nng.rootsimp)->nextsimp ) 
    {
      (nng.rootsimp)->nextsimp = new NNG_simp<T>; // IMakeSimp();
    }
  nng.lastsimp = nng.cursimp = (nng.rootsimp)->nextsimp;
  (nng.cursimp)->vert[0] = nng.datcnt;
  (nng.cursimp)->vert[1] = nng.datcnt + 1;
  (nng.cursimp)->vert[2] = nng.datcnt + 2;
  (nng.cursimp)->cent[0] = (nng.cursimp)->cent[1] = 0.5;
  (nng.cursimp)->cent[2] = nng.bignum;
  nng.numtri = 1;
  for (int i0=0; i0<nng.neicnt; i0++)
    {
      int j3 = -1;
      int j4 = 0;
      nng.lasttemp = nng.roottemp;
      nng.cursimp = nng.rootsimp;
      for (int i1=0; i1<nng.numtri; i1++)
      {
	nng.prevsimp = nng.cursimp;
	nng.cursimp = (nng.cursimp)->nextsimp;
	nng.xx = (nng.cursimp)->cent[2] -
	  (nng.joints(i0,0) - (nng.cursimp)->cent[0]) *
	  (nng.joints(i0,0) - (nng.cursimp)->cent[0]);
	if ( nng.xx > T(0) )
	  {
	    nng.xx -=
	      (nng.joints(i0,1) - (nng.cursimp)->cent[1]) *
	      (nng.joints(i0,1) - (nng.cursimp)->cent[1]);
            if ( nng.xx > T(0) )
	      {
		j4--;
		for (int i2=0; i2<3; i2++)
		  {
		    bool NextOne = false;
		    int j1 = 0;
		    if (j1 == i2)
		      {
			j1++;
		      }
		    int j2 = j1 + 1;
		    if ( j2 == i2 )
		      {
			j2++;
		      }
		    if ( j3 > 1 )
		      {
			int j5 = j3;
			nng.curtemp = nng.roottemp;
			for (int i3=0; i3<=j5; i3++)
			  {
			    nng.prevtemp = nng.curtemp;
			    nng.curtemp = (nng.curtemp)->nexttemp;
			    if ( (nng.cursimp)->vert[j1] == (nng.curtemp)->end[0] )
			      {
				if ( (nng.cursimp)->vert[j2] == (nng.curtemp)->end[1] )
				  {
				    if ( nng.curtemp == nng.lasttemp ) 
				      {
					nng.lasttemp = nng.prevtemp;
				      }
				    else
				      {
					(nng.prevtemp)->nexttemp = (nng.curtemp)->nexttemp;
					(nng.curtemp)->nexttemp = (nng.lasttemp)->nexttemp;
					(nng.lasttemp)->nexttemp = nng.curtemp;
				      }
				    j3--;
				    NextOne = true;
				    break;
				  }
			      }
			  }
			if ( NextOne )
			  {
			    continue;
			  }
		      }
		    if ( ! (nng.lasttemp)->nexttemp ) 
		      {
			(nng.lasttemp)->nexttemp = new NNG_temp; // IMakeTemp();
		      }
		    nng.lasttemp = (nng.lasttemp)->nexttemp;
		    j3++;
		    (nng.lasttemp)->end[0] = (nng.cursimp)->vert[j1];
		    (nng.lasttemp)->end[1] = (nng.cursimp)->vert[j2];
		  }
		if ( nng.cursimp == nng.lastsimp ) 
                  {
		    nng.lastsimp = nng.prevsimp;
		  }
		else
		  {
		    (nng.prevsimp)->nextsimp = (nng.cursimp)->nextsimp;
		    (nng.cursimp)->nextsimp = (nng.lastsimp)->nextsimp;
		    (nng.lastsimp)->nextsimp = nng.cursimp;
		    nng.cursimp = nng.prevsimp;
		  }
	      }
	  }
      }
      nng.curtemp = nng.roottemp;
      nng.cursimp = (nng.lastsimp)->nextsimp;
      for (int i1=0; i1<=j3; i1++)
	{
	  nng.curtemp = (nng.curtemp)->nexttemp;
	  for (int i2=0; i2<2; i2++)
	    {
	      nng.work3[3*i2] = nng.joints((nng.curtemp)->end[i2],0) - nng.joints(i0,0);
	      nng.work3[3*i2+1] = nng.joints((nng.curtemp)->end[i2],1) - nng.joints(i0,1);
	      nng.work3[3*i2+2] =
		nng.work3[3*i2] * (nng.joints((nng.curtemp)->end[i2],0) + nng.joints(i0,0)) / 2 +
		nng.work3[3*i2+1] * (nng.joints((nng.curtemp)->end[i2],1) + nng.joints(i0,1)) / 2;
	    }
	  nng.xx = nng.work3[0] * nng.work3[4] - nng.work3[3] * nng.work3[1];
	  (nng.cursimp)->cent[0] = (nng.work3[2] * nng.work3[4] - nng.work3[5] * nng.work3[1]) / nng.xx;
	  (nng.cursimp)->cent[1] = (nng.work3[0] * nng.work3[5] - nng.work3[3] * nng.work3[2]) / nng.xx;
	  (nng.cursimp)->cent[2] = 
            (nng.joints(i0,0) - (nng.cursimp)->cent[0]) * (nng.joints(i0,0) - (nng.cursimp)->cent[0]) +
            (nng.joints(i0,1) - (nng.cursimp)->cent[1]) * (nng.joints(i0,1) - (nng.cursimp)->cent[1]);
	  (nng.cursimp)->vert[0] = (nng.curtemp)->end[0];
	  (nng.cursimp)->vert[1] = (nng.curtemp)->end[1];
	  (nng.cursimp)->vert[2] = i0;
	  nng.lastsimp = nng.cursimp;
	  if ( ! (nng.cursimp)->nextsimp ) 
            {
	      (nng.cursimp)->nextsimp = new NNG_simp<T>; // IMakeSimp();
	    }
	  nng.cursimp = (nng.cursimp)->nextsimp;   
	  j4++;
	}
      nng.numtri += j4;
    }
}


template<class T> void TriNeigh(nngrid<T>& nng)
{
  if ( ! (nng.rootsimp)->nextsimp  )
    {
      (nng.rootsimp)->nextsimp = new NNG_simp<T>; // IMakeSimp();
    }
  nng.lastsimp = nng.cursimp = (nng.rootsimp)->nextsimp;
  (nng.cursimp)->vert[0] = nng.datcnt;
  (nng.cursimp)->vert[1] = nng.datcnt + 1;
  (nng.cursimp)->vert[2] = nng.datcnt + 2;
  (nng.cursimp)->cent[0] = (nng.cursimp)->cent[1] = 0.5;
  (nng.cursimp)->cent[2] = nng.bignum;
  nng.numtri = 1;
  for (int i0=0; i0<nng.datcnt; i0++)
    {
      if ( (nng.jndx)[i0] )
	{
	  int j3 = -1;
	  nng.lasttemp = nng.roottemp;
	  nng.cursimp = nng.rootsimp;
	  for (int i1=0; i1<nng.numtri; i1++)
	    {
	      nng.prevsimp = nng.cursimp;
	      nng.cursimp = (nng.cursimp)->nextsimp;
	      nng.xx = (nng.cursimp)->cent[2] -
		(nng.points(i0,0) - (nng.cursimp)->cent[0]) *
		(nng.points(i0,0) - (nng.cursimp)->cent[0]);
	      if ( nng.xx > T(0) )
		{
		  nng.xx -=
		    (nng.points(i0,1) - (nng.cursimp)->cent[1]) *
		    (nng.points(i0,1) - (nng.cursimp)->cent[1]);
		  if ( nng.xx > T(0) )
		    {
		      for (int i2=0; i2<3; i2++)
			{
			  bool NextOne = false;
			  int j1 = 0;
			  if ( j1 == i2 )
			    {
			      j1++;
			    }
			  int j2 = j1 + 1;
			  if ( j2 == i2 )
			    {
			      j2++;
			    }
			  if ( j3 > 1 )
			    {
			      int j5 = j3;
			      nng.curtemp = nng.roottemp;
			      for (int i3=0; i3<=j5; i3++)
				{
				  nng.prevtemp = nng.curtemp;
				  nng.curtemp = (nng.curtemp)->nexttemp;
				  if ( (nng.cursimp)->vert[j1] ==
				       (nng.curtemp)->end[0] )
				    {
				      if ( (nng.cursimp)->vert[j2] ==
					   (nng.curtemp)->end[1] )
					{
					  if ( nng.curtemp == nng.lasttemp )
					    {
					      nng.lasttemp = nng.prevtemp;
					    }
					  else
					    {
					      (nng.prevtemp)->nexttemp =
						(nng.curtemp)->nexttemp;
					      (nng.curtemp)->nexttemp =
						(nng.lasttemp)->nexttemp;
					      (nng.lasttemp)->nexttemp =
						nng.curtemp;
					    }
					  j3--;
					  NextOne = true;
					  break;
					}
				    }
				}
			      if ( NextOne )
				{
				  continue;
				}
			    }
			  if ( ! (nng.lasttemp)->nexttemp )
			    {
			      (nng.lasttemp)->nexttemp = new NNG_temp; // IMakeTemp();
			    }
			  nng.lasttemp = (nng.lasttemp)->nexttemp;
			  j3++;
			  (nng.lasttemp)->end[0] = (nng.cursimp)->vert[j1];
			  (nng.lasttemp)->end[1] = (nng.cursimp)->vert[j2];
			}
		      if ( nng.cursimp == nng.lastsimp )
			{ 
			  nng.lastsimp = nng.prevsimp;
			}
		      else
			{
			  (nng.prevsimp)->nextsimp = (nng.cursimp)->nextsimp;
			  (nng.cursimp)->nextsimp = (nng.lastsimp)->nextsimp;
			  (nng.lastsimp)->nextsimp = nng.cursimp;
			  nng.cursimp = nng.prevsimp;
			}
		    }
		}
	    }
	  nng.curtemp = nng.roottemp;
	  nng.cursimp = (nng.lastsimp)->nextsimp;
	  for (int i1=0; i1<=j3; i1++)
	    {
	      nng.curtemp = (nng.curtemp)->nexttemp;
	      for (int i2=0; i2<2; i2++)
		{
		  nng.work3[3*i2] = nng.points((nng.curtemp)->end[i2],0) - nng.points(i0,0);
		  nng.work3[3*i2+1] = nng.points((nng.curtemp)->end[i2],1) - nng.points(i0,1);
		  nng.work3[3*i2+2] =
		    nng.work3[3*i2] * (nng.points((nng.curtemp)->end[i2],0) + nng.points(i0,0)) / 2 +
		    nng.work3[3*i2+1] * (nng.points((nng.curtemp)->end[i2],1) + nng.points(i0,1)) / 2;
		}
	      nng.xx = nng.work3[0] * nng.work3[4] - nng.work3[3] * nng.work3[1];
	      (nng.cursimp)->cent[0] = (nng.work3[2] * nng.work3[4] - nng.work3[5] * nng.work3[1]) /nng.xx;
	      (nng.cursimp)->cent[1] = (nng.work3[0] * nng.work3[5] - nng.work3[3] * nng.work3[2]) /nng.xx;
	      (nng.cursimp)->cent[2] =
		(nng.points(i0,0) - (nng.cursimp)->cent[0]) * (nng.points(i0,0) - (nng.cursimp)->cent[0]) +
		(nng.points(i0,1) - (nng.cursimp)->cent[1]) * (nng.points(i0,1) - (nng.cursimp)->cent[1]);
	      (nng.cursimp)->vert[0] = (nng.curtemp)->end[0];
	      (nng.cursimp)->vert[1] = (nng.curtemp)->end[1];
	      (nng.cursimp)->vert[2] = i0;
	      nng.lastsimp = nng.cursimp;
	      if ( ! (nng.cursimp)->nextsimp ) 
		{
		  (nng.cursimp)->nextsimp = new NNG_simp<T>; // IMakeSimp();
		}
	      nng.cursimp = (nng.cursimp)->nextsimp;
	    }
	  nng.numtri += 2;
	}
    }
  nng.cursimp = nng.rootsimp;
  nng.asum = T(0);
  for (int i0=0; i0<nng.numtri; i0++) 
    {
      nng.cursimp = (nng.cursimp)->nextsimp;
      for (int i1=0; i1<2; i1++)
	{ 
	  nng.work3[i1] = nng.points((nng.cursimp)->vert[1],i1) - nng.points((nng.cursimp)->vert[0],i1);
	  nng.work3[3+i1] = nng.points((nng.cursimp)->vert[2],i1) - nng.points((nng.cursimp)->vert[0],i1);
	}
      nng.xx = nng.work3[0] * nng.work3[4] - nng.work3[1] * nng.work3[3];
      if ( nng.xx < T(0) )
	{
	  int j4 = (nng.cursimp)->vert[2];
	  (nng.cursimp)->vert[2] = (nng.cursimp)->vert[1];
	  (nng.cursimp)->vert[1] = j4;
	  if ( (nng.cursimp)->vert[0] < nng.datcnt ) 
            {
	      nng.asum -= nng.xx / 2;
	    }
	}
      else if ( (nng.cursimp)->vert[0] < nng.datcnt ) 
	{
	  nng.asum += nng.xx / 2;
	}
    }
}


template<class T> bool MakeGrid(nngrid<T>& nng,
				const nngridr_ini<T>& ini,
				std::vector< std::valarray<T> >& upoints,
				std::valarray<T>& udata)
{
  upoints.resize(ini.x_nodes * ini.y_nodes);
  for (size_t i=0; i<upoints.size(); i++)
    {
      upoints[i].resize(2);
    }
  udata.resize(ini.x_nodes * ini.y_nodes);
  size_t pdata = 0;
  std::fill_n(&((nng.jndx)[0]), nng.datcnt, 1);
  TriNeigh(nng);
  T x_increm = ( 0.00000001 + ini.xterm - ini.xstart ) / ( ini.x_nodes - 1 );
  T y_increm = ( 0.00000001 + ini.yterm - ini.ystart ) / ( ini.y_nodes - 1 );
  T wyd = ini.ystart - y_increm;
  if ( ini.arriba < 0 )
    {
      wyd = ini.yterm + y_increm;
    }
  T surf, surfe, surfn;
  for (int i=0; i<ini.y_nodes; i++) 
    {
      wyd += y_increm * ini.arriba;
      nng.points(nng.datcnt3,1) = wyd;
      T wxd = ini.xstart - x_increm;
      for (int j=0; j<ini.x_nodes; j++) 
	{
	  wxd += x_increm;
	  nng.points(nng.datcnt3,0) = wxd;
	  FindProp(nng,wxd,wyd);
	  if ( !ini.extrap && !nng.goodflag )
	    {
	      surf = ini.nuldat;
	    }
	  else
	    {
	      surf = Surface(nng);
	      if ( ini.igrad )
		{
		  surf = Meld(nng, ini, surf, wxd, wyd);
		}
	      if ( ini.non_neg )
		{
		  if ( surf < T(0) )
		    {
		      surf = T(0);
		    }
		}
	    }
	  // Here's the useful grid outputs (RF - May 2003)
	  upoints[i*ini.x_nodes+j][0] = wxd;
	  upoints[i*ini.x_nodes+j][1] = wyd;
	  udata[pdata++] = surf;
	}
    }
  return true;
}


template<class T> void FindProp(nngrid<T>& nng,
				const T& wxd,
				const T& wyd)
{
  const int scor[3][2] = {{1,2},{2,0},{0,1}};
  T work3[3][3], work4[3][2];
  nng.lastneig = nng.rootneig;
  nng.goodflag = false;
  nng.numnei = -1;
  nng.cursimp = nng.rootsimp;
  for (int i2=0; i2<nng.numtri; i2++)
    {
      nng.cursimp = (nng.cursimp)->nextsimp;
      T xx = (nng.cursimp)->cent[2] -
	(wxd - (nng.cursimp)->cent[0]) * (wxd - (nng.cursimp)->cent[0]);
      if ( xx > T(0) )
	{
	  xx -= (wyd - (nng.cursimp)->cent[1])*(wyd - (nng.cursimp)->cent[1]);
	  if ( xx > T(0) )
	    {
	      bool inside = false;
	      if ( (nng.cursimp)->vert[0] < nng.datcnt )
		{
		  inside = true;
		}
	      for (int i3=0; i3<3; i3++)
		{
		  for (int i4=0; i4<2; i4++)
		    {
		      work3[i4][0] = nng.points((nng.cursimp)->vert[scor[i3][i4]],0) - wxd;
		      work3[i4][1] = nng.points((nng.cursimp)->vert[scor[i3][i4]],1) - wyd;
		      work3[i4][2] = work3[i4][0] *
			(nng.points((nng.cursimp)->vert[scor[i3][i4]],0) + wxd) / 2 + work3[i4][1] *
			(nng.points((nng.cursimp)->vert[scor[i3][i4]],1) + wyd) / 2;
		    }
		  xx =  work3[0][0] * work3[1][1] - work3[1][0] * work3[0][1];
		  work4[i3][0] = ( work3[0][2] * work3[1][1] - work3[1][2] * work3[0][1] ) / xx;
		  work4[i3][1] = ( work3[0][0] * work3[1][2] - work3[1][0] * work3[0][2] ) / xx;
		}
	      int pos_count = 0;
	      for (int i3=0; i3<3; i3++)
		{
		  work3[2][i3] = 
		    ( (work4[scor[i3][0]][0] - (nng.cursimp)->cent[0]) * (work4[scor[i3][1]][1] - (nng.cursimp)->cent[1]) -
		      (work4[scor[i3][1]][0] - (nng.cursimp)->cent[0]) * (work4[scor[i3][0]][1] - (nng.cursimp)->cent[1]) ) / 2;
		  if ( work3[2][i3] > T(0) )
		    {
		      pos_count++;
		    }
		}
	      if ( pos_count > 2 && inside )
		{
		  nng.goodflag = true;
		}
	      for (int i3=0; i3<3; i3++)
		{
		  bool gotem = false;
		  if ( nng.numnei > 1 )
		    {
		      nng.curneig = nng.rootneig;
		      for (int i4=0; i4<=nng.numnei; i4++)
			{
			  nng.curneig = (nng.curneig)->nextneig;
			  if ( (nng.cursimp)->vert[i3] == (nng.curneig)->neinum )
			    {
			      (nng.curneig)->narea += work3[2][i3];
			      gotem = true;
			      break;
			    }
			}
		      if ( gotem )
			{
			  continue;
			}
		    }
		  if ( ! (nng.lastneig)->nextneig )
		    {
		      (nng.lastneig)->nextneig = new NNG_neig<T>; // IMakeNeig();
		    }
		  nng.lastneig = (nng.lastneig)->nextneig;
		  nng.numnei++;
		  (nng.lastneig)->neinum = (nng.cursimp)->vert[i3];
		  (nng.lastneig)->narea = work3[2][i3];
		}   
	    }
	}
    }
}


template<class T> T Surface(nngrid<T>& nng)
{
  T xx = T(0);
  nng.curneig = nng.rootneig;
  for (int i0=0; i0<=nng.numnei; i0++) 
    {
      nng.curneig = (nng.curneig)->nextneig;
      xx += (nng.curneig)->narea;
    }
  T asurf = T(0);
  nng.curneig = nng.rootneig;
  for (int i0=0; i0<=nng.numnei; i0++)
    {
      nng.curneig = (nng.curneig)->nextneig;
      (nng.curneig)->narea /= xx;
      asurf += (nng.curneig)->narea * nng.points((nng.curneig)->neinum,2);
    }
  return asurf;
}


template<class T> T Meld(nngrid<T>& nng,
			 const nngridr_ini<T>& ini,
			 T asurf,
			 const T& wxd,
			 const T& wyd)
{
  nng.curneig = nng.rootneig;
  for (int i0=0; i0<=nng.numnei; i0++)
    {
      nng.curneig = (nng.curneig)->nextneig;
      (nng.curneig)->coord = 0;
      if ( (nng.curneig)->narea > T(0.00001) && (nng.curneig)->narea < T(2) )
	{
	  if ( std::abs(nng.points((nng.curneig)->neinum,5)) > T(0.00001) )
	    {
	      T rS = std::abs(nng.points((nng.curneig)->neinum,5)) + ini.bI;
	      T rT = rS * ini.bJ;
	      T rB = 1 / rT;
	      T bD = std::pow((nng.curneig)->narea, rT);
	      T bB = bD * 2;
	      if ( bD > T(0.5) )
		{
		  bB = (1 - bD) * 2;
		}
	      bB = std::pow(bB, rS) / 2;
	      if ( bD > T(0.5) )
		{
		  bB = 1 - bB;
		}
	      T hP = std::pow(bB, rB);
	      (nng.curneig)->coord =
		( (nng.points((nng.curneig)->neinum,3) * nng.points((nng.curneig)->neinum,0) + 
		   nng.points((nng.curneig)->neinum,4) * nng.points((nng.curneig)->neinum,1) + 
		   nng.points((nng.curneig)->neinum,2) - nng.points((nng.curneig)->neinum,3) * wxd -
		   nng.points((nng.curneig)->neinum,4) * wyd) - asurf) * hP;
	    }
	}
    }
  nng.curneig = nng.rootneig;
  for (int i0=0; i0<=nng.numnei; i0++) 
    {
      nng.curneig = (nng.curneig)->nextneig;
      asurf += (nng.curneig)->coord;
    }
  return asurf; 
}


#endif /* _NNGRIDR_ */
