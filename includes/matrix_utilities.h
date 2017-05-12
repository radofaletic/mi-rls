/*
  matrix_utilities
  
  Rado Faletic
  Department of Physics
  Faculty of Science
  Australian National University  ACT  0200
  Australia
  
  Rado.Faletic@anu.edu.au
  5th August 2004
*/


#ifndef _MATRIX_UTILITIES_
#define _MATRIX_UTILITIES_


/* ---------- header files ---------- */
#include <valarray>
#include <vector>
/* ---------- user header files ---------- */
#include "front-end.h"
#include "grid.h"
#include "sparse_matrix.h"
/* ---------------------------------- */


/* -------------------------------------------------- */
/* ---------- class & function declaration ---------- */
/* -------------------------------------------------- */


template<class T> void AddDamping(SparseMatrix<T>&,
				  std::valarray<T>&,
				  const grid<T>&,
				  const T& = T(1));

template<class T> void AddSmoothing(SparseMatrix<T>&,
				    std::valarray<T>&,
				    const grid<T>&,
				    const T& = T(1));


/* ------------------------------------------ */
/* ---------- function definitions ---------- */
/* ------------------------------------------ */


/* ---------- AddDamping ---------- */
/* add damping to the equation Ax=b */
template<class T> void
AddDamping(SparseMatrix<T>& A, std::valarray<T>& b, const grid<T>& g, const T& damp)
{
  debug("AddDamping","adding damping to matrix system...");

  const std::valarray<bool>* refs = A.Referenced();

  if ( A.AddDiagonal(g.ncells(),damp*g.scale()) )
    {
      std::valarray<T> bv = b;
      size_t sor = 0;
      for (size_t i=0; i<refs->size(); i++)
	{
	  if ( (*refs)[i] )
	    {
	      sor++;
	    }
	}
      b.resize(bv.size()+sor,T(0));
      b[std::slice(0, bv.size(), 1)] = bv;
      debug(" - added a "+ntos(sor)+"*"+ntos(g.ncells())+" submatrix with diagonal "+ntos(damp*g.scale()));
    }
  else
    {
      debug(" - no damping to add");
    }
}
/* -------------------------------- */

/* ---------- AddSmoothing ---------- */
/* add smoothing to the equation Ax=b */
template<class T> void
AddSmoothing(SparseMatrix<T>& A, std::valarray<T>& b, const grid<T>& g, const T& smoothing)
{
  debug("AddSmoothing","adding smoothing to matrix system...");

  const std::valarray<bool>* refs = A.Referenced();

  std::vector< std::valarray<size_t> > s_neigh(0);
  std::vector<size_t> s_i(0);
  for (size_t i=0; i<g.ncells(); i++)
    {
      if ( !((*refs)[i]) )
	{
	  continue;
	}
      std::valarray<size_t> gn = g.get_neighbours(i);
      std::valarray<bool> gnref(true,gn.size());
      for (size_t j=0; j<gn.size(); j++)
	{
	  if ( !((*refs)[gn[j]]) )
	    {
	      gnref[j] = false;
	    }
	}
      s_neigh.push_back(gn[gnref]);
      s_i.push_back(i);
    }

  if ( A.AddNReverse(s_neigh, s_i, smoothing) )
    {
      std::valarray<T> bv = b;
      size_t sor = 0;
      for (size_t i=0; i<refs->size(); i++)
	{
	  if ( (*refs)[i] )
	    {
	      sor++;
	    }
	}
      b.resize(bv.size()+sor,T(0));
      b[std::slice(0, bv.size(), 1)] = bv;
      debug(" - added "+ntos(sor)+" smoothing rows, with smoothing "+ntos(smoothing));
    }
  else
    {
      debug(" - no smoothing to add");
    }
}
/* ---------------------------------- */


#endif /* _MATRIX_UTILITIES_ */
