/**
 lsqr
 This is a C version of LSQR, derived from the Fortran 77 implementation of C. C. Paige and M. A. Saunders.
 
 Contains auxiliary functions, data type definitions, and function prototypes for the iterative linear solver LSQR. This file defines functions for data type allocation and deallocation, lsqr itself (the main algorithm), and functions that scale, copy, and compute the Euclidean norm of a vector.
 
 08 Sep 1999: First version from James W. Howse <jhowse@lanl.gov>
 07 Nov 2004: updated for C++ by Rado Faletic <Rado.Faletic@anu.edu.au>
 18 Apr 2022: updated by Rado Faletic
 */





#ifndef _LSQR_
#define _LSQR_





/*
 *------------------------------------------------------------------------------
 *
 *     LSQR  finds a solution x to the following problems:
 *
 *     1. Unsymmetric equations --    solve  A*x = b
 *
 *     2. Linear least squares  --    solve  A*x = b
 *                                    in the least-squares sense
 *
 *     3. Damped least squares  --    solve  (   A    )*x = ( b )
 *                                           ( damp*I )     ( 0 )
 *                                    in the least-squares sense
 *
 *     where 'A' is a matrix with 'm' rows and 'n' columns, 'b' is an
 *     'm'-vector, and 'damp' is a scalar.  (All quantities are real.)
 *     The matrix 'A' is intended to be large and sparse.
 *
 *
 *     Notation
 *     --------
 *
 *     The following quantities are used in discussing the subroutine
 *     parameters:
 *
 *     'Abar'   =  (   A    ),          'bbar'  =  ( b )
 *                 ( damp*I )                      ( 0 )
 *
 *     'r'      =  b  -  A*x,           'rbar'  =  bbar  -  Abar*x
 *
 *     'rnorm'  =  sqrt( norm(r)**2  +  damp**2 * norm(x)**2 )
 *              =  norm( rbar )
 *
 *     'rel_prec'  =  the relative precision of floating-point arithmetic
 *                    on the machine being used.  For example, on the IBM 370,
 *                    'rel_prec' is about 1.0E-6 and 1.0D-16 in single and double
 *                    precision respectively.
 *
 *     LSQR  minimizes the function 'rnorm' with respect to 'x'.
 *
 *------------------------------------------------------------------------------
 */





/* ---------- standard header files ---------- */
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>
#include <valarray>
#include <vector>





/* ------------ user header files ------------ */
#include "sparse_matrix.h"





/* ---------- function definitions ---------- */





template<class T> inline T sqr(const T& x) { return x * x; };





/*-------------------------------------------------------------------------*/
/*                                                                         */
/*  Define the function 'nrm2()'.  This function takes a vector arguement  */
/*  and computes the Euclidean or L2 norm of this vector.  Note  that this */
/*  is a version of the BLAS function 'dnrm2()' rewritten to use the       */
/*  current data structures.                                               */
/*                                                                         */
/*-------------------------------------------------------------------------*/
template<class T> inline T nrm2(const std::valarray<T>& vec)
{
    return std::sqrt( ( vec * vec ).sum() );
}





/*
 *------------------------------------------------------------------------------
 *
 *     Input Quantities
 *     ----------------
 *
 *     damp_val     input  The damping parameter for problem 3 above.
 *                         ('damp_val' should be 0.0 for problems 1 and 2.)
 *                         If the system A*x = b is incompatible, values
 *                         of 'damp_val' in the range
 *                            0 to sqrt('rel_prec')*norm(A)
 *                         will probably have a negligible effect.
 *                         Larger values of 'damp_val' will tend to decrease
 *                         the norm of x and reduce the number of
 *                         iterations required by LSQR.
 *
 *                         The work per iteration and the storage needed
 *                         by LSQR are the same for all values of 'damp_val'.
 *
 *     rel_mat_err  input  An estimate of the relative error in the data
 *                         defining the matrix 'A'.  For example,
 *                         if 'A' is accurate to about 6 digits, set
 *                         rel_mat_err = 1.0e-6 .
 *
 *     rel_rhs_err  input  An extimate of the relative error in the data
 *                         defining the right hand side (rhs) vector 'b'.  For
 *                         example, if 'b' is accurate to about 6 digits, set
 *                         rel_rhs_err = 1.0e-6 .
 *
 *     cond_lim     input  An upper limit on cond('Abar'), the apparent
 *                         condition number of the matrix 'Abar'.
 *                         Iterations will be terminated if a computed
 *                         estimate of cond('Abar') exceeds 'cond_lim'.
 *                         This is intended to prevent certain small or
 *                         zero singular values of 'A' or 'Abar' from
 *                         coming into effect and causing unwanted growth
 *                         in the computed solution.
 *
 *                         'cond_lim' and 'damp_val' may be used separately or
 *                         together to regularize ill-conditioned systems.
 *
 *                         Normally, 'cond_lim' should be in the range
 *                         1000 to 1/rel_prec.
 *                         Suggested value:
 *                         cond_lim = 1/(100*rel_prec)  for compatible systems,
 *                         cond_lim = 1/(10*sqrt(rel_prec)) for least squares.
 *
 *             Note:  If the user is not concerned about the parameters
 *             'rel_mat_err', 'rel_rhs_err' and 'cond_lim', any or all of them
 *             may be set to zero.  The effect will be the same as the values
 *             'rel_prec', 'rel_prec' and 1/rel_prec respectively.
 *
 *     max_iter     input  An upper limit on the number of iterations.
 *                         Suggested value:
 *                         max_iter = n/2   for well-conditioned systems
 *                                          with clustered singular values,
 *                         max_iter = 4*n   otherwise.
 *
 *------------------------------------------------------------------------------
 */
template<class T>
class lsqr_input
{
public:
    T                damp_val;
    T                rel_mat_err;
    T                rel_rhs_err;
    T                cond_lim;
    std::size_t      max_iter;
    bool             all_positive;
    lsqr_input(const T& dv = T(0), const T& rme = T(0), const T& rre = T(0), const T& cl = T(0), const std::size_t& mi = 100, const bool& ap = false)
    : damp_val(dv), rel_mat_err(rme), rel_rhs_err(rre), cond_lim(cl), max_iter(mi), all_positive(ap) { };
};





/*
 *------------------------------------------------------------------------------
 *
 *     Output Quantities
 *     -----------------
 *
 *     term_flag       output  An integer giving the reason for termination:
 *
 *                     0       x = x0  is the exact solution.
 *                             No iterations were performed.
 *
 *                     1       The equations A*x = b are probably compatible.
 *                             Norm(A*x - b) is sufficiently small, given the
 *                             values of 'rel_mat_err' and 'rel_rhs_err'.
 *
 *                     2       The system A*x = b is probably not
 *                             compatible.  A least-squares solution has
 *                             been obtained that is sufficiently accurate,
 *                             given the value of 'rel_mat_err'.
 *
 *                     3       An estimate of cond('Abar') has exceeded
 *                             'cond_lim'.  The system A*x = b appears to be
 *                             ill-conditioned.  Otherwise, there could be an
 *                             error in subroutine APROD.
 *
 *                     4       The equations A*x = b are probably
 *                             compatible.  Norm(A*x - b) is as small as
 *                             seems reasonable on this machine.
 *
 *                     5       The system A*x = b is probably not
 *                             compatible.  A least-squares solution has
 *                             been obtained that is as accurate as seems
 *                             reasonable on this machine.
 *
 *                     6       Cond('Abar') seems to be so large that there is
 *                             no point in doing further iterations,
 *                             given the precision of this machine.
 *                             There could be an error in subroutine APROD.
 *
 *                     7       The iteration limit 'max_iter' was reached.
 *
 *     num_iters       output  The number of iterations performed.
 *
 *     frob_mat_norm   output  An estimate of the Frobenius norm of 'Abar'.
 *                             This is the square-root of the sum of squares
 *                             of the elements of 'Abar'.
 *                             If 'damp_val' is small and if the columns of 'A'
 *                             have all been scaled to have length 1.0,
 *                             'frob_mat_norm' should increase to roughly
 *                             sqrt('n').
 *                             A radically different value for 'frob_mat_norm'
 *                             may indicate an error in subroutine APROD (there
 *                             may be an inconsistency between modes 1 and 2).
 *
 *     mat_cond_num    output  An estimate of cond('Abar'), the condition
 *                             number of 'Abar'.  A very high value of
 *                             'mat_cond_num'
 *                             may again indicate an error in APROD.
 *
 *     resid_norm      output  An estimate of the final value of norm('rbar'),
 *                             the function being minimized (see notation
 *                             above).  This will be small if A*x = b has
 *                             a solution.
 *
 *     mat_resid_norm  output  An estimate of the final value of
 *                             norm( Abar(transpose)*rbar ), the norm of
 *                             the residual for the usual normal equations.
 *                             This should be small in all cases.
 *                             ('mat_resid_norm'
 *                             will often be smaller than the true value
 *                             computed from the output vector 'x'.)
 *
 *     sol_norm        output  An estimate of the norm of the final
 *                             solution vector 'x'.
 *
 *     std_err_vec     output  The vector which returns the standard error
 *                             estimates  for the components of 'x'.
 *                             This vector has a length of 'num_cols'.
 *                             For each i, std_err_vec(i) is set to the value
 *                             'resid_norm' * sqrt( sigma(i,i) / T ),
 *                             where sigma(i,i) is an estimate of the i-th
 *                             diagonal of the inverse of Abar(transpose)*Abar
 *                             and  T = 1      if  m <= n,
 *                                  T = m - n  if  m > n  and  damp_val = 0,
 *                                  T = m      if  damp_val != 0.
 *
 *------------------------------------------------------------------------------
 */
template<class T> class lsqr_output
{
public:
    unsigned short   term_flag;
    std::size_t      num_iters;
    T                frob_mat_norm;
    T                mat_cond_num;
    T                resid_norm;
    T                mat_resid_norm;
    T                b_norm;
    T                sol_norm;
    std::valarray<T> std_err_vec;
    lsqr_output(const unsigned short& tf = 0, const std::size_t& ni = 0,
                const T& fmn = T(0), const T& mcn = T(0), const T& rn = T(0),
                const T& mrn = T(0), const T& bn = T(0), const T& sn = T(0),
                const std::size_t& sevs = 0)
    : term_flag(tf), num_iters(ni), frob_mat_norm(fmn), mat_cond_num(mcn), resid_norm(rn), mat_resid_norm(mrn), b_norm(bn), sol_norm(sn), std_err_vec(std::valarray<T>(T(0), sevs)) { };
    void write(const std::string& filename)
    {
        std::vector<std::string> term_msg(8);
        term_msg[0] = "The exact solution is x = x0";
        term_msg[1] = "The residual Ax - b is small enough, given ATOL and BTOL";
        term_msg[2] = "The least squares error is small enough, given ATOL";
        term_msg[3] = "The estimated condition number has exceeded CONLIM";
        term_msg[4] = "The residual Ax - b is small enough, given machine precision";
        term_msg[5] = "The least squares error is small enough, given machine precision";
        term_msg[6] = "The estimated condition number has exceeded machine precision";
        term_msg[7] = "The iteration limit has been reached";
        
        std::ofstream lsqrfile(filename.c_str());
        
        lsqrfile.setf(std::ios_base::scientific, std::ios_base::floatfield);
        lsqrfile << "\n\tISTOP = " << std::setw(3) << this->term_flag
        << "\t\t\tITER = " << std::setw(9) << this->num_iters << "\n"
        << "	|| A ||_F = " << std::setw(13) << std::setprecision(5) << this->frob_mat_norm
        << "\tcond( A ) =      " << std::setw(13) << std::setprecision(5) << this->mat_cond_num << "\n"
        << "	|| r ||_2 = " << std::setw(13) << std::setprecision(5) << this->resid_norm
        << "\t|| A^T r ||_2 =  " << std::setw(13) << std::setprecision(5) << this->mat_resid_norm << "\n"
        << "	|| b ||_2 = " << std::setw(13) << std::setprecision(5) << this->b_norm
        << "\t|| x - x0 ||_2 = " << std::setw(13) << std::setprecision(5) << this->sol_norm << "\n"
        << std::endl;
        lsqrfile << "  " << term_msg[this->term_flag] << "\n" << std::endl;
        
        lsqrfile.close();
    }
};





/*-------------------------------------------------------------------------*/
/*                                                                         */
/*  Define the LSQR function.                                              */
/*                                                                         */
/*-------------------------------------------------------------------------*/
template<class T> lsqr_output<T> LSQR(SparseMatrix<T>&, std::valarray<T>&, const std::valarray<T>&, const lsqr_input<T>&);

/*
 *------------------------------------------------------------------------------
 *
 *     LSQR  finds a solution x to the following problems:
 *
 *     1. Unsymmetric equations --    solve  A*x = b
 *
 *     2. Linear least squares  --    solve  A*x = b
 *                                    in the least-squares sense
 *
 *     3. Damped least squares  --    solve  (   A    )*x = ( b )
 *                                           ( damp*I )     ( 0 )
 *                                    in the least-squares sense
 *
 *     where 'A' is a matrix with 'm' rows and 'n' columns, 'b' is an
 *     'm'-vector, and 'damp' is a scalar.  (All quantities are real.)
 *     The matrix 'A' is intended to be large and sparse.
 *
 *
 *     Notation
 *     --------
 *
 *     The following quantities are used in discussing the subroutine
 *     parameters:
 *
 *     'Abar'   =  (   A    ),          'bbar'  =  ( b )
 *                 ( damp*I )                      ( 0 )
 *
 *     'r'      =  b  -  A*x,           'rbar'  =  bbar  -  Abar*x
 *
 *     'rnorm'  =  sqrt( norm(r)**2  +  damp**2 * norm(x)**2 )
 *              =  norm( rbar )
 *
 *     'rel_prec'  =  the relative precision of floating-point arithmetic
 *                    on the machine being used.  Typically 2.22e-16
 *                    with 64-bit arithmetic.
 *
 *     LSQR  minimizes the function 'rnorm' with respect to 'x'.
 *
 *
 *     References
 *     ----------
 *
 *     C.C. Paige and M.A. Saunders,  LSQR: An algorithm for sparse
 *          linear equations and sparse least squares,
 *          ACM Transactions on Mathematical Software 8, 1 (March 1982),
 *          pp. 43-71.
 *
 *     C.C. Paige and M.A. Saunders,  Algorithm 583, LSQR: Sparse
 *          linear equations and least-squares problems,
 *          ACM Transactions on Mathematical Software 8, 2 (June 1982),
 *          pp. 195-209.
 *
 *     C.L. Lawson, R.J. Hanson, D.R. Kincaid and F.T. Krogh,
 *          Basic linear algebra subprograms for Fortran usage,
 *          ACM Transactions on Mathematical Software 5, 3 (Sept 1979),
 *          pp. 308-323 and 324-325.
 *
 *------------------------------------------------------------------------------
 */
template<class T> lsqr_output<T> LSQR(SparseMatrix<T>& A, std::valarray<T>& x, const std::valarray<T>& b, const lsqr_input<T>& input)
{
    std::vector<std::string> term_msg(8);
    term_msg[0] = "The exact solution is x = x0";
    term_msg[1] = "The residual Ax - b is small enough, given ATOL and BTOL";
    term_msg[2] = "The least squares error is small enough, given ATOL";
    term_msg[3] = "The estimated condition number has exceeded CONLIM";
    term_msg[4] = "The residual Ax - b is small enough, given machine precision";
    term_msg[5] = "The least squares error is small enough, given machine precision";
    term_msg[6] = "The estimated condition number has exceeded machine precision";
    term_msg[7] = "The iteration limit has been reached";
    
    std::ostream& lsqr_fp_out = std::cout;
#if DEBUG > 0
    //lsqr_fp_out = std::cout;
#endif /* DEBUG */
    
    if ( lsqr_fp_out )
    {
        lsqr_fp_out.setf(std::ios_base::scientific, std::ios_base::floatfield);
        lsqr_fp_out << "  Least Squares Solution of A*x = b\n"
        << "	The matrix A has " << std::setw(7) << A.num_rows()
        << " rows and " << std::setw(7) << A.num_cols() << " columns\n"
        << "	The damping parameter is\tDAMP = " << std::setw(10) << std::setprecision(2) << input.damp_val << "\n"
        << "	ATOL = " << std::setw(10) << std::setprecision(2) << input.rel_mat_err
        << "\t\tCONDLIM = " << std::setw(10) << std::setprecision(2) << input.cond_lim << "\n"
        << "	BTOL = " << std::setw(10) << std::setprecision(2) << input.rel_rhs_err
        << "\t\tITERLIM = " << std::setw(10) << input.max_iter << "\n"
        << std::endl;
    }
    
    lsqr_output<T> output(0, 0, T(0), T(0), T(0), T(0), T(0), A.num_cols());
    
    long term_iter = 0;
    
    /*
     *     bidiag_wrk_vec  workspace  This float vector is a workspace for the
     *                                current iteration of the
     *                                Lanczos bidiagonalization.
     *                                This vector has length 'num_cols'.
     *
     *     srch_dir_vec    workspace  This float vector contains the search direction
     *                                at the current iteration.  This vector has a
     *                                length of 'num_cols'.
     */
    std::valarray<T> work_bidiag_wrk_vec(T(0), A.num_cols());
    std::valarray<T> work_srch_dir_vec(T(0), A.num_cols());
    
    T bbnorm = T(0);
    T ddnorm = T(0);
    T xxnorm = T(0);
    
    T cs2 = T(-1);
    T sn2 = T(0);
    T zeta = T(0);
    T res = T(0);
    
    T cond_tol = ( input.cond_lim > T(0) ) ? T(1) / input.cond_lim : std::numeric_limits<T>::epsilon();
    
    
    /*
     *  Set up the initial vectors u and v for bidiagonalization.  These satisfy
     *  the relations
     *             BETA*u = b - A*x0
     *             ALPHA*v = A^T*u
     */
    if ( x.size() != A.num_cols() )
    {
        x.resize(A.num_cols());
    }
    /* Compute b - A*x0 and store in vector u which initially held vector b */
    std::valarray<T> rhs_vec = -b;
    A.multiply(rhs_vec, x);
    rhs_vec *= T(-1);
    
    /* compute Euclidean length of u and store as BETA */
    T alpha = T(0);
    T beta = nrm2( rhs_vec );
    
    if ( beta > T(0) )
    {
        /* scale vector u by the inverse of BETA */
        rhs_vec /= beta;
        
        /* Compute matrix-vector product A^T*u and store it in vector v */
        A.multiplyT(work_bidiag_wrk_vec, rhs_vec);
        
        /* compute Euclidean length of v and store as ALPHA */
        alpha = nrm2( work_bidiag_wrk_vec );
    }
    
    if ( alpha > T(0) )
    {
        /* scale vector v by the inverse of ALPHA */
        work_bidiag_wrk_vec /= alpha;
        
        /* copy vector v to vector w */
        work_srch_dir_vec = work_bidiag_wrk_vec;
    }
    
    output.mat_resid_norm = alpha * beta;
    output.resid_norm = beta;
    T bnorm = beta;
    /*
     *  If the norm || A^T r || is zero, then the initial guess is the exact
     *  solution.  Exit and report this.
     */
    if ( ( output.mat_resid_norm == T(0) ) && lsqr_fp_out )
    {
        output.b_norm = bnorm;
        lsqr_fp_out.setf(std::ios_base::scientific, std::ios_base::floatfield);
        lsqr_fp_out << "\tISTOP = " << std::setw(3) << output.term_flag
        << "\t\t\tITER = " << std::setw(9) << output.num_iters << "\n"
        << "	|| A ||_F = " << std::setw(13) << std::setprecision(5) << output.frob_mat_norm
        << "\tcond( A ) = " << std::setw(13) << std::setprecision(5) << output.mat_cond_num << "\n"
        << "	|| r ||_2 = " << std::setw(13) << std::setprecision(5) << output.resid_norm
        << "\t|| A^T r ||_2 = " << std::setw(13) << std::setprecision(5) << output.mat_resid_norm << "\n"
        << "	|| b ||_2 = " << std::setw(13) << std::setprecision(5) << output.b_norm
        << "\t|| x - x0 ||_2 = " << std::setw(13) << std::setprecision(5) << output.sol_norm << "\n"
        << std::endl;
        
        lsqr_fp_out << "  " << term_msg[output.term_flag] << "\n" << std::endl;
        
        return output;
    }
    
    T rhobar = alpha;
    T phibar = beta;
    /*
     *  If statistics are printed at each iteration, print a header and the initial
     *  values for each quantity.
     */
    if ( lsqr_fp_out )
    {
        lsqr_fp_out << "  ITER     || r ||    Compatible  ||A^T r|| / ||A|| ||r||  || A ||   cond( A )\n" << std::endl;
        lsqr_fp_out.setf(std::ios_base::scientific, std::ios_base::floatfield);
        lsqr_fp_out << std::setw(6) << output.num_iters << " "
        << std::setw(13) << std::setprecision(5) << output.resid_norm << " "
        << std::setw(10) << std::setprecision(2) << T(1) << " \t"
        << std::setw(10) << std::setprecision(2) << alpha / beta << " \t"
        << std::setw(10) << std::setprecision(2) << output.frob_mat_norm << " "
        << std::setw(10) << std::setprecision(2) << output.mat_cond_num << std::endl;
    }
    
    /*
     *  The main iteration loop is continued as long as no stopping criteria
     *  are satisfied and the number of total iterations is less than some upper
     *  bound.
     */
    while ( !output.term_flag )
    {
        output.num_iters++;
        /*
         *     Perform the next step of the bidiagonalization to obtain
         *     the next vectors u and v, and the scalars ALPHA and BETA.
         *     These satisfy the relations
         *                BETA*u  =  A*v  -  ALPHA*u,
         *                ALFA*v  =  A^T*u  -  BETA*v.
         */
        /* scale vector u by the negative of ALPHA */
        rhs_vec *= -alpha;
        
        /* compute A*v - ALPHA*u and store in vector u */
        A.multiply(rhs_vec, work_bidiag_wrk_vec);
        
        /* compute Euclidean length of u and store as BETA */
        beta = nrm2( rhs_vec );
        
        /* accumulate this quantity to estimate Frobenius norm of matrix A */
        bbnorm += sqr(alpha) + sqr(beta) + sqr(input.damp_val);
        
        if ( beta > T(0) )
        {
            /* scale vector u by the inverse of BETA */
            rhs_vec /= beta;
            
            /* scale vector v by the negative of BETA */
            work_bidiag_wrk_vec *= -beta;
            
            /* compute A^T*u - BETA*v and store in vector v */
            A.multiplyT(work_bidiag_wrk_vec, rhs_vec);
            
            /* compute Euclidean length of v and store as ALPHA */
            alpha = nrm2( work_bidiag_wrk_vec );
            
            if ( alpha > T(0) )
            {
                /* scale vector v by the inverse of ALPHA */
                work_bidiag_wrk_vec /= alpha;
            }
        }
        /*
         *     Use a plane rotation to eliminate the damping parameter.
         *     This alters the diagonal (RHOBAR) of the lower-bidiagonal matrix.
         */
        T cs1 = rhobar / std::sqrt( sqr(rhobar) + sqr(input.damp_val) );
        
        T psi = phibar * input.damp_val / std::sqrt( sqr(rhobar) + sqr(input.damp_val) );
        phibar *= cs1;
        /*
         *     Use a plane rotation to eliminate the subdiagonal element (BETA)
         *     of the lower-bidiagonal matrix, giving an upper-bidiagonal matrix.
         */
        T rho = std::sqrt( sqr(rhobar) + sqr(input.damp_val) + sqr(beta) );
        T cs = std::sqrt( sqr(rhobar) + sqr(input.damp_val) ) / rho;
        T sn = beta / rho;
        
        T theta = sn * alpha;
        rhobar = -cs * alpha;
        T phi = cs * phibar;
        phibar *= sn;
        T tau = sn * phi;
        /*
         *     Update the solution vector x, the search direction vector w, and the
         *     standard error estimates vector se.
         */
        /* update the solution vector x */
        x += ( phi / rho ) * work_srch_dir_vec;
        if ( input.all_positive )
        {
            for (std::size_t i=0; i<x.size(); i++)
            {
                if ( x[i] < T(0) ) x[i] = T(0);
            }
            //x = std::abs(x);
        }
        
        /* update the standard error estimates vector se */
        output.std_err_vec += sqr( std::valarray<T>(work_srch_dir_vec / rho) );
        
        /* accumulate this quantity to estimate condition number of A */
        ddnorm += (sqr( std::valarray<T>(work_srch_dir_vec / rho) )).sum();
        
        /* update the search direction vector w */
        work_srch_dir_vec = work_bidiag_wrk_vec - ( theta / rho ) * work_srch_dir_vec;
        
        /*
         *     Use a plane rotation on the right to eliminate the super-diagonal element
         *     (THETA) of the upper-bidiagonal matrix.  Then use the result to estimate
         *     the solution norm || x ||.
         */
        T delta = sn2 * rho;
        T gammabar = -cs2 * rho;
        T zetabar = (phi - delta * zeta) / gammabar;
        
        /* compute an estimate of the solution norm || x || */
        output.sol_norm = std::sqrt( xxnorm + sqr(zetabar) );
        
        T gamma = std::sqrt( sqr(gammabar) + sqr(theta) );
        cs2 = gammabar / gamma;
        sn2 = theta / gamma;
        zeta = (phi - delta * zeta) / gamma;
        
        /* accumulate this quantity to estimate solution norm || x || */
        xxnorm += sqr(zeta);
        /*
         *     Estimate the Frobenius norm and condition of the matrix A, and the
         *     Euclidean norms of the vectors r and A^T*r.
         */
        output.frob_mat_norm = std::sqrt( bbnorm );
        output.mat_cond_num = output.frob_mat_norm * std::sqrt( ddnorm );
        
        res += sqr(psi);
        output.resid_norm = std::sqrt( sqr(phibar) + res );
        
        output.mat_resid_norm = alpha * fabs( tau );
        /*
         *     Use these norms to estimate the values of the three stopping criteria.
         */
        T stop_crit_1 = output.resid_norm / bnorm;
        
        T stop_crit_2 = ( output.resid_norm > T(0) ) ? output.mat_resid_norm / ( output.frob_mat_norm * output.resid_norm ) : T(0);
        
        T stop_crit_3 = T(1) / output.mat_cond_num;
        
        T resid_tol = input.rel_rhs_err + input.rel_mat_err * output.mat_resid_norm * output.sol_norm / bnorm;
        
        T resid_tol_mach = std::numeric_limits<T>::epsilon()
        + std::numeric_limits<T>::epsilon() * output.mat_resid_norm * output.sol_norm / bnorm;
        /*
         *     Check to see if any of the stopping criteria are satisfied.
         *     First compare the computed criteria to the machine precision.
         *     Second compare the computed criteria to the the user specified precision.
         */
        /* iteration limit reached */
        if ( output.num_iters >= input.max_iter )
        {
            output.term_flag = 7;
        }
        
        /* condition number greater than machine precision */
        if ( stop_crit_3 <= std::numeric_limits<T>::epsilon() )
        {
            output.term_flag = 6;
        }
        /* least squares error less than machine precision */
        if ( stop_crit_2 <= std::numeric_limits<T>::epsilon() )
        {
            output.term_flag = 5;
        }
        /* residual less than a function of machine precision */
        if ( stop_crit_1 <= resid_tol_mach )
        {
            output.term_flag = 4;
        }
        
        /* condition number greater than CONLIM */
        if ( stop_crit_3 <= cond_tol )
        {
            output.term_flag = 3;
        }
        /* least squares error less than ATOL */
        if ( stop_crit_2 <= input.rel_mat_err )
        {
            output.term_flag = 2;
        }
        /* residual less than a function of ATOL and BTOL */
        if ( stop_crit_1 <= resid_tol )
        {
            output.term_flag = 1;
        }
        /*
         *  If statistics are printed at each iteration, print a header and the initial
         *  values for each quantity.
         */
        if ( lsqr_fp_out )
        {
            lsqr_fp_out.setf(std::ios_base::scientific, std::ios_base::floatfield);
            lsqr_fp_out << std::setw(6) << output.num_iters << " "
            << std::setw(13) << std::setprecision(5) << output.resid_norm << " "
            << std::setw(10) << std::setprecision(2) << stop_crit_1 << " \t"
            << std::setw(10) << std::setprecision(2) << stop_crit_2 << " \t"
            << std::setw(10) << std::setprecision(2) << output.frob_mat_norm << " "
            << std::setw(10) << std::setprecision(2) << output.mat_cond_num << std::endl;
        }
        /*
         *     The convergence criteria are required to be met on NCONV consecutive
         *     iterations, where NCONV is set below.  Suggested values are 1, 2, or 3.
         */
        if ( output.term_flag == 0 )
        {
            term_iter = -1;
        }
        
        long term_iter_max = 1;
        term_iter++;
        
        if ( ( term_iter < term_iter_max ) && ( output.num_iters < input.max_iter ) )
        {
            output.term_flag = 0;
        }
    } /* end while loop */
    /*
     *  Finish computing the standard error estimates vector se.
     */
    T temp = (  A.num_rows() > A.num_cols() ) ? (T)( A.num_rows() - A.num_cols() ) : T(1);
    
    if ( sqr(input.damp_val) > T(0) )
    {
        temp = (T) ( A.num_rows() );
    }
    
    temp = output.resid_norm / std::sqrt( temp );
    
    /* update the standard error estimates vector se */
    output.std_err_vec = temp * std::sqrt(output.std_err_vec);
    
    /*
     *  If statistics are printed at each iteration, print the statistics for the
     *  stopping condition.
     */
    if ( lsqr_fp_out )
    {
        output.b_norm = bnorm;
        lsqr_fp_out.setf(std::ios_base::scientific, std::ios_base::floatfield);
        lsqr_fp_out << "\n\tISTOP = " << std::setw(3) << output.term_flag
        << "\t\t\tITER = " << std::setw(9) << output.num_iters << "\n"
        << "	|| A ||_F = " << std::setw(13) << std::setprecision(5) << output.frob_mat_norm
        << "\tcond( A ) =      " << std::setw(13) << std::setprecision(5) << output.mat_cond_num << "\n"
        << "	|| r ||_2 = " << std::setw(13) << std::setprecision(5) << output.resid_norm
        << "\t|| A^T r ||_2 =  " << std::setw(13) << std::setprecision(5) << output.mat_resid_norm << "\n"
        << "	|| b ||_2 = " << std::setw(13) << std::setprecision(5) << output.b_norm
        << "\t|| x - x0 ||_2 = " << std::setw(13) << std::setprecision(5) << output.sol_norm << "\n"
        << std::endl;
        
        lsqr_fp_out << "  " << term_msg[output.term_flag] << "\n" << std::endl;
    }
    
    return output;
}





#endif /* _LSQR_ */
