/**
 for testing sparse matrix routines
 */

#include <valarray>

#if DEBUG >= 5
#ifndef USE_MESSAGES
#define USE_MESSAGES
#endif /* USE_MESSAGES */
#endif /* DEBUG */

#include "angles.h"
#include "conversions.h"
#include "file.h"
#include "front-end.h"
#include "matrix.h"
#include "sparse_matrix.h"
#include "lsqr.h"

#ifdef USE_MESSAGES
#include <fstream>
std::ofstream messages;
#endif /* USE_MESSAGES */

int main(int argc, char* argv[])
{
#ifdef USE_MESSAGES
	messages.open((std::string(argv[0])+std::string(".messages")).c_str());
#endif /* USE_MESSAGES */
	
	message("\n     ----- ----- -----     \n");
	
	SparseMatrix<float> A(3,4);
	message("matrix size = "+std::to_string(A.row_dim())+"x"+std::to_string(A.col_dim()));
	message(indent(A.print()));
	Matrix<float> M(3,4,float(1));
	message("should be\n"+indent(M.print()));
	A.write("mtest_output_3x3.matrix",Unformatted);
	
	message("\n     ----- ----- -----     \n");
	
	A.clear();
	message("matrix size = "+std::to_string(A.row_dim())+"x"+std::to_string(A.col_dim()));
	message(indent(A.print()));
	Matrix<float> Z;
	message("should be\n"+indent(Z.print()));
	
	message("\n     ----- ----- -----     \n");
	
	A.read("mtest_output_3x3.matrix",Unformatted);
	message("matrix size = "+std::to_string(A.row_dim())+"x"+std::to_string(A.col_dim()));
	message(indent(A.print()));
	M.init(3, 4, float(1));
	message("should be\n"+indent(M.print()));
	
	message("\n     ----- ----- -----     \n");
	
	std::valarray<std::size_t> el(3);
	el[0] = 5; el[1] = 6; el[2] = 8;
	std::valarray<float> da(3);
	da[0] = 6.6667; da[1] = 5.5; da[2] = Angle::pi;
	
	A.AddRow(el,da);
	message("matrix size = "+std::to_string(A.row_dim())+"x"+std::to_string(A.col_dim()));
	message(indent(A.print()));
	M.init(4,9,float(0));
	M(0,0) = M(1,1) = M(2,2) = float(1);
	M(3,5) = da[0];
	M(3,6) = da[1];
	M(3,8) = da[2];
	message("should be\n"+indent(M.print()));
	A.write("mtest_output_4x9.matrix",Formatted);
	
	message("\n     ----- ----- -----     \n");
	
	A.clear();
	message("matrix size = "+std::to_string(A.row_dim())+"x"+std::to_string(A.col_dim()));
	message(indent(A.print()));
	message("should be\n"+indent(Z.print()));
	
	message("\n     ----- ----- -----     \n");
	
	M.init(3,3,float(0));
	el.resize(3); el[0] = 0; el[1] = 1; el[2] = 2;
	da.resize(3); da[0] = 1; da[1] = 2; da[2] = 8;
	A.AddRow(el,da);
	M(0,0) = da[0];
	M(0,1) = da[1];
	M(0,2) = da[2];
	el.resize(2); el[0] = 0; el[1] = 2;
	da.resize(2); da[0] = 5; da[1] = Angle::pi;
	A.AddRow(el,da);
	M(1,0) = da[0];
	M(1,2) = da[1];
	el.resize(2); el[0] = 1; el[1] = 2;
	da.resize(2); da[0] = -2; da[1] = -2;
	A.AddRow(el,da);
	M(2,1) = da[0];
	M(2,2) = da[1];
	message("matrix size = "+std::to_string(A.row_dim())+"x"+std::to_string(A.col_dim()));
	message(indent(A.print()));
	message("should be\n"+indent(M.print()));
	A.write("mtest_output_3x3solve.matrix",Unformatted);
	
	message("\n     ----- ----- -----     \n");
	
	A.clear();
	message("matrix size = "+std::to_string(A.row_dim())+"x"+std::to_string(A.col_dim()));
	message(indent(A.print()));
	message("should be\n"+indent(Z.print()));
	
	message("\n     ----- ----- -----     \n");
	
	A.read("mtest_output_3x3solve.matrix",Unformatted);
	message("matrix size = "+std::to_string(A.row_dim())+"x"+std::to_string(A.col_dim()));
	message(indent(A.print()));
	message("should be\n"+indent(M.print()));
	
	message("\n     ----- ----- -----     \n");
	
	da.resize(3); da[0] = 9.2; da[1] = -50; da[2] = 0;
	message("compute v = A*"+vtos(da));
	std::valarray<float> v = A * da;
	message(indent("v = "+vtos(v)));
	da[0] = -90.8; da[1] = 46; da[2] = 100;
	message("should be\n"+indent(vtos(da)));
	
	message("\n     ----- ----- -----     \n");
	
	da[0] = 9.2; da[1] = -50; da[2] = 0;
	v[0] = 2; v[1] = 2; v[2] = Angle::pi;
	message("compute u = "+vtos(v)+" + A*"+vtos(da));
	A.multiply(v, da);
	message(indent("u = "+vtos(v)));
	da[0] = -88.8; da[1] = 48; da[2] = 100+Angle::pi;
	message("should be\n"+indent(vtos(da)));
	
	message("\n     ----- ----- -----     \n");
	
	da[0] = 9.2; da[1] = -50; da[2] = 0;
	v[0] = 2; v[1] = 2; v[2] = Angle::pi;
	message("compute v = "+vtos(v)+" + A(T)*"+vtos(da));
	A.multiplyT(v, da);
	message(indent("u = "+vtos(v)));
	da[0] = -238.8; da[1] = 20.4; da[2] = 73.6 - 50*Angle::pi + Angle::pi;
	message("should be\n"+indent(vtos(da)));
	
	message("\n     ----- ----- -----     \n");
	
	da[0] = 9.2; da[1] = -50; da[2] = 0;
	message("solve A*x="+vtos(da));
	//std::valarray<float> x = A.iterate(da,question("how many iterations to perform",100),question("order of convergence desired",10));
	std::valarray<float> x(A.cols());
	lsqr_input<float> linput(0,0,0,0,question("how many iterations to perform",100));
	LSQR(A, x, da, linput);
	message("x = "+vtos(x));
	x[0] = -12.24579668; x[1] = -3.57429945; x[2] = 3.574299447;
	message("should be\n    "+vtos(x));
	message("FINISHED running "+std::string(argv[0]));
}
