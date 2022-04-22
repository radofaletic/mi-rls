/**
 for testing ray tracing AND inversion
 */

#include <functional>
#include <cmath>
#include <string>
#include <valarray>
#include <vector>

#if DEBUG >= 5
#ifndef USE_MESSAGES
#define USE_MESSAGES
#endif /* USE_MESSAGES */
#endif /* DEBUG */

#include "conversions.h"
#include "fortran_io.h"
#include "front-end.h"
#include "grid.h"
#include "line.h"
#include "matrix_utilities.h"
#include "sparse_matrix.h"
#include "tomography.h"

#ifdef USE_MESSAGES
#include <fstream>
std::ofstream messages;
#endif /* USE_MESSAGES */

int main(int argc, char* argv[])
{
#ifdef USE_MESSAGES
	messages.open((std::string(argv[0])+std::string(".messages")).c_str());
#endif /* USE_MESSAGES */
	
	/*
	 read in grid and a-priori data
	 */
	grid_input gridinputs;
	gridinputs.type() = structured;
	gridinputs.load_grid() = true;
	gridinputs.format() = Binary;
	gridinputs.precision() = Single;
	gridinputs.multidomain() = false;
	gridinputs.blanking() = false;
	gridinputs.gridfile() = "data/tests/stark_3d.PBG";
	gridinputs.datafile() = "data/tests/stark_3d.PBS";
	gridinputs.qdata() = 1;
	
	grid<float> test_grid(gridinputs);
	test_grid.read_data(gridinputs);
	test_grid.give_dataname("stark_3d");
	std::string the_dataname = test_grid.get_dataname();
	std::string the_filename = test_grid.get_filename();
	SparseMatrix<float> A(0,test_grid.ncells());
	std::valarray<float> b;
	std::valarray<float> cellw;
	
	test_grid.set_basis(l2norm);
	const std::string mname = "saved_data.A";
	const std::string dname = "saved_data.b";
	const std::string cname = "saved_data.c";
	const int rotations = 12;
	//const int rotations = 3;
	
	if ( !yesno("read existing matrix",false) )
	{
		
		//
		// large objects
		//
		std::valarray<std::size_t> cellnumbers;
		std::valarray<float> weights;
		
		//
		// generate all of the rays
		//
		message("generating all rays");
        std::size_t meshsize = 32;
		std::valarray<float>* ratt = new std::valarray<float>(test_grid.scale()/10.0f,test_grid.dim());
		std::valarray<float>* ratu = new std::valarray<float>(test_grid.min() + *ratt);
		std::valarray<float>* ratv = new std::valarray<float>(test_grid.max() - *ratt);
        std::size_t Nrows;
        std::size_t Ncols;
		std::vector< line<float> > ray = Tomography::synthetic_rays(rotations, meshsize, *ratu, *ratv, Nrows, Ncols);
		delete ratt;
		delete ratu;
		delete ratv;
		
		//
		// do the tomographic projections
		//
		b.resize(ray.size(), 0.0f);
		std::valarray<bool> blanks;
		test_grid.give_dataname(the_dataname+"-projected");
		test_grid.write(the_filename+"-projected",SaveGrid,Binary);
		Tomography::projection(test_grid, ray, blanks, b, A, walkfast, true, true, true);
		
		//
		// write tomographic data to files
		//
		message("writing final projections");
		test_grid.write(the_filename+"-projected", SaveData, Binary);
		test_grid.write(the_filename+"-projected", SaveGrid, Binary);
		message("writing PNG file");
		Tomography::pngwrite(the_filename+"-projected.png", Nrows, Ncols, b);
		
		cellw.resize(test_grid.ncells(), 0.0f);
		for (std::size_t i=0; i<test_grid.ncells(); i++)
		{
			cellw[i] = test_grid[i];
		}
		test_grid.clear_data();
		
		message("writing matrix to file");
		A.write(mname, Unformatted);
		message("writing residuals to file");
		Fortran::fdump_vector(dname, b, false, Unformatted);
		message("writing projected weights to file");
		Fortran::fdump_vector(cname, cellw, false, Unformatted);
	}
	else
	{
		message("loading existing matrix");
		A.read(mname,Unformatted);
		bool bs = false;
		Fortran::fretrieve_vector(dname, b, bs, Unformatted);
		Fortran::fretrieve_vector(cname, cellw, bs, Unformatted);
	}
	
	//
	// damping
	//
	message("adding damping equations");
	AddDamping(A, b, test_grid, question("matrix damping value",float(0.001)));
	
	//
	// smoothing
	//
	message("adding neighbour smoothing equations");
	AddSmoothing(A, b, test_grid, question("matrix smoothing value",float(0.1)));
	
	//
	// get rid of now irrelevant neighbour and surface data
	//
	test_grid.clear_aux();
	
	//
	// solve Ax=b
	//
	message("inverting " + std::string(typeid(A).name()));
	std::valarray<float> x(float(0), test_grid.ncells());
	x = A.iterate(b, 25000, 5);
	
	for (std::size_t i=0; i<test_grid.ncells(); i++)
	{
		test_grid[i] = test_grid.applybasis(i, x[i]);
	}
	
	message("writing solution");
	test_grid.give_dataname(the_dataname);
	test_grid.write_data(the_filename+"-testoutput", Binary, cellw);
	test_grid.write(the_filename+"-testoutput", SaveGrid, Binary);
	
	message("FINISHED running " + std::string(argv[0]));
	message("please wait until the programme exits");
	return 0;
}
