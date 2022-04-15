/*
 A programme to perform a linear tomographic inversion on any greyscale PNG file
 
 Rado Faletic
 8th July 2004
 */

/*
 #undef DEBUG
 */

#include <cmath>
#include <string>
#include <valarray>
#include <vector>
#include "angles.h"
#include "argv.h"
#include "front-end.h"
#include "grid.h"
#include "line.h"
#include "lsqr.h"
#include "matrix_utilities.h"
#include "sparse_matrix.h"
#include "tomography.h"

typedef double real;

int main(int argc, char* argv[])
{
	std::string ipngfilename = "input.png";
	std::string opngfilename = "output.png";
	size_t resolution = 256;
	size_t nangles = 180;
	bool write_intermediates = false;
	
	std::vector<args> fswitch = get_args(argc, argv);
	if ( !fswitch.size() )
	{
		fswitch.resize(1);
		fswitch[0].var() = "help";
		fswitch[0].val() = "";
	}
	for (size_t i=0; i<fswitch.size(); i++)
	{
		if ( fswitch[i].var("help") )
		{
			message("\nby Rado Faletic 2004\n");
			message("below is a list of flags:\n");
			message("--input=<inputfile>\n\tthe input PNG file ("+ipngfilename+")");
			message("--output=<outputfile>\n\tthe output PNG file ("+opngfilename+")");
			message("--resolution=<res>\n\tthe number of data points in each projection ("+ntos(resolution)+")");
			message("--angles=<nangles>\n\tthe number of angles, or projections ("+ntos(nangles)+")");
			message("--intermediades=on/off\n\twrite PNG files of intermediate steps (off)");
			return 1;
		}
		else if ( fswitch[i].var("input") || fswitch[i].var("i") )
		{
			ipngfilename = fswitch[i].val();
		}
		else if ( fswitch[i].var("output") || fswitch[i].var("o") )
		{
			opngfilename = fswitch[i].val();
		}
		else if ( fswitch[i].var("resolution") || fswitch[i].var("res") )
		{
			std::istringstream choice(fswitch[i].val());
			choice >> resolution;
		}
		else if ( fswitch[i].var("angles") )
		{
			std::istringstream choice(fswitch[i].val());
			choice >> nangles;
			if ( !nangles )
			{
				message("\"angles\" must be non-zero.\nUse the --help option to learn more.");
				throw; return 1;
			}
		}
		else if ( fswitch[i].var("intermediates") )
		{
			if ( fswitch[i].val() == std::string("on") || fswitch[i].val() == std::string("ON") ||
				fswitch[i].val() == std::string("yes") || fswitch[i].val() == std::string("YES") )
			{
				write_intermediates = true;
			}
			else
			{
				write_intermediates = false;
			}
		}
		else
		{
			message("Unrecognised option '"+fswitch[i].var()+"'.\nUse the --help option to learn more.");
			throw; return 1;
		}
	}
	if ( ipngfilename.substr(ipngfilename.size()-4,4) != ".png" && ipngfilename.substr(ipngfilename.size()-4,4) != ".PNG" )
	{
		ipngfilename += ".png";
	}
	if ( opngfilename.substr(opngfilename.size()-4,4) != ".png" && opngfilename.substr(opngfilename.size()-4,4) != ".PNG" )
	{
		opngfilename += ".png";
	}
	
	// read PNG image file
	message("reading '"+ipngfilename+"'");
	size_t Nrows, Ncols;
	std::valarray<real> input_data;
	Angle::axes sino_axis;
	std::valarray<real> angles;
	real scale;
	bool realdata;
	std::valarray<bool> blanks;
	Tomography::pngread(ipngfilename, Nrows, Ncols, input_data, blanks, sino_axis, angles, scale, realdata);
	if ( Nrows*Ncols != input_data.size() )
	{
		message("'"+ipngfilename+"' seems to be a multi-image tomographic file, this software will only read a single-image file.");
		throw; return 1;
	}
	
	// set up grid
	message("setting up the grid");
	grid_input mygrid;
	mygrid.type() = structured;
	mygrid.load_grid() = false;
	mygrid.g_nX() = Ncols;
	mygrid.g_nY() = Nrows;
	mygrid.g_nZ() = 0;
	mygrid.g_scale() = scale;
	grid<real> the_grid(mygrid);
	
	the_grid.put_adata(input_data);
	std::string gdn = opngfilename.substr(0,opngfilename.size()-4);
	the_grid.give_dataname(gdn);
	
	// generate projection lines
	message("setting up projection rays");
	if ( !resolution )
	{
		resolution = size_t(std::ceil(std::sqrt((real)(Ncols*Ncols+Nrows*Nrows))));
	}
	real dlength = scale * std::sqrt((real)(Ncols*Ncols+Nrows*Nrows));
	scale = dlength / ( resolution + 1 );
	dlength /= 2;
	std::vector< std::valarray<real> > ipoints(resolution+2);
	for (size_t i=0; i<ipoints.size(); i++)
	{
		ipoints[i].resize(2);
		ipoints[i][0] = i * scale - dlength;
		ipoints[i][1] = 0;
		ipoints[i] += the_grid.center();
	}
	ipoints.erase(ipoints.end()-1);
	ipoints.erase(ipoints.begin());
	angles.resize(nangles);
	for (size_t i=0; i<angles.size(); i++)
	{
		angles[i] = ((real)(180*i))/((real)(nangles));
	}
	
	Rotation<real> rot(2);
	rot.set_origin(the_grid.center());
	std::valarray<real> islope(2);
	islope[0] = 0;
	islope[1] = 1;
	std::vector< line<real> > rays(0);
	for (size_t i=0; i<angles.size(); i++)
	{
		rot.reset(angles[i], Angle::XY);
		std::valarray<real> slope = rot.O(islope);
		for (size_t j=0; j<ipoints.size(); j++)
		{
			line<real> tline(slope, rot(ipoints[j]));
			rays.push_back(tline);
		}
	}
	
	// create tomograms
	message("creating tomograms");
	SparseMatrix<real> A(0,the_grid.ncells());
	std::valarray<real> b;
	Tomography::projection(the_grid, rays, blanks, b, A, walkfast, true, true, true);
	
	if ( write_intermediates ) // save sinogram
	{
		std::string o = opngfilename.substr(0,opngfilename.size()-4);
		o += "_1-tomograms.png";
		message("saving '"+o+"'");
		std::valarray<real> otb(real(0), blanks.size());
		otb[blanks] = b;
		Tomography::pngwrite(o, nangles, resolution, otb, sino_axis, angles, scale, false, true);
	}
	
	// add smoothing
	//message("adding smoothing equations");
	//AddSmoothing(A, b, the_grid, real(0.1));
	the_grid.clear_aux();
	
	// doing matrix inversion
	std::valarray<real> x(real(0), A.cols());
	message("compressing matrix");
	A.Compress(x);
	message("inverting matrix");
	lsqr_input<real> linput(0,0,0,0,1000);
	LSQR(A, x, b, linput);
	message("uncompressing matrix");
	A.Uncompress(x);
	
	// write output
	message("saving '"+opngfilename+"'");
	Tomography::pngwrite(opngfilename, Nrows, Ncols, x, sino_axis, angles, scale, true, true);
	
	return 0;
}
