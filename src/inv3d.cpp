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
	std::string ip3dgfilename = "input.PBG";
	std::string ip3dqfilename = "input.PBS";
	short ip3dqi = 1;
	std::string op3dfilename = "output";
	size_t xresolution = 32;
	size_t yresolution = 32;
	size_t nangles = 12;
	Angle::axes rot_axis = Angle::X;
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
			message("--input_g=<inputfile>\n\tthe Plot3D grid file ("+ip3dgfilename+")");
			message("--input_q=<inputfile>\n\tthe Plot3D Q file ("+ip3dqfilename+")");
			message("--input_i=<n>\n\tthe Plot3D Q variable number between 1 and 5 ("+ntos(ip3dqi)+")");
			message("--output=<outputfile>\n\tthe output Plot3D file ("+op3dfilename+")");
			message("--xresolution=<res>\n\tx resolution of the projection plane ("+ntos(xresolution)+")");
			message("--yresolution=<res>\n\ty resolution of the projection plane ("+ntos(yresolution)+")");
			message("--angles=<nangles>\n\tthe number of angles, or projections ("+ntos(nangles)+")");
			message("--axis=<axis>\n\taxis of rotation (X)");
			message("--intermediates=on/off\n\twrite PNG files of intermediate steps (off)");
			return 1;
		}
		else if ( fswitch[i].var("input_g") )
		{
			ip3dgfilename = fswitch[i].val();
		}
		else if ( fswitch[i].var("input_q") )
		{
			ip3dqfilename = fswitch[i].val();
		}
		else if ( fswitch[i].var("input_i") )
		{
			std::istringstream itmp(fswitch[i].val());
			itmp >> ip3dqi;
			if ( ip3dqi < 1 || 5 < ip3dqi )
			{
				message("\"input_i\" must be between 1 and 5.\nUse the --help option to learn more.");
				throw; return 1;
			}
		}
		else if ( fswitch[i].var("output") || fswitch[i].var("o") )
		{
			op3dfilename = fswitch[i].val();
		}
		else if ( fswitch[i].var("xresolution") || fswitch[i].var("xres") )
		{
			std::istringstream choice(fswitch[i].val());
			choice >> xresolution;
		}
		else if ( fswitch[i].var("yresolution") || fswitch[i].var("yres") )
		{
			std::istringstream choice(fswitch[i].val());
			choice >> yresolution;
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
		else if ( fswitch[i].var("axis") )
		{
			if ( fswitch[i].val() == std::string("X") || fswitch[i].val() == std::string("YZ") )
			{
				rot_axis = Angle::X;
			}
			else if ( fswitch[i].val() == std::string("Y") || fswitch[i].val() == std::string("ZX") )
			{
				rot_axis = Angle::Y;
			}
			//else if ( fswitch[i].val() == std::string("Z") || fswitch[i].val() == std::string("XY") )
			//  {
			//    rot_axis = Angle::Z;
			//  }
			else
			{
				message("\"axis\" must be one of X or Y.\nUse the --help option to learn more.");
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
	// read grid
	message("reading the grid");
	grid_input mygrid;
	mygrid.type() = structured;
	mygrid.format() = Binary;
	mygrid.precision() = Single;
	mygrid.multidomain() = false;
	mygrid.blanking() = false;
	mygrid.load_grid() = true;
	mygrid.gridfile() = ip3dgfilename;
	mygrid.datafile() = ip3dqfilename;
	mygrid.qdata() = ip3dqi;
	grid<real> the_grid(mygrid);
	the_grid.read_data(mygrid);
	
	std::string gdn = op3dfilename;
	the_grid.give_dataname(op3dfilename);
	
	std::valarray<real> angles;
	
	// generate projection lines
	message("setting up rays");
	real dlength = 1;
	real dlength_s = 1;
	std::valarray<real> min2(2);
	std::valarray<real> max2(2);
	std::valarray<real> min3 = the_grid.min();
	std::valarray<real> max3 = the_grid.max();
	switch(rot_axis)
	{
		case Angle::X: case Angle::YZ:
			min2[0] = min3[1];
			min2[1] = min3[2];
			max2[0] = max3[1];
			max2[1] = max3[2];
			dlength = norm(&min2, &max2);
			dlength_s = std::abs(max3[0]-min3[0]);
			break;
		case Angle::Y: case Angle::ZX:
			min2[0] = min3[2];
			min2[1] = min3[0];
			max2[0] = max3[2];
			max2[1] = max3[0];
			dlength = norm(&min2, &max2);
			dlength_s = std::abs(max3[1]-min3[1]);
			break;
	}
	real scale = dlength / (xresolution + 1);
	real scale_s = dlength_s / (yresolution + 1);
	dlength /= 2;
	dlength_s /= 2;
	std::vector< std::valarray<real> > ipoints(xresolution*yresolution, std::valarray<real>(3));
	size_t counter = 0;
	for (size_t j=1; j<xresolution+1; j++)
	{
		for (size_t i=1; i<yresolution+1; i++)
		{
			ipoints[counter][0] = i * scale_s - dlength_s;
			ipoints[counter][1] = j * scale - dlength;
			ipoints[counter][2] = 0;
			ipoints[counter] += the_grid.center();
			counter++;
		}
	}
	
	angles.resize(nangles);
	for (size_t i=0; i<angles.size(); i++)
	{
		angles[i] = ((real)(180*i))/((real)(nangles));
	}
	
	Rotation<real> rot(3);
	rot.set_origin(the_grid.center());
	std::valarray<real> islope(3);
	islope[0] = 0;
	islope[1] = 0;
	islope[2] = 1;
	std::vector< line<real> > rays(0);
	for (size_t i=0; i<angles.size(); i++)
	{
		rot.reset(angles[i], rot_axis);
		std::valarray<real> slope = rot.O(islope);
		for (size_t j=0; j<ipoints.size(); j++)
		{
			line<real> tline(slope, rot(ipoints[j]));
			rays.push_back(tline);
		}
	}
	
	// project rays
	message("projecting rays");
	SparseMatrix<real> A(0,the_grid.ncells());
	std::valarray<real> b;
	std::valarray<bool> blanks(true, rays.size());
	Tomography::projection(the_grid, rays, blanks, b, A, walkfast, true, true, true);
	
	if ( write_intermediates ) // save projections
	{
		std::string o = op3dfilename;
		o += "_1-projections.png";
		message("saving '"+o+"'");
		std::valarray<real> otb(real(0), blanks.size());
		otb[blanks] = b;
		Tomography::pngwrite(o, yresolution, xresolution, otb, rot_axis, angles, scale, false, true);
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
	lsqr_input<real> linput(0,0,0,0,100);
	LSQR(A, x, b, linput);
	message("uncompressing matrix");
	A.Uncompress(x);
	
	for (size_t i=0; i<x.size(); i++)
	{
		the_grid[i] = x[i];
	}
	
	// write output
	message("saving '"+op3dfilename+"'");
	the_grid.write(op3dfilename, SaveGrid, Binary);
	the_grid.write(op3dfilename, SaveData, Binary);
	
	return 0;
}
