/**
 A programme to create a sinogram from any PNG image file. If the PNG contains a "data" tag filled with IEEE format numbers, as long as it isn't a multiple image file, then these values will be used. If there is no "data" tag then this programme looks for the appropriate pCAL tag, and scales the data appropriately. Failing either of these two options, this programme will simply use the scalar image values.
 
 Rado Faletic
 8th August 2004
 22nd April 2022
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
#include "sparse_matrix.h"
#include "tomography.h"

int main(int argc, char* argv[])
{
	std::string ipngfilename = "input.png";
	std::string opngfilename = "output.png";
    std::size_t resolution = 0;
    std::size_t nangles = 180;
	std::string stretch = "outer";
	
	std::vector<args> fswitch = get_args(argc, argv);
	if ( fswitch.size() == 0 )
	{
		fswitch.resize(1);
		fswitch[0].var() = "help";
		fswitch[0].val() = "";
	}
	for (std::size_t i=0; i<fswitch.size(); i++)
	{
		if ( fswitch[i].var("help") )
		{
			message("\nby Rado Faletic  2003, 2004\n");
			message("below is a list of flags:\n");
			message("--input=<inputfile>\n\tthe input PNG file ("+ipngfilename+")");
			message("--output=<outputfile>\n\tthe output PNG file ("+opngfilename+")");
			message("--resolution=<res>\n\tthe number of data points in each projection (determined by file)");
			message("--angles=<nangles>\n\tthe number of angles, or projections ("+std::to_string(nangles)+")");
			message("--stretch=inner/middle/outer\n\thow to structure the projections ("+stretch+")");
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
			if ( !resolution )
			{
				message("\"resolution\" must be non-zero.\nUse the --help option to learn more.");
				throw;
                return 1;
			}
		}
		else if ( fswitch[i].var("angles") )
		{
			std::istringstream choice(fswitch[i].val());
			choice >> nangles;
			if ( !nangles )
			{
				message("\"angles\" must be non-zero.\nUse the --help option to learn more.");
				throw;
                return 1;
			}
		}
		else if ( fswitch[i].var("stretch") || fswitch[i].var("s") )
		{
			stretch = fswitch[i].val();
		}
		else
		{
			message("Unrecognised option '"+fswitch[i].var()+"'.\nUse the --help option to learn more.");
			throw;
            return 1;
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
    std::size_t Nrows, Ncols;
	std::valarray<double> input_data;
	Angle::axes sino_axis;
	std::valarray<double> angles;
    double scale;
	bool realdata;
	std::valarray<bool> blanks;
	Tomography::pngread(ipngfilename, Nrows, Ncols, input_data, blanks, sino_axis, angles, scale, realdata);
	if ( Nrows*Ncols != input_data.size() )
	{
		message("'"+ipngfilename+"' seems to be a multi-image tomographic file, this software will only read a single-image file.");
		throw;
        return 1;
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
	grid<double> the_grid(mygrid);
	
	the_grid.put_adata(input_data);
	the_grid.give_dataname(opngfilename.substr(0,opngfilename.size()-4));
	
	// generate projection lines
	message("setting up projection rays");
    double scale_d = scale;
    std::size_t restmp = std::size_t(std::sqrt(double(Nrows*Nrows+Ncols*Ncols))+0.5);
	if ( stretch == "inner" )
	{
		restmp = std::min(Ncols, Nrows);
	}
	else if ( stretch == "middle" )
	{
		restmp = std::max(Ncols, Nrows);
	}
	else if ( stretch == "outer" )
	{
		restmp = std::size_t(std::sqrt(double(Nrows*Nrows+Ncols*Ncols))+0.5);
	}
	if ( !resolution )
	{
		resolution = restmp;
	}
	else
	{
		scale_d *= double(restmp) / resolution;
	}
	std::vector< std::valarray<double> > ipoints(0);
	for (std::size_t i=0; i<resolution; i++)
	{
		std::valarray<double> inode = the_grid.center();
		inode[0] += ( double(i) - double(resolution-1) / 2 ) * scale_d;
		ipoints.push_back(inode);
	}
	angles.resize(nangles);
	for (std::size_t i=0; i<angles.size(); i++)
	{
		angles[i] = ((double)(180*i))/((double)(nangles));
	}
	
	Rotation<double> rot(2);
	rot.set_origin(the_grid.center());
	std::valarray<double> islope(2);
	islope[0] = 0;
	islope[1] = 1;
	std::vector< line<double> > rays(0);
	for (std::size_t i=0; i<angles.size(); i++)
	{
		rot.reset(angles[i], Angle::XY);
		std::valarray<double> slope = rot.O(islope);
		for (std::size_t j=0; j<ipoints.size(); j++)
		{
			line<double> tline(slope, rot(ipoints[j]));
			rays.push_back(tline);
		}
	}
	
	// create tomograms
	message("creating tomograms from "+std::to_string(rays.size())+" rays");
	SparseMatrix<double> A(0,the_grid.ncells());
	std::valarray<double> otb;
	Tomography::projection(the_grid, rays, blanks, otb, A, walkfast, true, false, true);
	std::valarray<double> b(double(0), blanks.size());
	b[blanks] = otb;
	otb.resize(0);
	
	// save sinogram
	message("saving '"+opngfilename+"'");
	Tomography::pngwrite(opngfilename, angles.size(), ipoints.size(), b, sino_axis, angles, scale, true, true);
	
	return 0;
}
