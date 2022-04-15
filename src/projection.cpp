/*
 A programme to create projections through any single domain Plot3D file
 
 Rado Faletic
 11th August 2004
 */

/*
 #undef DEBUG
 */

#include <cmath>
#include <fftw3.h>
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

typedef double real;

int main(int argc, char* argv[])
{
	std::string ip3dgfilename = "input.PBG";
	std::string ip3dqfilename = "input.PBS";
	short ip3dqi = 1;
	std::string output = "output.png";
	size_t resolution = 0;
	size_t nangles = 12;
	Angle::axes rot_axis = Angle::X;
	std::string stretch = "middle";
	
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
			message("--input=<inputfile>\n\tthe Plot3D file ("+ip3dgfilename+","+ip3dqfilename+","+ntos(ip3dqi)+")");
			message("--output=<outputfile>\n\tthe output PNG file ("+output+")");
			message("--resolution=<res>\n\tresolution of the projection plane ("+ntos(resolution)+")");
			message("--angles=<nangles>\n\tthe number of angles, or projections ("+ntos(nangles)+")");
			message("--axis=<axis>\n\taxis of rotation (X)");
			message("--stretch=inner/middle/outer\n\taxis of rotation ("+stretch+")");
			return 1;
		}
		else if ( fswitch[i].var("input") )
		{
			std::string tmp = fswitch[i].val();
			std::string::size_type get_c = tmp.find(",");
			if ( get_c != tmp.npos )
			{
				ip3dgfilename = tmp.substr(0,get_c);
				get_c++;
				tmp = tmp.substr(get_c, tmp.size()-get_c);
				get_c = tmp.find(",");
				if ( get_c != tmp.npos )
				{
					ip3dqfilename = tmp.substr(0,get_c);
					get_c++;
					tmp = tmp.substr(get_c, tmp.size()-get_c);
					std::istringstream itmp(tmp);
					itmp >> ip3dqi;
					if ( ip3dqi < 1 || 5 < ip3dqi )
					{
						message("\"input_i\" must be between 1 and 5.\nUse the --help option to learn more.");
						throw; return 1;
					}
				}
				else
				{
					ip3dqfilename = tmp;
				}
			}
			else
			{
				ip3dgfilename = tmp;
				ip3dqfilename = tmp.substr(0,tmp.size()-3) + "PBS";
			}
		}
		else if ( fswitch[i].var("output") || fswitch[i].var("o") )
		{
			output = fswitch[i].val();
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
		else if ( fswitch[i].var("stretch") )
		{
			stretch = fswitch[i].val();
		}
		else
		{
			message("Unrecognised option '"+fswitch[i].var()+"'.\nUse the --help option to learn more.");
			throw; return 1;
		}
	}
	
	if ( output.substr(output.size()-4,4) != ".png" && output.substr(output.size()-4,4) != ".PNG" )
	{
		output += ".png";
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
	the_grid.give_dataname(output.substr(0,output.size()-4));
	
	std::valarray<real> angles(real(0), nangles);
	for (size_t i=0; i<angles.size(); i++)
	{
		angles[i] = ((real)(180*i))/((real)(nangles));
	}
	
	// write grid to PNG
	message("writing the grid to PNG");
	std::valarray<real> x = the_grid.get_adata();
	std::valarray<real> xo(real(0), x.size());
	size_t Nx = the_grid.nex();
	size_t Ny = the_grid.ney();
	size_t Nz = the_grid.nez();
	Nx = ( Nx <= 1 ) ? 1 : Nx - 1;
	Ny = ( Ny <= 1 ) ? 1 : Ny - 1;
	Nz = ( Nz <= 1 ) ? 1 : Nz - 1;
	size_t N1 = 0;
	size_t N2 = 0;
	switch(rot_axis)
	{
		case Angle::X: case Angle::YZ:
			for (size_t k=0; k<Nx; k++)
			{
				for (size_t j=0; j<Ny; j++)
				{
					xo[std::slice(k*Ny*Nz+j*Nz, Nz, 1)] = x[std::slice(j*Nx+k, Nz, Nx*Ny)];
				}
			}
			N1 = Ny;
			N2 = Nz;
			break;
		case Angle::Y: case Angle::ZX:
			for (size_t k=0; k<Ny; k++)
			{
				for (size_t j=0; j<Nz; j++)
				{
					xo[std::slice(k*Nz*Nx+j*Nx, Nx, 1)] = x[std::slice(k*Nx+j*Nx*Ny, Nx, 1)];
				}
			}
			N1 = Nz;
			N2 = Nx;
			break;
	}
	x.resize(0);
	Tomography::pngwrite(ip3dgfilename.substr(0,ip3dgfilename.size()-4)+".png", N1, N2, xo, rot_axis, angles, the_grid.scale(), true, false);
	xo.resize(0);
	
	// generate projection lines
	message("setting up rays");
	real scale_d = the_grid.scale();
	size_t restmp = std::max(N1, N2);
	if ( stretch == "inner" )
	{
		restmp = std::min(N1, N2);
	}
	else if ( stretch == "middle" )
	{
		restmp = std::max(N1, N2);
	}
	else if ( stretch == "outer" )
	{
		restmp = size_t(std::sqrt(real(N1*N1+N2*N2))+0.5);
	}
	if ( !resolution )
	{
		resolution = restmp;
	}
	else
	{
		scale_d *= real(restmp) / resolution;
	}
	std::vector< std::valarray<real> > ipoints(0);
	for (size_t j=0; j<resolution; j++)
	{
		for (size_t i=0; i<resolution; i++)
		{
			std::valarray<real> inode = the_grid.center();
			inode[0] += ( real(i) - real(resolution-1) / 2 ) * scale_d;
			inode[1] += ( real(j) - real(resolution-1) / 2 ) * scale_d;
			ipoints.push_back(inode);
		}
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
	message("projecting "+ntos(rays.size())+" rays");
	SparseMatrix<real> A(0,the_grid.ncells());
	std::valarray<real> b;
	std::valarray<bool> blanks(true, rays.size());
	Tomography::projection(the_grid, rays, blanks, b, A, walkfast, true, false, false);
	
	message("clearing elements no longer needed");
	A.clear();
	rays.clear();
	the_grid.clear();
	
	// save projections
	message("saving '"+output+"'");
	std::valarray<real> otb(real(0), blanks.size());
	otb[blanks] = b;
	Tomography::pngwrite(output, resolution, resolution, otb, rot_axis, angles, scale_d, true, true);
	
	return 0;
}
