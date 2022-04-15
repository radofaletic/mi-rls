/*
 A programme to convert any CFDRC-created Plot3D file (grid, scalar,
 and vector) into a simple Plot3D file format. This consists of a
 Binary file with no I-blanking and no multidomains wherever possible.
 
 Rado Faletic
 12th February 2004
 */

#include <string>
#include <valarray>
#include <vector>
#include "argv.h"
#include "file.h"
#include "front-end.h"
#include "plot3d.h"

int main(int argc, char* argv[])
{
	
	std::string gfilename = "";
	std::valarray<double> X, Y, Z;
	std::valarray<bool> B;
	std::valarray<size_t> Nx, Ny, Nz;
	dataformat format = Formatted;
	dataprecision precision = Single;
	bool multidomain = true;
	bool blanking = true;
	std::string nfilename = "output";
	std::string ext = "";
	std::string dext = "";
	
	std::vector<args> fswitch = get_args(argc, argv);
	if ( fswitch.size() == 0 )
	{
		fswitch.resize(1);
		fswitch[0].var() = "help";
		fswitch[0].val() = "";
	}
	for (size_t i=0; i<fswitch.size(); i++)
	{
		if ( fswitch[i].var("help") )
		{
			message("by Rado Faletic  2003, 2004\n");
			message("this programme converts any CFDRC-created Plot3D file (grid, scalar, and vector) into a simple Plot3D file format. This consists of a Binary file with no I-blanking and no multidomains wherever possible.\n");
			message("below is a list of flags:");
			message("--input=<inputfilename>\n\tthe CFDRC Plot3D grid file (full file name)");
			message("--output=<outputfilename>\n\toutput filename (without extensions)");
		}
		else if ( fswitch[i].var("input") || fswitch[i].var("file") )
		{
			gfilename = fswitch[i].val();
		}
		else if ( fswitch[i].var("output") )
		{
			nfilename = fswitch[i].val();
		}
		else
		{
			message("unrecognised option '"+fswitch[i].var()+"'\nuse the --help option to learn more");
			return 1;
		}
	}
	
	if ( gfilename.size() == 0 )
	{
		return 1;
	}
	if ( gfilename[gfilename.size()-1] == 'D' )
	{
		precision = Double;
		dext = "D";
	}
	switch (precision)
	{
		case Single:
			if ( gfilename[gfilename.size()-2] == 'F' )
			{
				format = Formatted;
				ext += "F";
			}
			else if ( gfilename[gfilename.size()-2] == 'U' )
			{
				format = Unformatted;
				ext += "U";
			}
			else if ( gfilename[gfilename.size()-2] == 'B' )
			{
				format = Binary;
				ext += "B";
			}
			break;
		case Double:
			if ( gfilename[gfilename.size()-3] == 'F' )
			{
				format = Formatted;
				ext += "F";
			}
			else if ( gfilename[gfilename.size()-3] == 'U' )
			{
				format = Unformatted;
				ext += "U";
			}
			else if ( gfilename[gfilename.size()-3] == 'B' )
			{
				format = Binary;
				ext += "B";
			}
			break;
	}
	
	
	Plot3D::read(X, Y, Z, B, Nx, Ny, Nz, gfilename, format, precision, multidomain, blanking);
	
	switch (precision)
	{
		case Single:
			gfilename.resize(gfilename.size()-2);
			break;
		case Double:
			gfilename.resize(gfilename.size()-3);
			break;
	}
	
	std::vector< std::valarray<double> > data(10);
	Plot3D::read_data(data[0], data[1], data[2], data[3], data[4], gfilename+ext+"S"+dext, format, precision, multidomain);
	Plot3D::read_data(data[5], data[6], data[7], data[8], data[9], gfilename+ext+"V"+dext, format, precision, multidomain);
	
	// create a simple Plot3D file that MayaVi can read
	format = Binary;
	precision = Single;
	if ( Nx.size() > 2 )
	{
		multidomain = true;
	}
	else
	{
		multidomain = false;
	}
	blanking = false;
	Plot3D::write(X, Y, Z, B, Nx, Ny, Nz, nfilename+".PBG", format, precision, multidomain, blanking);
	Plot3D::write_data(data[5], std::valarray<double>(data[5]*data[6]), std::valarray<double>(data[5]*data[7]), std::valarray<double>(data[5]*data[8]), data[1], Nx, Ny, Nz, nfilename+".PBS", format, precision, multidomain);
	Plot3D::write_var(nfilename+".VAR", "rho", "um", "vm", "wm", "T");
	
	return 0;
}
