/*
 A programme to convert any structured Plot3D file (XYZ and Q)
 into a simple Plot3D file format. This consists of a Binary file
 with no I-blanking and no multidomains wherever possible.
 
 Rado Faletic
 15th May 2004
 */

#include <string>
#include <valarray>
#include <vector>
#include "argv.h"
#include "file.h"
#include "front-end.h"
#include "plot3d.h"

typedef double real;

int main(int argc, char* argv[])
{
	
	std::string gfilename = "";
	std::string qfilename = "";
	std::valarray<real> X, Y, Z;
	std::valarray<bool> B;
	std::valarray<size_t> Nx, Ny, Nz;
	dataformat format = Formatted;
	dataprecision precision = Single;
	bool multidomain = false;
	bool blanking = false;
	std::string nfilename = "output";
	
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
			message("below is a list of flags:");
			message("--xyzfile=<inputfilename>\n\tthe Plot3D XYZ file ("+gfilename+")");
			message("--qfile=<inputfilename>\n\tthe Plot3D Q file ("+qfilename+"");
			message("--multidomain=on|off\n\tturns multidomain on for <inputfilename> (off)");
			message("--iblanking=on|off\n\tturns iblanking on for <inputfilename> (off)");
			message("--format=Formatted|Unformatted|Binary\n\tthe file format for the input file (Formatted)");
			message("--precision=Single|Double\n\tthe file precision for the input file (Single)");
			message("--output=<outputfilename>\n\toutput filename, without extensions ("+nfilename+")");
			return 1;
		}
		else if ( fswitch[i].var("xyzfile") )
		{
			gfilename = fswitch[i].val();
		}
		else if ( fswitch[i].var("qfile") )
		{
			qfilename = fswitch[i].val();
		}
		else if ( fswitch[i].var("multidomain") )
		{
			if ( fswitch[i].val("on") )
			{
				multidomain = true;
			}
		}
		else if ( fswitch[i].var("iblanking") )
		{
			if ( fswitch[i].val("on") )
			{
				blanking = true;
			}
		}
		else if ( fswitch[i].var("format") )
		{
			if ( fswitch[i].val("Unformatted") || fswitch[i].val("u") )
			{
				format = Unformatted;
			}
			else if ( fswitch[i].val("Binary") || fswitch[i].val("b") )
			{
				format = Binary;
			}
		}
		else if ( fswitch[i].var("precision") )
		{
			if ( fswitch[i].val("Double") || fswitch[i].val("d") )
			{
				precision = Double;
			}
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
	Plot3D::read(X, Y, Z, B, Nx, Ny, Nz, gfilename, format, precision, multidomain, blanking);
	
	std::vector< std::valarray<real> > data(5);
	if ( qfilename.size() > 0 )
	{
		Plot3D::read_data(data[0], data[1], data[2], data[3], data[4], qfilename, format, precision, multidomain);
	}
	
	// create a simple Plot3D that MayaVi can read
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
	if ( qfilename.size() > 0 )
	{
		Plot3D::write_data(data[0], data[1], data[2], data[3], data[4], Nx, Ny, Nz, nfilename+".PBS", format, precision, multidomain);
	}
	
	return 0;
}
