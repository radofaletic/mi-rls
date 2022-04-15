/*
 A programme to convert any 2D structured Plot3D file (XYZ and Q)
 into a PNG image file.
 
 Rado Faletic
 6th July 2004
 */

#include <sstream>
#include <string>
#include <valarray>
#include "argv.h"
#include "file.h"
#include "front-end.h"
#include "plot3d.h"

int main(int argc, char* argv[])
{
	
	std::string gfilename = "";
	std::string dfilename = "";
	std::valarray<double> X, Y, Z;
	std::valarray<bool> B;
	std::valarray<size_t> Nx, Ny, Nz;
	dataformat format = Binary;
	dataprecision precision = Single;
	bool multidomain = false;
	bool blanking = false;
	short qnum = 1;
	Plot3D::values get_switch = Plot3D::NODES;
	std::string pngfilename = "output.png";
	
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
			message("--xyzfile=<inputfilename>\n\tthe Plot3D XYZ file");
			message("--qfile=<inputfilename>\n\tthe Plot3D Q file");
			message("--multidomain=on|off\n\tturns multidomain on for <inputfilename> (off)");
			message("--iblanking=on|off\n\tturns iblanking on for <inputfilename> (off)");
			message("--format=Formatted|Unformatted|Binary\n\tthe file format for the input file (Binary)");
			message("--precision=Single|Double\n\tthe file precision for the input file (Single)");
			message("--qnumber=<n>\n\tn=1..5 denotes the variable in the qfile to be displayed in the image (1)");
			message("--method=NODES|CENTROIDS\n\timage the value of the Plot3D nodes or the cell centres (NODES)");
			message("--output=<outputfilename>\n\toutput filename (output.png)");
			return 1;
		}
		else if ( fswitch[i].var("xyzfile") )
		{
			gfilename = fswitch[i].val();
		}
		else if ( fswitch[i].var("qfile") )
		{
			dfilename = fswitch[i].val();
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
			if ( fswitch[i].val("Formatted") || fswitch[i].val("f") )
			{
				format = Formatted;
			}
			else if ( fswitch[i].val("Unformatted") || fswitch[i].val("u") )
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
		else if ( fswitch[i].var("qnumber") )
		{
			std::istringstream choice(fswitch[i].val());
			choice >> qnum;
		}
		else if ( fswitch[i].var("method") )
		{
			if ( fswitch[i].val("CENTROIDS") || fswitch[i].val("centroids") || fswitch[i].val("c") )
			{
				get_switch = Plot3D::CENTROIDS;
			}
		}
		else if ( fswitch[i].var("output") )
		{
			pngfilename = fswitch[i].val();
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
	
	std::valarray<double> data;
	if ( dfilename.size() > 0 )
	{
		Plot3D::extract_data(data, dfilename, format, precision, multidomain, qnum);
	}
	
	Plot3D::to_PNG(X, Y, Z, Nx, Ny, Nz, data, pngfilename, get_switch);
	
	return 0;
}
