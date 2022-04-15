/*
 A programme to convert a Real-Imaginary image from x-y to r-theta
 
 Rado Faletic
 23rd October 2004
 */

/*
 #undef DEBUG
 */

#include <complex>
#include <string>
#include <valarray>
#include "angles.h"
#include "argv.h"
#include "front-end.h"
#include "tomography.h"

typedef double real;

int main(int argc, char* argv[])
{
	std::string ipngfilename = "input.png";
	std::string opngfilename = "output.png";
	
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
			message("--input=<inputfile>\n\tthe input PNG sinogram file ("+ipngfilename+")");
			message("--output=<outputfile>\n\tthe output PNG file ("+opngfilename+")");
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
	std::valarray<real> data;
	Angle::axes sino_axis;
	std::valarray<real> angles;
	real scale;
	bool realdata;
	std::valarray<bool> blanks;
	Tomography::pngread(ipngfilename, Nrows, Ncols, data, blanks, sino_axis, angles, scale, realdata);
	
	if ( 2*Nrows*Ncols != data.size() )
	{
		message("'"+ipngfilename+"' is not a 2-image file");
		throw; return 1;
	}
	
	real rmin = real(0);
	real rmax = real(0);
	real tmin = real(0);
	real tmax = real(0);
	for (size_t i=0; i<Nrows*Ncols; i++)
	{
		std::complex<real> d(data[i],data[Nrows*Ncols+i]);
		data[i] = std::abs(d);
		data[Nrows*Ncols+i] = std::arg(d);
		rmin = ( !i ) ? data[i] : std::min(rmin, data[i]);
		rmax = ( !i ) ? data[i] : std::max(rmax, data[i]);
		tmin = ( !i ) ? data[i] : std::min(tmin, data[Nrows*Ncols+i]);
		tmax = ( !i ) ? data[i] : std::max(tmax, data[Nrows*Ncols+i]);
	}
	real m = ( rmax - rmin ) / ( tmax - tmin );
	real c = rmax - m * tmax;
	for (size_t i=0; i<Nrows*Ncols; i++)
	{
		//data[i] = tmin;
		//data[Nrows*Ncols+i] = rmin;
		data[Nrows*Ncols+i] = m * data[Nrows*Ncols+i] + c;
	}
	
	Tomography::pngwrite(opngfilename, Nrows, Ncols, data, sino_axis, angles, scale, realdata, false);
	
	return 0;
}
