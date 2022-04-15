/*
 adds transparency to a PNG image
 
 Rado Faletic
 14th November 2004
 */

/*
 #undef DEBUG
 */

#include <iostream>
#include <sstream>
#include <string>
#include <valarray>
#include "angles.h"
#include "tomography.h"

typedef double real;

int main(int argc, char* argv[])
{
	std::string ifilename = "input.png";
	std::string ofilename = "";
	real spread = real(75);
	real error = real(10);
	bool fullmode = true;
	
	// read command-line switches
	if ( argc == 1 )
	{
		std::cerr << "\ause the --help option to learn more" << std::endl;
		return 1;
	}
	for (unsigned int i=1; i<argc; i++)
	{
		std::ostringstream input(argv[i]);
		if ( input.str().substr(0,6) == std::string("--help") )
		{
			std::cout << "by Rado Faletic 2004\n"
			<< "below is a list of flags:\n"
			<< "--input=<filename>\n\tthe input PNG filename (" << ifilename << ")\n"
			<< "--output=<filename>\n\tthe output filename (determined from input)" << std::endl;
			return 1;
		}
		else if ( input.str().substr(0,8) == std::string("--input=") )
		{
			ifilename = input.str().substr(8,input.str().size()-8);
		}
		else if ( input.str().substr(0,9) == std::string("--output=") )
		{
			ofilename = input.str().substr(9,input.str().size()-9);
		}
		else
		{
			std::cerr << "\aUnrecognised option '"+input.str()+"'.\nUse the --help option to learn more.";
			return 1;
		}
	}
	if ( ifilename.substr(ifilename.size()-4,4) != ".png" && ifilename.substr(ifilename.size()-4,4) != ".PNG" )
	{
		ifilename += ".png";
	}
	if ( ofilename == "" )
	{
		ofilename = ifilename.substr(0, ifilename.size()-4)+"_trans.png";
	}
	
	size_t Nrows = 0;
	size_t Ncols = 0;
	std::valarray<real> data;
	std::valarray<bool> blanks;
	Angle::axes raxis = Angle::X;
	std::valarray<real> angles;
	real scale = 1;
	bool realdata = true;
	bool tpcy = false;
	std::cout << "reading " << ifilename << std::endl;
	Tomography::pngread(ifilename, Nrows, Ncols, data, blanks, raxis, angles, scale, realdata);
	
	Tomography::pngwrite(ofilename, Nrows, Ncols, data, raxis, angles, scale, realdata, true);
	
	return 0;
}
