/**
 cuts a specific number of data slices from the PNG file
 
 Rado Faletic
 28th November 2004
 22nd April 2022
 */

/*
 #undef DEBUG
 */

#include <iostream>
#include <sstream>
#include <string>
#include <valarray>

#include "angles.h"
#include "plot3d.h"
#include "tomography.h"

int main(int argc, char* argv[])
{
	std::string ifilename = "input.png";
	std::string ofilename = "";
    std::size_t cutter = 0;
	double cut_min = double(0);
	bool cut_min_tf = false;
	
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
			<< "--output=<filename>\n\tthe output PNG filename (determined from input)\n"
			<< "--cut=<n>\n\tthe number of slices to discard (" << cutter << ")\n"
			<< "--min=<n>\n\tthe filler value to use (min_of_image)" << std::endl;
			return 1;
		}
		else if ( input.str().substr(0, 8) == std::string("--input=") )
		{
			ifilename = input.str().substr(8, input.str().size() - 8);
		}
		else if ( input.str().substr(0, 9) == std::string("--output=") )
		{
			ofilename = input.str().substr(9, input.str().size() - 9);
		}
		else if ( input.str().substr(0, 6) == std::string("--cut=") )
		{
			std::istringstream choice(input.str().substr(6, input.str().size() - 6));
			choice >> cutter;
		}
		else if ( input.str().substr(0, 6) == std::string("--min=") )
		{
			std::istringstream choice(input.str().substr(6, input.str().size() - 6));
			choice >> cut_min;
			cut_min_tf = true;
		}
		else
		{
			std::cerr << "\aUnrecognised option '" + input.str() + "'.\nUse the --help option to learn more.";
			return 1;
		}
	}
	if ( ifilename.substr(ifilename.size() - 4, 4) != ".png" && ifilename.substr(ifilename.size() - 4, 4) != ".PNG" )
	{
		ifilename += ".png";
	}
	if ( ofilename == "" )
	{
		ofilename = ifilename.substr(0, ifilename.size() - 4) + "_cut.png";
	}
	else if ( ofilename.substr(ofilename.size() - 4, 4) != ".png" && ofilename.substr(ofilename.size() - 4, 4) != ".PNG" )
	{
		ofilename += ".png";
	}
	
	// search for PNG files
    std::size_t Nrows = 0;
    std::size_t Ncols = 0;
	std::valarray<double> data;
	std::valarray<bool> blanks;
	Angle::axes raxis = Angle::X;
	std::valarray<double> angles;
    double scale = 1;
	bool realdata = true;
	std::cout << "reading " << ifilename << std::endl;
	Tomography::pngread(ifilename, Nrows, Ncols, data, blanks, raxis, angles, scale, realdata);
	
	std::valarray<double> sdata = data;
	if ( cut_min_tf )
	{
		sdata[std::slice(0,cutter * Nrows * Ncols, 1)] = cut_min;
	}
	else
	{
		sdata.resize(data.size() - cutter * Nrows * Ncols);
		sdata = data[std::slice(cutter * Nrows * Ncols, data.size() - cutter * Nrows * Ncols, 1)];
	}
	data.resize(0);
	
	std::cout << "writing " << ofilename << std::endl;
	Tomography::pngwrite(ofilename, Nrows, Ncols, sdata, raxis, angles, scale, realdata, true);
	
	return 0;
}
