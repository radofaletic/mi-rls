/*
 A programme to convert from IEEE real format image file to PNG image file
 
 Rado Faletic
 14th July 2004
 */

/*
 #undef DEBUG
 */

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <valarray>
#include "angles.h"
#include "file.h"
#include "fortran_io.h"
#include "tomography.h"

typedef double real;

int main(int argc, char* argv[])
{
	std::string ifilename = "input.bin";
	std::string ofilename = "output.png";
	dataprecision precision = Single;
	size_t xresolution = 640;
	size_t yresolution = 480;
	size_t xchop = 0;
	bool xdir = true;
	bool byte_swapping = false;
	bool trans = true;
	bool addflip = false;
	
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
			<< "--input=<filename>\n\tthe input filename (" << ifilename << ")\n"
			<< "--output=<filename>\n\tthe output filename (" << ofilename << ")\n"
			<< "--xres=<n>\n\tnumber of columns (" << xresolution << ")\n"
			<< "--yres=<n>\n\tnumber of rows (" << yresolution << ")\n"
			<< "--xchop=<n>\n\tchop the x resolution (" << xchop << ")\n"
			<< "--xdir=left|right\n\tthe x chop direction (right)\n"
			<< "--yres=<n>\n\tnumber of rows (" << yresolution << ")\n"
			<< "--precision=Single|Double\n\tprecision of the input data (Single)\n"
			<< "--swap=on|off\n\tperform byte swapping (off)\n"
			<< "--trans=on|off\n\twrite transparent pixels (on)\n"
			<< "--addlfip\n\tcreate a mirror image and add it to the original" << std::endl;
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
		else if ( input.str().substr(0,7) == std::string("--xres=") )
		{
			std::istringstream choice(input.str().substr(7,input.str().size()-7));
			choice >> xresolution;
		}
		else if ( input.str().substr(0,7) == std::string("--yres=") )
		{
			std::istringstream choice(input.str().substr(7,input.str().size()-7));
			choice >> yresolution;
		}
		else if ( input.str().substr(0,8) == std::string("--xchop=") )
		{
			std::istringstream choice(input.str().substr(8,input.str().size()-8));
			choice >> xchop;
		}
		else if ( input.str().substr(0,7) == std::string("--xdir=") )
		{
			if ( input.str().substr(7,1) == std::string("l") )
			{
				xdir = false;
			}
			else if ( input.str().substr(7,1) == std::string("r") )
			{
				xdir = true;
			}
		}
		else if ( input.str().substr(0,12) == std::string("--precision=") )
		{
			if ( input.str().substr(12,input.str().size()-12) == std::string("Single") ||
				input.str().substr(12,input.str().size()-12) == std::string("single") )
			{
				precision = Single;
			}
			else if ( input.str().substr(12,input.str().size()-12) == std::string("Double") ||
					 input.str().substr(12,input.str().size()-12) == std::string("double") )
			{
				precision = Double;
			}
			else
			{
				std::cerr << "\a\"precision\" must be one of Single or Double.\nUse the --help option to learn more.";
				return 1;
			}
		}
		else if ( input.str().substr(0,7) == std::string("--swap=") )
		{
			if ( input.str().substr(7,input.str().size()-7) == std::string("on") )
			{
				byte_swapping = true;
			}
			else
			{
				byte_swapping = false;
			}
		}
		else if ( input.str().substr(0,8) == std::string("--trans=") )
		{
			if ( input.str().substr(8,input.str().size()-8) == std::string("off") )
			{
				trans = false;
			}
			else
			{
				trans = true;
			}
		}
		else if ( input.str() == std::string("--addflip") )
		{
			addflip = true;
		}
		else
		{
			std::cerr << "\aUnrecognised option '"+input.str()+"'.\nUse the --help option to learn more.";
			return 1;
		}
	}
	
	std::valarray<real> data(xresolution*yresolution);
	
	std::ifstream file;
	file.open(ifilename.c_str(),std::ios_base::binary);
	Fortran::fread_vector(file, &(data[0]), &(data[data.size()]), byte_swapping, Binary, precision);
	
	if ( xchop )
	{
		size_t oxresolution = xresolution;
		xresolution = xchop;
		std::valarray<real> odata = data;
		data.resize(xresolution*yresolution);
		for (size_t j=0; j<yresolution; j++)
		{
			if ( xdir )
			{
				data[std::slice(j*xresolution,xresolution,1)] =
				odata[std::slice(j*oxresolution+(oxresolution-xchop),xresolution,1)];
			}
			else
			{
				data[std::slice(j*xresolution,xresolution,1)] = odata[std::slice(j*oxresolution,xresolution,1)];
			}
		}
	}
	
	real dmin = 0.0;
	real dmax = data.max();
	for (size_t i=0; i<data.size(); i++)
	{
		if ( -500.0 < data[i] && data[i] < dmin )
		{
			dmin = data[i];
		}
	}
	for (size_t i=0; i<data.size(); i++)
	{
		if ( data[i] < dmin )
		{
			data[i] = dmin;
		}
	}
	bool ec = true; // edge cut
	for (size_t i=0; i<xresolution; i++)
	{
		if ( data[i] != dmin )
		{
			ec = false;
		}
	}
	if ( ec )
	{
		std::valarray<real> odata = data;
		yresolution--;
		data.resize(xresolution*yresolution);
		data[std::slice(0,data.size(),1)] = odata[std::slice(xresolution, data.size(),1)];
	}
	ec = true;
	for (size_t i=data.size()-xresolution; i<data.size(); i++)
	{
		if ( data[i] != dmin )
		{
			ec = false;
		}
	}
	if ( ec )
	{
		std::valarray<real> odata = data;
		yresolution--;
		data.resize(xresolution*yresolution);
		data[std::slice(0,data.size(),1)] = odata[std::slice(0, data.size(),1)];
	}
	ec = true;
	for (size_t i=0; i<data.size(); i+=xresolution)
	{
		if ( data[i] != dmin )
		{
			ec = false;
		}
	}
	if ( ec )
	{
		std::valarray<real> odata = data;
		size_t oxresolution = xresolution;
		xresolution--;
		data.resize(xresolution*yresolution);
		for (size_t j=0; j<yresolution; j++)
		{
			data[std::slice(j*xresolution, xresolution, 1)] = odata[std::slice(j*oxresolution+1,xresolution,1)];
		}
	}
	ec = true;
	for (size_t i=xresolution-1; i<data.size(); i+=xresolution)
	{
		if ( data[i] != dmin )
		{
			ec = false;
		}
	}
	if ( ec )
	{
		std::valarray<real> odata = data;
		size_t oxresolution = xresolution;
		xresolution--;
		data.resize(xresolution*yresolution);
		for (size_t j=0; j<yresolution; j++)
		{
			data[std::slice(j*xresolution, xresolution, 1)] = odata[std::slice(j*oxresolution,xresolution,1)];
		}
	}
	if ( addflip )
	{
		std::valarray<real> odata = data;
		data.resize(2*data.size());
		data[std::slice(odata.size(), odata.size(), 1)] = odata;
		for (size_t j=0; j<yresolution; j++)
		{
			data[std::slice(odata.size()-(j+1)*xresolution, xresolution, 1)] = odata[std::slice(j*xresolution, xresolution, 1)];
		}
		yresolution *= 2;
	}
	
	Tomography::pngwrite(ofilename, yresolution, xresolution, data, Angle::X, std::valarray<real>(0), real(1), true, trans);
	
	return 0;
}
