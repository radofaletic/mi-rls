/**
 A programme to collate a series of PNG files
 
 Rado Faletic
 15th June 2004
 22nd April 2022
 */

/*
 #undef DEBUG
 */

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <valarray>
#include <vector>
#include "angles.h"
#include "extra_math.h"
#include "tomography.h"

int main(int argc, char* argv[])
{
	std::string ifilename = "input";
	std::string ofilename = "";
    std::size_t trans = 0;
	bool tpcy = false;
    std::size_t xchop = 0;
	bool xdir = true;
	Angle::axes rot_axis = Angle::X;
	
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
			<< "--input=<filename>\n\tinput the filename part of filename_xxx.png (" << ifilename << ")\n"
			<< "--output=<filename>\n\toutput filename (" << ifilename << "_proj.png)\n"
			<< "--transparency=<n>\n\tvalue for transparency (none)\n"
			<< "--xchop=<n>\n\tchop the x resolution (" << xchop << ")\n"
			<< "--xdir=left|right\n\tthe x chop direction (right)\n"
			<< "--axis=<axis>\n\taxis of rotation (X)" << std::endl;
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
		else if ( input.str().substr(0,15) == std::string("--transparency=") )
		{
			std::istringstream choice(input.str().substr(15,input.str().size()-15));
			choice >> trans;
			tpcy = true;
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
		else if ( input.str().substr(0,7) == std::string("--axis=") )
		{
			if ( input.str().substr(7,input.str().size()-7) == std::string("X") ||
				input.str().substr(7,input.str().size()-7) == std::string("YZ") )
			{
				rot_axis = Angle::X;
			}
			else if ( input.str().substr(7,input.str().size()-7) == std::string("Y") ||
					 input.str().substr(7,input.str().size()-7) == std::string("ZX") )
			{
				rot_axis = Angle::Y;
			}
			//else if ( input.str().substr(7,input.str().size()-7) == std::string("Z") ||
			//	    input.str().substr(7,input.str().size()-7) == std::string("XY") )
			//  {
			//    rot_axis = Angle::Z;
			//  }
			else
			{
				std::cerr << "\a\"axis\" must be one of X or Y.\nUse the --help option to learn more.";
				return 1;
			}
		}
		else
		{
			std::cerr << "\aUnrecognised option '"+input.str()+"'.\nUse the --help option to learn more.";
			return 1;
		}
	}
	
	// search for PNG files
    std::size_t Nrows = 0;
    std::size_t Ncols = 0;
	std::valarray<double> angles;
	std::vector< std::valarray<double> > inputs;
    double scale = 1;
	bool realdata = false;
	for (std::size_t i=0; i<180; i++)
	{
		std::string thisfilename = ifilename + "_";
		if ( i < 100 )
		{
			thisfilename += "0";
		}
		if ( i < 10 )
		{
			thisfilename += "0";
		}
		std::ostringstream tmp;
		tmp << i;
		thisfilename += tmp.str() + ".png";
		
		std::ifstream tryfile(thisfilename.c_str());
		if ( !tryfile )
		{
			continue;
		}
		else
		{
			tryfile.close();
		}
		std::cout << "reading " << thisfilename;
        std::size_t tNrows, tNcols, rNrows, rNcols;
		std::valarray<double> tinput;
		Angle::axes taxis;
		std::valarray<double> tangle;
        double tscale;
		std::valarray<bool> blanks;
		if ( Tomography::pngread(thisfilename, tNrows, tNcols, tinput, blanks, taxis, tangle, tscale, realdata) )
		{
			if ( !Nrows && !Ncols )
			{
				Nrows = tNrows;
				Ncols = tNcols;
				rNrows = tNrows;
				rNcols = tNcols;
			}
			else if ( rNrows != tNrows || rNcols != tNcols )
			{
				std::cerr << "\a\nimage dimensions differ" << std::endl;
				return 1;
			}
			if ( xchop )
			{
				Ncols = xchop;
				std::valarray<double> otinput = tinput;
				tinput.resize(Ncols*Nrows);
				for (std::size_t j=0; j<Nrows; j++)
				{
					if ( xdir )
					{
						tinput[std::slice(j*Ncols,Ncols,1)] = otinput[std::slice(j*tNcols+(tNcols-xchop),Ncols,1)];
					}
					else
					{
						tinput[std::slice(j*Ncols,Ncols,1)] = otinput[std::slice(j*tNcols,Ncols,1)];
					}
				}
			}
			inputs.push_back(tinput);
			std::cout << ".";
			if ( tangle.size() > 1 )
			{
				std::cerr << "\a\nstored angles are incompatible" << std::endl;
				return 1;
			}
			tangle.resize(angles.size());
			tangle = angles;
			angles.resize(angles.size()+1);
			angles[std::slice(0,tangle.size(),1)] = tangle;
			angles[angles.size()-1] = i;
			std::cout << ".";
			scale = tscale;
			std::cout << ".";
		}
		std::cout << " done." << std::endl;
	}
	
	// write projection file
	if ( !ofilename.size() )
	{
		ofilename = ifilename + "_proj.png";
	}
	std::valarray<double> output(inputs.size()*Nrows*Ncols);
	for (std::size_t i=0; i<inputs.size(); i++)
	{
		output[std::slice(i*Nrows*Ncols, Nrows*Ncols, 1)] = inputs[i];
	}
	if ( tpcy )
	{
		output += 1.0;
        double tv = trans + 1.0;
		for (std::size_t i=0; i<output.size(); i++)
		{
			if ( eq(output[i], tv) )
			{
				output[i] = 0.0;
			}
		}
	}
	if ( !Tomography::pngwrite(ofilename, Nrows, Ncols, output, rot_axis, angles, scale, realdata, tpcy) )
	{
		std::cerr << "\a\nerror writing " << ofilename << std::endl;
		return 1;
	}
	std::cout << "completed creating " << ofilename << std::endl;
	
	return 0;
}
