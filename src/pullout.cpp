/*
  Pull out one slice from the series of slices
 
  Rado Faletic
  7th December 2004
*/
 
/*
#undef DEBUG
*/

#include <iostream>
#include <sstream>
#include <string>
#include <valarray>
#include "angles.h"
#include "conversions.h"
#include "plot3d.h"
#include "tomography.h"

typedef double real;
                                                                                
int main(int argc, char* argv[])
{
  std::string ifilename = "input.png";
  std::string ofilename = "";
  std::string direction = "slice";
  size_t slice = 1;
  bool halve = false;

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
		    << "--direction=slice|vslice|hslice\n\tpull out a flat slice, vertical slice, or horizontal slice (" << direction << ")\n"
		    << "--slice=<n>\n\tthe number of the slice to pull out (" << slice << ")" << std::endl;
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
      else if ( input.str().substr(0,12) == std::string("--direction=") )
        {
          direction = input.str().substr(12,input.str().size()-12);
        }
      else if ( input.str().substr(0,8) == std::string("--slice=") )
        {
	  if ( input.str().substr(8,input.str().size()-8) == "halve" )
	    {
	      halve = true;
	    }
	  else
	    {
	      std::istringstream choice(input.str().substr(8,input.str().size()-8));
	      choice >> slice;
	    }
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
      ofilename = ifilename.substr(0, ifilename.size()-4)+"_slice"+ntos(slice)+".png";
    }
  else if ( ofilename.substr(ofilename.size()-4,4) != ".png" && ofilename.substr(ofilename.size()-4,4) != ".PNG" )
    {
      ofilename += ".png";
    }

  // search for PNG files
  size_t Nrows = 0;
  size_t Ncols = 0;
  std::valarray<real> data;
  std::valarray<bool> blanks;
  Angle::axes raxis = Angle::X;
  std::valarray<real> angles;
  real scale = 1;
  bool realdata = true;
  std::cout << "reading " << ifilename << std::endl;
  Tomography::pngread(ifilename, Nrows, Ncols, data, blanks, raxis, angles, scale, realdata);

  std::valarray<real> sdata;
  size_t npics = data.size()/(Nrows*Ncols);
  if ( halve )
    {
      sdata.resize(Nrows*Ncols*npics/2);
      for (size_t i=0; i<npics/2; i++)
	{
	  sdata[std::slice(i*Nrows*Ncols,Nrows*Ncols,1)] = data[std::slice(2*i*Nrows*Ncols,Nrows*Ncols,1)];
	}
    }
  else if ( direction == "vslice" || direction == "hslice" )
    {
      sdata.resize(Nrows*angles.size());

      std::vector< std::valarray<real> > ts(angles.size(), std::valarray<real>(Nrows*Ncols));
      for (size_t i=0; i<ts.size(); i++)
	{
	  ts[i] = data[std::slice(i*Nrows*Ncols,Nrows*Ncols,1)];
	}
      Ncols = angles.size();
    }
  else if ( direction == "hslice" )
    {
      sdata.resize(angles.size()*Ncols);
      Nrows = angles.size();
    }
  else
    {
      sdata.resize(Nrows*Ncols);
      sdata = data[std::slice((slice-1)*Nrows*Ncols,Nrows*Ncols,1)];
    }

  std::cout << "writing " << ofilename << std::endl;
  Tomography::pngwrite(ofilename, Nrows, Ncols, sdata, raxis, angles, scale, true, true);

  return 0;
}
