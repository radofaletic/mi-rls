/*
  Generate statistics about a PNG file
 
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
#include "tomography.h"

typedef double real;
                                                                                
int main(int argc, char* argv[])
{
  std::string ifilename = "input.png";

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
		    << "--input=<filename>\n\tthe input PNG filename (" << ifilename << ")" << std::endl;
          return 1;
        }
      else if ( input.str().substr(0,8) == std::string("--input=") )
        {
          ifilename = input.str().substr(8,input.str().size()-8);
        }
      else if ( input.str().substr(0,2) != std::string("--") )
	{
	  ifilename = input.str();
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

  std::valarray<real> datab = data[blanks];
  std::cout << " SIZE: " << data.size() << std::endl
	    << "   NI: " << data.size()/(Nrows*Ncols) << std::endl
	    << "SCALE: " << scale << std::endl
	    << "  MAX: " << data.max() << std::endl
	    << "  MIN: " << data.min() << std::endl
	    << " MINB: " << datab.min() << std::endl
	    << "  AVE: " << data.sum()/data.size() << std::endl;

  size_t counter = 0;
  for (size_t i=0; i<data.size(); i++)
    {
      if ( data[i] < real(0) )
	{
	  counter++;
	}
    }
  std::cout << "  NEG: " << real(100)*real(counter)/data.size() << "%" << std::endl;

  return 0;
}
