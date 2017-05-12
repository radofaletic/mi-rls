/*
  linearly scale real data in a PNG file
 
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
  std::string ofilename = "";
  real scale_mul = real(1);
  real scale_add = real(0);

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
		    << "--mul=<m>\n\tthe linear multiplier (" << scale_mul << ")\n"
		    << "--add=<c>\n\tthe linear additive (" << scale_add << ")" << std::endl;
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
      else if ( input.str().substr(0,6) == std::string("--mul=") )
        {
          std::istringstream choice(input.str().substr(6,input.str().size()-6));
          choice >> scale_mul;
        }
      else if ( input.str().substr(0,6) == std::string("--add=") )
        {
          std::istringstream choice(input.str().substr(6,input.str().size()-6));
          choice >> scale_add;
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
      ofilename = ifilename.substr(0, ifilename.size()-4)+"_scale.png";
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

  if ( !realdata )
    {
      std::cerr << "PNG file '" << ifilename << "' does not contain real data." << std::endl;
      //return 1;
    }

  if ( scale_mul != real(1) )
    {
      data *= scale_mul;
    }
  if ( scale_add != real(0) )
    {
      data += scale_add;
    }

  std::cout << "writing " << ofilename << std::endl;
  Tomography::pngwrite(ofilename, Nrows, Ncols, data, raxis, angles, scale, true, true);

  return 0;
}
