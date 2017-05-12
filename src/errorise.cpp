/*
  create errors in an image
 
  Rado Faletic
  21st April 2005
*/
 
/*
#undef DEBUG
*/

#include <iostream>
#include <sstream>
#include <string>
#include <valarray>
#include "angles.h"
#include "gaussian.h"
#include "tomography.h"

typedef double real;
                                                                                
int main(int argc, char* argv[])
{
  std::string ifilename = "input.png";
  std::string ofilename = "";
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
		    << "--output=<filename>\n\tthe output filename (determined from input)\n"
		    << "--error=<%>\n\tvariance applied to the pixels (" << error << "%)\n"
		    << "--mode=full|trans\n\terror the full image, or non-transparent (full)" << std::endl;
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
      else if ( input.str().substr(0,8) == std::string("--error=") )
        {
	  std::istringstream choice(input.str().substr(8,input.str().size()-8));
	  choice >> error;
        }
      else if ( input.str().substr(0,7) == std::string("--mode=") )
        {
	  fullmode = ( input.str().substr(7,input.str().size()-7) == "trans" ) ? false : true;
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
      ofilename = ifilename.substr(0, ifilename.size()-4)+"_errorise.png";
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
  bool tpcy = false;
  std::cout << "reading " << ifilename << std::endl;
  Tomography::pngread(ifilename, Nrows, Ncols, data, blanks, raxis, angles, scale, realdata);

  // create the errors
  for (size_t t=0; t<angles.size(); t++)
    {
      std::valarray<real> tdata = data[std::slice(t*Nrows*Ncols,Nrows*Ncols,1)];
      real dave = real(tdata.sum()) / real(tdata.size());
      real variance = error / real(100);
      for (size_t i=0; i<Nrows*Ncols; i++)
	{
	  if ( blanks[i] || fullmode )
	    {
	      data[t*Nrows*Ncols+i] += dave * gaussian(variance);
	      if ( data[t*Nrows*Ncols+i] < real(0) ) data[t*Nrows*Ncols+i] = real(0);
	    }
	  else
	    {
	      tpcy = true;
	    }
	}
    }

  Tomography::pngwrite(ofilename, Nrows, Ncols, data, raxis, angles, scale, realdata, tpcy);
 
  return 0;
}
