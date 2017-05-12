/*
  A programme to return the real component of a complex data set

  Rado Faletic
  27th March 2005
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

  std::valarray<real> d = data[std::slice(0,Nrows*Ncols,1)];

  Tomography::pngwrite(opngfilename, Nrows, Ncols, d, sino_axis, angles, scale, realdata, false);

  return 0;
}
