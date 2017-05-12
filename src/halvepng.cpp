/*
  A programme to perform a linear tomographic inversion on any greyscale PNG file
 
  Rado Faletic
  13th July 2004
*/
 
/*
#undef DEBUG
*/
 
#include <cmath>
#include <iostream>
#include <sstream>
#include <string>
#include <valarray>
#include <vector>
#include "angles.h"
#include "tomography.h"
                                                                                
typedef double real;
                                                                                
int main(int argc, char* argv[])
{
  std::string ifilename = "input.png";
  std::string ofilename = "output.png";

  if ( argc != 3 )
    {
      std::cerr << "\ause the --help option to learn more" << std::endl;
      return 1;
    }
  for (unsigned short i=1; i<argc; i++)
    {
      std::ostringstream input(argv[i]);
      if ( input.str().substr(0,2) == std::string("--") )
	{
	  if ( input.str().substr(0,6) == std::string("--help") )
	    {
	      std::cout << "by Rado Faletic 2004\n"
			<< "below is a list of flags:\n"
			<< "--input=<filename>\n\tinput filename (" << ifilename << ")\n"
			<< "--output=<filename>\n\toutput filename (" << ofilename << ")" << std::endl;
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
      else
	{
	  switch(i)
	    {
	    case 1:
	      ifilename = input.str();
	      break;
	    case 2:
	      ofilename = input.str();
	      break;
	    }
	}
    }

  size_t oNrows, oNcols;
  std::valarray<real> opngdata;
  std::valarray<bool> oblanks;
  Angle::axes axis;
  std::valarray<real> angles;
  real scale;
  bool realdata;
  Tomography::pngread(ifilename, oNrows, oNcols, opngdata, oblanks, axis, angles, scale, realdata);
  if ( !angles.size() )
    {
      angles.resize(1, real(0));
    }
  real dmin = opngdata.min();

  size_t Nrows = ( oNrows%2 ) ? ( oNrows / 2 ) + 1 : oNrows / 2;
  size_t Ncols = ( oNcols%2 ) ? ( oNcols / 2 ) + 1 : oNcols / 2;
  std::valarray<real> pngdata(angles.size()*Nrows*Ncols);
  std::valarray<bool> blanks(angles.size()*Nrows*Ncols);

  std::valarray<bool> blank(4);
  std::valarray<real> value(4);
  for (size_t angle=0; angle<angles.size(); angle++)
    {
      size_t oz = angle * oNrows * oNcols;
      size_t z = angle * Nrows * Ncols;
      for (size_t j=0; j<Nrows; j++)
	{
	  size_t oj1 = 2*j;
	  size_t oj2 = ( oj1 == oNrows - 1 ) ? oj1 : oj1 + 1;
	  for (size_t i=0; i<Ncols; i++)
	    {
	      size_t oi1 = 2*i;
	      size_t oi2 = ( oi1 == oNcols - 1 ) ? oi1 : oi1 + 1;
	      blank[0] = oblanks[oz+oj1*oNcols+oi1];
	      blank[1] = oblanks[oz+oj1*oNcols+oi2];
	      blank[2] = oblanks[oz+oj2*oNcols+oi1];
	      blank[3] = oblanks[oz+oj2*oNcols+oi2];
	      value[0] = opngdata[oz+oj1*oNcols+oi1];
	      value[1] = opngdata[oz+oj1*oNcols+oi2];
	      value[2] = opngdata[oz+oj2*oNcols+oi1];
	      value[3] = opngdata[oz+oj2*oNcols+oi2];
	      unsigned short counter = 0;
	      for (size_t b=0; b<4; b++)
		{
		  if ( blank[b] )
		    {
		      counter++;
		    }
		}
	      if ( counter < 2 )
		{
		  pngdata[z+j*Ncols+i] = dmin;
		  blanks[z+j*Ncols+i] = false;
		}
	      else if ( counter == 2 )
		{
		  if ( blank[0] && blank[1] )
		    {
		      pngdata[z+j*Ncols+i] = ( value[0] + value[1] ) / 2;
		    }
		  else if ( blank[0] && blank[2] )
		    {
		      pngdata[z+j*Ncols+i] = ( value[0] + value[2] ) / 2;
		    }
		  else if ( blank[0] && blank[3] )
		    {
		      pngdata[z+j*Ncols+i] = ( value[0] + value[3] ) / 2;
		    }
		  else if ( blank[1] && blank[2] )
		    {
		      pngdata[z+j*Ncols+i] = ( value[1] + value[2] ) / 2;
		    }
		  else if ( blank[1] && blank[3] )
		    {
		      pngdata[z+j*Ncols+i] = ( value[1] + value[3] ) / 2;
		    }
		  else if ( blank[2] && blank[3] )
		    {
		      pngdata[z+j*Ncols+i] = ( value[2] + value[3] ) / 2;
		    }
		  blanks[z+j*Ncols+i] = true;
		}
	      else if ( counter == 3 )
		{
		  if ( blank[0] && blank[1] && blank[2] )
		    {
		      pngdata[z+j*Ncols+i] = ( value[0] + value[1] + value[2] ) / 3;
		    }
		  else if ( blank[0] && blank[1] && blank[3] )
		    {
		      pngdata[z+j*Ncols+i] = ( value[0] + value[1] + value[3] ) / 3;
		    }
		  else if ( blank[0] && blank[2] && blank[3] )
		    {
		      pngdata[z+j*Ncols+i] = ( value[0] + value[2] + value[3] ) / 3;
		    }
		  else if ( blank[1] && blank[2] && blank[3] )
		    {
		      pngdata[z+j*Ncols+i] = ( value[1] + value[2] + value[3] ) / 3;
		    }
		  blanks[z+j*Ncols+i] = true;
		}
	      else
		{
		  pngdata[z+j*Ncols+i] = ( value[0] + value[1] + value[2] + value[3] ) / 4;
		  blanks[z+j*Ncols+i] = true;
		}
	    }
	}
    }
  scale *= 2;
  bool tpcy = false;
  for (size_t i=0; i<blanks.size(); i++)
    {
      if ( !blanks[i] )
	{
	  tpcy = true;
	  break;
	}
    }
  if ( angles.size() == 1 && angles[0] == real(0) )
    {
      angles.resize(0);
    }
  Tomography::pngwrite(ofilename, Nrows, Ncols, pngdata, axis, angles, scale, realdata, tpcy);

  return 0;
}
