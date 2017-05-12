/*
  A programme to convert a PNG tomographic file into a Plot3D file
 
  Rado Faletic
  7th December 2004
*/
 
/*
#undef DEBUG
*/

#include <iostream>
#include <string>
#include <valarray>
#include "angles.h"
#include "plot3d.h"
#include "tomography.h"

typedef double real;
                                                                                
int main(int argc, char* argv[])
{
  std::string ifilename = "input.png";
  std::string ofilename = "";
  real xw = real(0);
  bool xwr = false;
  real yw = real(0);
  bool ywr = false;
  real l = real(1);
  real gd = real(1);

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
		    << "--output=<filename>\n\tthe output filename with any extensions (determined from input)" << "\n"
		    << "--l=<wavelength>\n\twavelength of the rays (" << l << ")" << "\n"
		    << "--gd=<number>\n\tGladstone-Dale coefficient (" << gd << ")" << "\n"
		    << "--xw=<x-width>\n\tthe x-width of the output file (NULL)" << "\n"
		    << "--yw=<y-width>\n\tthe y-width of the output file (NULL)" << std::endl;
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
      else if ( input.str().substr(0,4) == std::string("--l=") )
        {
          std::istringstream choice(input.str().substr(4,input.str().size()-4));
          choice >> l;
        }
      else if ( input.str().substr(0,5) == std::string("--gd=") )
        {
          std::istringstream choice(input.str().substr(5,input.str().size()-5));
          choice >> gd;
        }
      else if ( input.str().substr(0,5) == std::string("--xw=") )
        {
          std::istringstream choice(input.str().substr(5,input.str().size()-5));
          choice >> xw;
	  xwr = true;
        }
      else if ( input.str().substr(0,5) == std::string("--yw=") )
        {
          std::istringstream choice(input.str().substr(5,input.str().size()-5));
          choice >> yw;
	  ywr = true;
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
      ofilename = ifilename.substr(0, ifilename.size()-4);
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

  if ( data.size() == Nrows * Ncols )
    {
      raxis = Angle::Z;
    }

  if ( l != real(1) || gd != real(1) )
    {
      data *= l / gd;
    }

  // make the grid
  size_t Nx = 1;
  size_t Ny = 1;
  size_t Nz = 1;

  switch(raxis)
    {
    case Angle::X: case Angle::YZ:
      Nx = data.size()/(Nrows*Ncols);
      Ny = Nrows;
      Nz = Ncols;
      break;
    case Angle::Y: case Angle::ZX:
      Nx = Ncols;
      Ny = data.size()/(Nrows*Ncols);
      Nz = Nrows;
      break;
    case Angle::Z: case Angle::XY:
      Nx = Ncols;
      Ny = Nrows;
      Nz = 1;
      break;
    }

  std::valarray<real> X(scale, Nx*Ny*Nz);
  std::valarray<real> Y(scale, Nx*Ny*Nz);
  std::valarray<real> Z(scale, Nx*Ny*Nz);
  std::valarray<bool> B(true, Nx*Ny*Nz);
  std::valarray<real> D(real(0), Nx*Ny*Nz);
  if ( xwr )
    {
      X = ( Nx <= 1 ) ? 1 : xw / ( Nx - 1);
      Y = X[0];
      Z = X[0];
    }
  if ( ywr )
    {
      Y = ( Ny <= 1 ) ? 1 : yw / ( Ny - 1 );
      Z = ( Nz <= 1 ) ? 1 : yw / ( Nz - 1 );
      if ( !xwr ) X = Y[0];
    }
  for (size_t k=0; k<Nz; k++)
    {
      for (size_t j=0; j<Ny; j++)
	{
	  for (size_t i=0; i<Nx; i++)
	    {
	      size_t pos = k*(Ny*Nx) + j*Nx + i;
	      X[pos] *= i;
	      Y[pos] *= j;
	      Z[pos] *= k;
	      if ( !blanks[pos] ) B[pos] = false;
	      switch(raxis)
		{
		case Angle::X: case Angle::YZ:
		  D[pos] = data[k+j*Ncols+i*Nrows*Ncols];
		  break;
		case Angle::Y: case Angle::ZX:
		  D[pos] = data[k*Ncols+j*Nrows*Ncols+i];
		  break;
		case Angle::Z: case Angle::XY:
		  D[pos] = data[j*Ncols+i];
		  break;
		}
	    }
	}
    }
  std::valarray<real> O(real(0), Nx*Ny*Nz);

  message("saving '"+ofilename+"'");
  Plot3D::write(X, Y, Z, B, Nx, Ny, Nz, ofilename+".PBG", Binary);
  Plot3D::write_data(D, O, O, O, O, Nx, Ny, Nz, ofilename+".PBS", Binary);
  Plot3D::write_var(ofilename+".VAR", "ofilename");

  return 0;
}
