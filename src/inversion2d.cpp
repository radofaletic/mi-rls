/**
 A programme to perform a tomographic inversion of any greyscale PNG file
 
 Rado Faletic
 10th March 2005
 22nd April 2022
 */

/*
 #undef DEBUG
 */

#include <cmath>
#include <fftw3.h>
#include <string>
#include <valarray>
#include <vector>
#include "angles.h"
#include "argv.h"
#include "front-end.h"
#include "grid.h"
#include "interpolation.h"
#include "line.h"
#include "lsqr.h"
#include "matrix_utilities.h"
#include "sparse_matrix.h"
#include "tomography.h"

typedef double real;

int main(int argc, char* argv[])
{
	std::string ipngfilename = "input.png";
	std::string opngfilename = "output.png";
	std::string matrixfile = "";
	std::size_t my_ires = 0;
	interpolation_method iterp = BILINEAR;
	bool write_intermediates = false;
	bool fft_on = true;
	bool inv_on = true;
    std::size_t iterations = 1000;
    double damping = 0.0;
    double smoothing = 1.0;
	std::string stretch = "outer";
	
	std::vector<args> fswitch = get_args(argc, argv);
	if ( !fswitch.size() )
	{
		fswitch.resize(1);
		fswitch[0].var() = "help";
		fswitch[0].val() = "";
	}
	for (std::size_t i=0; i<fswitch.size(); i++)
	{
		if ( fswitch[i].var("help") )
		{
			message("\nby Rado Faletic 2004\n");
			message("below is a list of flags:\n");
			message("--input=<inputfile>\n\tthe input PNG sinogram file ("+ipngfilename+")");
			message("--output=<outputfile>\n\tthe output PNG file ("+opngfilename+")");
			message("--matrix=<matrixfile>\n\tread from this file ("+matrixfile+")");
			message("--resolution=<res>\n\tthe resolution to use, by default is determined from input file ("+std::to_string(my_ires)+")");
			message("--interpolation=bilinear/nn\n\tchoose the interpolation method (bilinear)");
			message("--intermediates=on/off\n\twrite PNG files of intermediate steps (off)");
			message("--method=fft/inv\n\twhich inversion method to use (fft & inv)");
			message("--iterations=<iterations>\n\tnumber of interations to use for inv ("+std::to_string(iterations)+")");
			message("--damping=<damping>\n\tmatrix damping value for inv ("+std::to_string(damping)+")");
			message("--smoothing=<smoothing>\n\tmatrix smoothing value for inv ("+std::to_string(smoothing)+")");
			message("--stretch=inner/middle/outer\n\tthe inversion domain to use ("+stretch+")");
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
		else if ( fswitch[i].var("matrix") || fswitch[i].var("m") )
		{
			matrixfile = fswitch[i].val();
		}
		else if ( fswitch[i].var("resolution") || fswitch[i].var("res") )
		{
			std::istringstream choice(fswitch[i].val());
			choice >> my_ires;
		}
		else if ( fswitch[i].var("interpolation") || fswitch[i].var("interpolate") )
		{
			if ( fswitch[i].val() == std::string("NN") || fswitch[i].val() == std::string("nn") )
			{
				iterp = NN;
			}
			else
			{
				iterp = BILINEAR;
			}
		}
		else if ( fswitch[i].var("intermediates") )
		{
			if ( fswitch[i].val() == std::string("on") || fswitch[i].val() == std::string("ON") ||
				fswitch[i].val() == std::string("yes") || fswitch[i].val() == std::string("YES") )
			{
				write_intermediates = true;
			}
			else
			{
				write_intermediates = false;
			}
		}
		else if ( fswitch[i].var("method") )
		{
			if ( fswitch[i].val() == std::string("fft") || fswitch[i].val() == std::string("FFT") )
			{
				fft_on = true;
				inv_on = false;
			}
			else if ( fswitch[i].val() == std::string("inv") || fswitch[i].val() == std::string("INV") )
			{
				fft_on = false;
				inv_on = true;
			}
		}
		else if ( fswitch[i].var("iterations") || fswitch[i].var("iter") )
		{
			std::istringstream choice(fswitch[i].val());
			choice >> iterations;
		}
		else if ( fswitch[i].var("damping") || fswitch[i].var("damp") )
		{
			std::istringstream choice(fswitch[i].val());
			choice >> damping;
		}
		else if ( fswitch[i].var("smoothing") || fswitch[i].var("smth") )
		{
			std::istringstream choice(fswitch[i].val());
			choice >> smoothing;
		}
		else if ( fswitch[i].var("stretch") )
		{
			stretch = fswitch[i].val();
		}
		else
		{
			message("Unrecognised option '"+fswitch[i].var()+"'.\nUse the --help option to learn more.");
			throw;
            return 1;
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
    std::size_t Nrows, Ncols;
	std::valarray<double> sinogram;
	Angle::axes sino_axis;
	std::valarray<double> angles;
    double scale;
	bool realdata;
	std::valarray<bool> blanks;
	Tomography::pngread(ipngfilename, Nrows, Ncols, sinogram, blanks, sino_axis, angles, scale, realdata);
	
    std::size_t resolution = Ncols;
	if ( Nrows*Ncols != sinogram.size() )
	{
		message("'"+ipngfilename+"' seems to be a multi-image tomographic file, this software will only read a single-image file.");
		throw;
        return 1;
	}
	if ( !(angles.size() != Nrows) )
	{
		angles.resize(Nrows);
		for (std::size_t i=0; i<angles.size(); i++)
		{
			angles[i] = (double(180*i)) / (double(Nrows));
		}
	}
    std::size_t nangles = angles.size();
	
	if ( fft_on )
	{
        std::size_t ires = ( my_ires ) ? my_ires : resolution;
		// perform 1D Fourier Transform on each tomogram
		message("performing 1D FFT on tomograms");
		fftw_complex* one_d_in = (fftw_complex*)fftw_malloc(resolution * sizeof(fftw_complex));
		fftw_complex* one_d_out = (fftw_complex*)fftw_malloc(resolution * sizeof(fftw_complex));
		fftw_plan p1 = fftw_plan_dft_1d(resolution, one_d_in, one_d_out, FFTW_FORWARD, FFTW_MEASURE);
		std::valarray<double> fft_re(sinogram.size());
		std::valarray<double> fft_im(sinogram.size());
        std::size_t DC = ( !(resolution%2) ) ? (resolution-2)/2 : (resolution-1)/2;
        double sres = std::sqrt(double(resolution));
		for (std::size_t j=0; j<nangles; j++)
		{
			for (std::size_t i=0; i<resolution; i++)
			{
				// shift the data so that negative parts are tagged on the end, rather than the beginning
				one_d_in[i][0] = sinogram[j*resolution+(i+DC)%resolution];
				one_d_in[i][1] = 0;
			}
			fftw_execute(p1);
			for (std::size_t i=0; i<resolution; i++)
			{
				// be sure to shift the DFT output so that negative frequencies come before the DC (ready for 2D interpolation)
				fft_re[j*resolution+i] = one_d_out[(i+(resolution-DC))%resolution][0];
				fft_im[j*resolution+i] = one_d_out[(i+(resolution-DC))%resolution][1];
			}
		}
		fftw_destroy_plan(p1);
		fftw_free(one_d_in);
		fftw_free(one_d_out);
		
		if ( write_intermediates ) // saving 1d Fourier Transforms
		{
			std::string o1fft = opngfilename.substr(0,opngfilename.size()-4);
			o1fft += "_2-fft_tomograms.png";
			message("saving '"+o1fft+"'");
			std::valarray<double> fout(fft_re.size()+fft_im.size());
			fout[std::slice(0,fft_re.size(),1)] = fft_re;
			fout[std::slice(fft_re.size(),fft_im.size(),1)] = fft_im;
			Tomography::pngwrite(o1fft, nangles, resolution, fout, sino_axis, angles, scale, false, true);
		}
		
		// interpolate in 2D Fourier space
		switch (iterp)
		{
			case NN:
				message("natural neighbour interpolating in 2D Fourier space, "+std::to_string(ires*ires)+" points");
				message("real values");
				nn_interpolate(fft_re, angles, ires, ires);
				message("imaginary values");
				nn_interpolate(fft_im, angles, ires, ires);
				break;
			case BILINEAR:
				message("bilinear interpolating in 2D Fourier space, "+std::to_string(ires*ires)+" points");
				message("real values");
				bilinear_interpolate(fft_re, angles, ires, ires);
				message("imaginary values");
				bilinear_interpolate(fft_im, angles, ires, ires);
				break;
		}
		
		if ( write_intermediates ) // saving 2D Fourier space image
		{
			std::string o2fft = opngfilename.substr(0,opngfilename.size()-4);
			o2fft += "_3-fft_2d.png";
			message("saving '"+o2fft+"'");
			std::valarray<double> fout(fft_re.size()+fft_im.size());
			fout[std::slice(0,fft_re.size(),1)] = fft_re;
			fout[std::slice(fft_re.size(),fft_im.size(),1)] = fft_im;
			Tomography::pngwrite(o2fft, ires, ires, fout, sino_axis, angles, scale, false, true);
		}
		
		// inverse 2D Fourier Transform
		message("performing reverse 2D FFT on plane");
		fftw_complex* two_d_in = (fftw_complex*)fftw_malloc(ires * ires * sizeof(fftw_complex));
		fftw_complex* two_d_out = (fftw_complex*)fftw_malloc(ires * ires * sizeof(fftw_complex));
		fftw_plan p2 = fftw_plan_dft_2d(ires, ires, two_d_in, two_d_out, FFTW_BACKWARD, FFTW_MEASURE);
		DC = ( !(ires%2) ) ? (ires-2)/2 : (ires-1)/2;
		for (std::size_t i=0; i<ires; i++) // NOTE: swap x-y format from Column-major to Row-major for FFTW
		{
			for (std::size_t j=0; j<ires; j++)
			{
				// be sure to shift the x,y frequencies so that the negative frequencies are reversed
				two_d_in[i*ires+j][0] = fft_re[((j+DC)%ires)*ires+(i+DC)%ires];
				two_d_in[i*ires+j][1] = fft_im[((j+DC)%ires)*ires+(i+DC)%ires];
			}
		}
		fftw_execute(p2);
		DC = ires - DC;
		for (std::size_t j=0; j<ires; j++) // NOTE: swap x-y format from Row-major back to Column-major
		{
			for (std::size_t i=0; i<ires; i++)
			{
				// swap the x-y components so that negative frequencies come before DC
				fft_re[j*ires+i] = two_d_out[((i+DC)%ires)*ires+(j+DC)%ires][0] / (ires*ires);
				fft_im[j*ires+i] = two_d_out[((i+DC)%ires)*ires+(j+DC)%ires][1] / (ires*ires);
			}
		}
		fftw_destroy_plan(p2);
		fftw_free(two_d_in);
		fftw_free(two_d_out);
		
		if ( write_intermediates ) // saving reverse images
		{
			std::string orf = opngfilename.substr(0,opngfilename.size()-4);
			orf += "_4-2dfft.png";
			message("saving '"+orf+"'");
			std::valarray<double> fout(fft_re.size()+fft_im.size());
			fout[std::slice(0,fft_re.size(),1)] = fft_re;
			fout[std::slice(fft_re.size(),fft_im.size(),1)] = fft_im;
			Tomography::pngwrite(orf, ires, ires, fout, sino_axis, angles, scale, false, true);
		}
		
		std::valarray<double> sqrfft(ires*ires);
		sqrfft = std::sqrt(fft_re*fft_re + fft_im*fft_im);
		fft_re.resize(0);
		fft_im.resize(0);
		
		if ( write_intermediates ) // saving reverse images
		{
			std::string offtrn = opngfilename.substr(0,opngfilename.size()-4);
			offtrn += "_5-fft.png";
			message("saving '"+offtrn+"'");
			Tomography::pngwrite(offtrn, ires, ires, sqrfft, sino_axis, angles, scale, false, false);
		}
		
        std::size_t fres = std::size_t(double(ires)*std::sqrt(double(0.5))+double(0.5)); // "outer"
		if ( stretch == "inner" || stretch == "middle" )
		{
			fres = ires;
		}
		std::valarray<double> cfft(0.0, fres*fres);
		for (std::size_t j=0; j<fres; j++)
		{
			cfft[std::slice(j*fres,fres,1)] = sqrfft[std::slice(ires*((ires-fres)/2)+((ires-fres)/2)+j*ires,fres,1)];
		}
		sqrfft.resize(0);
		
		// write output
		std::string offtn = opngfilename.substr(0,opngfilename.size()-4);
		offtn += "_fft.png";
		message("saving '"+offtn+"'");
		Tomography::pngwrite(offtn, fres, fres, cfft, sino_axis, angles, scale, true, false);
	}
	
	if ( inv_on )
	{
        std::size_t gres = ( my_ires ) ? my_ires : std::size_t(double(resolution)*std::sqrt(double(0.5))+double(0.5));
		if ( stretch == "inner" || stretch == "middle" )
		{
			gres = ( my_ires ) ? my_ires : resolution;
		}
		
		// set up grid
		message("setting up the grid");
		grid_input mygrid;
		mygrid.type() = structured;
		mygrid.load_grid() = false;
		mygrid.g_nX() = gres;
		mygrid.g_nY() = gres;
		mygrid.g_nZ() = 0;
		mygrid.g_scale() = ( my_ires ) ? scale * double(gres) / my_ires : scale;
		grid<double> the_grid(mygrid);
		
		the_grid.give_dataname(opngfilename.substr(0,opngfilename.size()-4));
		
		SparseMatrix<double> A(0,the_grid.ncells());
		std::valarray<double> b(0);
		
		if ( matrixfile.size() )
		{
			message("reading matrix "+matrixfile);
			A.read(matrixfile, b);
		}
		else
		{
			// generate projection lines
			message("setting up projection rays");
			std::vector< std::valarray<double> > ipoints(0);
			for (std::size_t i=0; i<resolution; i++)
			{
				std::valarray<double> tmp_p = the_grid.center();
				tmp_p[0] += ( double(i) - double(resolution-1) / 2 ) * mygrid.g_scale();
				ipoints.push_back(tmp_p);
			}
			Rotation<double> rot(2);
			rot.set_origin(the_grid.center());
			std::valarray<double> islope(2);
			islope[0] = 0;
			islope[1] = 1;
			std::vector< line<double> > rays(0);
			for (std::size_t i=0; i<angles.size(); i++)
			{
				rot.reset(angles[i], Angle::XY);
				std::valarray<double> slope = rot.O(islope);
				for (std::size_t j=0; j<ipoints.size(); j++)
				{
					line<double> tline(slope, rot(ipoints[j]));
					rays.push_back(tline);
				}
			}
			
			message("projecting "+std::to_string(rays.size())+" rays");
			std::valarray<double> btmp;
			Tomography::projection(the_grid, rays, blanks, btmp, A, walkfast, true, true, true);
			
			if ( true ) // projections
			{
				std::string ma = opngfilename.substr(0,opngfilename.size()-4);
				ma += "_projection.png";
				message("saving '"+ma+"'");
				std::valarray<double> otb(double(0), blanks.size());
				otb[blanks] = btmp;
				Tomography::pngwrite(ma, angles.size(), resolution, otb, sino_axis, angles, scale, true, false);
			}
			
			b.resize(btmp.size());
			b = sinogram[blanks];
			
			std::string oma = opngfilename.substr(0,opngfilename.size()-4);
			oma = opngfilename.substr(0,opngfilename.size()-4);
			oma += "_inv.matrix";
			message("saving '"+oma+"'");
			A.write(oma, b);
		}
		
		// add smoothing
		if ( smoothing )
		{
			message("adding smoothing equations, with smoothing = "+std::to_string(smoothing));
			AddSmoothing(A, b, the_grid, smoothing);
		}
		message("clearing grid");
		the_grid.clear();
		
		// doing matrix inversion
		std::valarray<double> x(double(0), A.cols());
		message("inverting matrix");
		lsqr_input<double> linput(damping,0,0,0,iterations,true);
		lsqr_output<double> Aoutput = LSQR(A, x, b, linput);
		
		// write output
		std::string oinvn = opngfilename.substr(0,opngfilename.size()-4);
		oinvn += "_inv.png";
		message("saving '"+oinvn+"'");
		Tomography::pngwrite(oinvn, gres, gres, x, sino_axis, angles, scale, true, false);
		oinvn = oinvn.substr(0,opngfilename.size()-4);
		oinvn += "_inv.lsqr";
		message("saving '"+oinvn+"'");
		Aoutput.write(oinvn);
		
	}
	
	return 0;
}
