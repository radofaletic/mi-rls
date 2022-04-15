/*
 A programme to perform an FFT tomographic inversion on any greyscale PNG file
 
 Rado Faletic
 7th July 2004
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
#include "sparse_matrix.h"
#include "tomography.h"

typedef double real;

int main(int argc, char* argv[])
{
	std::string ipngfilename = "input.png";
	std::string opngfilename = "output.png";
	size_t resolution = 256;
	size_t nangles = 180;
	interpolation_method iterp = BILINEAR;
	bool write_intermediates = false;
	
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
			message("\nby Rado Faletic 2003, 2004\n");
			message("below is a list of flags:\n");
			message("--input=<inputfile>\n\tthe input PNG file ("+ipngfilename+")");
			message("--output=<outputfile>\n\tthe output PNG file ("+opngfilename+")");
			message("--resolution=<res>\n\tthe number of data points in each projection ("+ntos(resolution)+")");
			message("--angles=<nangles>\n\tthe number of angles, or projections ("+ntos(nangles)+")");
			message("--interpolation=<interpolation_method>\n\t(bilinear)");
			message("--intermediates=on/off\n\twrite PNG files of intermediate steps (off)");
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
		else if ( fswitch[i].var("resolution") || fswitch[i].var("res") )
		{
			std::istringstream choice(fswitch[i].val());
			choice >> resolution;
		}
		else if ( fswitch[i].var("angles") )
		{
			std::istringstream choice(fswitch[i].val());
			choice >> nangles;
			if ( !nangles )
			{
				message("\"angles\" must be non-zero.\nUse the --help option to learn more.");
				throw; return 1;
			}
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
	std::valarray<real> input_data;
	Angle::axes sino_axis;
	std::valarray<real> angles;
	real scale;
	bool realdata;
	std::valarray<bool> blanks;
	Tomography::pngread(ipngfilename, Nrows, Ncols, input_data, blanks, sino_axis, angles, scale, realdata);
	if ( Nrows*Ncols != input_data.size() )
	{
		message("'"+ipngfilename+"' seems to be a multi-image tomographic file, this software will only read a single-image file.");
		throw; return 1;
	}
	
	// set up grid
	message("setting up the grid");
	grid_input mygrid;
	mygrid.type() = structured;
	mygrid.load_grid() = false;
	mygrid.g_nX() = Ncols;
	mygrid.g_nY() = Nrows;
	mygrid.g_nZ() = 0;
	mygrid.g_scale() = scale;
	grid<real> the_grid(mygrid);
	
	the_grid.put_adata(input_data);
	std::string gdn = opngfilename.substr(0,opngfilename.size()-4);
	the_grid.give_dataname(gdn);
	
	// generate projection lines
	message("setting up projection rays");
	if ( !resolution )
	{
		resolution = size_t(std::ceil(std::sqrt((real)(Ncols*Ncols+Nrows*Nrows))));
	}
	real dlength = scale * std::sqrt((real)(Ncols*Ncols+Nrows*Nrows));
	scale = dlength / ( resolution + 1 );
	dlength /= 2;
	std::vector< std::valarray<real> > ipoints(resolution+2);
	for (size_t i=0; i<ipoints.size(); i++)
	{
		ipoints[i].resize(2);
		ipoints[i][0] = i * scale - dlength;
		ipoints[i][1] = 0;
		ipoints[i] += the_grid.center();
	}
	ipoints.erase(ipoints.end()-1);
	ipoints.erase(ipoints.begin());
	angles.resize(nangles);
	for (size_t i=0; i<angles.size(); i++)
	{
		angles[i] = ((real)(180*i))/((real)(nangles));
	}
	
	Rotation<real> rot(2);
	rot.set_origin(the_grid.center());
	std::valarray<real> islope(2);
	islope[0] = 0;
	islope[1] = 1;
	std::vector< line<real> > rays(0);
	for (size_t i=0; i<angles.size(); i++)
	{
		rot.reset(angles[i], Angle::XY);
		std::valarray<real> slope = rot.O(islope);
		for (size_t j=0; j<ipoints.size(); j++)
		{
			line<real> tline(slope, rot(ipoints[j]));
			rays.push_back(tline);
		}
	}
	
	// create tomograms
	message("creating tomograms");
	SparseMatrix<real> A(0,the_grid.ncells());
	Tomography::projection(the_grid, rays, blanks, input_data, A, walkfast, true, false, true);
	std::valarray<real> sinogram(real(0), blanks.size());
	sinogram[blanks] = input_data;
	
	if ( write_intermediates ) // save sinogram
	{
		std::string o = opngfilename.substr(0,opngfilename.size()-4);
		o += "_1-tomograms.png";
		message("saving '"+o+"'");
		Tomography::pngwrite(o, nangles, resolution, sinogram, sino_axis, angles, scale, false, true);
	}
	
	// perform 1D Fourier Transform on each tomogram
	message("performing 1D FFT on tomograms");
	fftw_complex* one_d_in = (fftw_complex*)fftw_malloc(resolution * sizeof(fftw_complex));
	fftw_complex* one_d_out = (fftw_complex*)fftw_malloc(resolution * sizeof(fftw_complex));
	fftw_plan p1 = fftw_plan_dft_1d(resolution, one_d_in, one_d_out, FFTW_FORWARD, FFTW_MEASURE);
	std::valarray<real> fft_re(sinogram.size());
	std::valarray<real> fft_im(sinogram.size());
	size_t DC = ( !(resolution%2) ) ? (resolution-2)/2 : (resolution-1)/2;
	real sres = std::sqrt(real(resolution));
	for (size_t j=0; j<nangles; j++)
	{
		for (size_t i=0; i<resolution; i++)
		{
			// shift the data so that negative parts are tagged on the end, rather than the beginning
			one_d_in[i][0] = sinogram[j*resolution+(i+DC)%resolution];
			one_d_in[i][1] = 0;
		}
		fftw_execute(p1);
		for (size_t i=0; i<resolution; i++)
		{
			// be sure to shift the DFT output so that negative frequencies come before the DC (ready for 2D interpolation)
			fft_re[j*resolution+i] = one_d_out[(i+(resolution-DC))%resolution][0]/sres;
			fft_im[j*resolution+i] = one_d_out[(i+(resolution-DC))%resolution][1]/sres;
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
		std::valarray<real> fout(fft_re.size()+fft_im.size());
		fout[std::slice(0,fft_re.size(),1)] = fft_re;
		fout[std::slice(fft_re.size(),fft_im.size(),1)] = fft_im;
		Tomography::pngwrite(o1fft, nangles, resolution, fout, sino_axis, angles, scale, false, true);
	}
	
	// interpolate in 2D Fourier space
	switch (iterp)
	{
		case NN:
			message("natural neighbour interpolating in 2D Fourier space");
			nn_interpolate(fft_re, angles);
			nn_interpolate(fft_im, angles);
			break;
		case BILINEAR:
			message("bilinear interpolating in 2D Fourier space");
			bilinear_interpolate(fft_re, angles);
			bilinear_interpolate(fft_im, angles);
			break;
	}
	
	if ( write_intermediates ) // saving 2D Fourier space image
	{
		std::string o2fft = opngfilename.substr(0,opngfilename.size()-4);
		o2fft += "_3-fft_2d.png";
		message("saving '"+o2fft+"'");
		std::valarray<real> fout(fft_re.size()+fft_im.size());
		fout[std::slice(0,fft_re.size(),1)] = fft_re;
		fout[std::slice(fft_re.size(),fft_im.size(),1)] = fft_im;
		Tomography::pngwrite(o2fft, resolution, resolution, fout, sino_axis, angles, scale, false, true);
	}
	
	// inverse 2D Fourier Transform
	message("performing reverse 2D FFT on plane");
	fftw_complex* two_d_in = (fftw_complex*)fftw_malloc(resolution * resolution * sizeof(fftw_complex));
	fftw_complex* two_d_out = (fftw_complex*)fftw_malloc(resolution * resolution * sizeof(fftw_complex));
	fftw_plan p2 = fftw_plan_dft_2d(resolution, resolution, two_d_in, two_d_out, FFTW_BACKWARD, FFTW_MEASURE);
	for (size_t i=0; i<resolution; i++) // NOTE: swap x-y format from Column-major to Row-major for FFTW
	{
		for (size_t j=0; j<resolution; j++)
		{
			// be sure to shift the x,y frequencies so that the negative frequencies are reversed
			two_d_in[j*resolution+i][0] = fft_re[((i+DC)%resolution)*resolution+(j+DC)%resolution];
			two_d_in[j*resolution+i][1] = fft_im[((i+DC)%resolution)*resolution+(j+DC)%resolution];
		}
	}
	fftw_execute(p2);
	for (size_t i=0; i<resolution; i++) // NOTE: swap x-y format from Row-major back to Column-major
	{
		for (size_t j=0; j<resolution; j++)
		{
			// swap the x-y components so that negative frequencies come before DC
			fft_re[j*resolution+i] = two_d_out[((i+(resolution-DC))%resolution)*resolution+(j+(resolution-DC))%resolution][0] / resolution;
			fft_im[j*resolution+i] = two_d_out[((i+(resolution-DC))%resolution)*resolution+(j+(resolution-DC))%resolution][1] / resolution;
		}
	}
	fftw_destroy_plan(p2);
	fftw_free(two_d_in);
	fftw_free(two_d_out);
	
	if ( write_intermediates ) // saving reverse images
	{
		std::string orf = opngfilename.substr(0,opngfilename.size()-4);
		orf += "_4-2d.png";
		message("saving '"+orf+"'");
		std::valarray<real> fout(fft_re.size()+fft_im.size());
		fout[std::slice(0,fft_re.size(),1)] = fft_re;
		fout[std::slice(fft_re.size(),fft_im.size(),1)] = fft_im;
		Tomography::pngwrite(orf, resolution, resolution, fout, sino_axis, angles, scale, false, true);
	}
	
	// write output
	message("saving '"+opngfilename+"'");
	std::valarray<real> sqrfft(resolution*resolution);
	sqrfft = std::sqrt(fft_re*fft_re + fft_im*fft_im);
	Tomography::pngwrite(opngfilename, resolution, resolution, sqrfft, sino_axis, angles, scale, true, true);
	
	return 0;
}
