/*
 A programme to perform a linear tomographic inversion on any greyscale PNG file
 
 Rado Faletic
 10th July 2004
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
	std::string ip3dgfilename = "input.PBG";
	std::string ip3dqfilename = "input.PBS";
	short ip3dqi = 1;
	std::string op3dfilename = "output";
	size_t xresolution = 32;
	size_t yresolution = 32;
	size_t nangles = 12;
	Angle::axes rot_axis = Angle::X;
	interpolation_method iterp = BILINEAR;
	bool write_intermediates = false;
	bool fft_on = true;
	bool inv_on = true;
	size_t iterations = 1000;
	real damping = 0.0;
	real smoothing = 0.0;
	
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
			message("--input_g=<inputfile>\n\tthe Plot3D grid file ("+ip3dgfilename+")");
			message("--input_q=<inputfile>\n\tthe Plot3D Q file ("+ip3dqfilename+")");
			message("--input_i=<n>\n\tthe Plot3D Q variable number between 1 and 5 ("+ntos(ip3dqi)+")");
			message("--method=fft/inv\n\twhich inversion method to use (fft & inv)");
			message("--output=<outputfile>\n\tthe output Plot3D/PNG file, without file extension ("+op3dfilename+")");
			message("--xresolution=<res>\n\tx resolution of the projection plane ("+ntos(xresolution)+")");
			message("--yresolution=<res>\n\ty resolution of the projection plane ("+ntos(yresolution)+")");
			message("--angles=<nangles>\n\tthe number of angles, or projections ("+ntos(nangles)+")");
			message("--iterations=<iterations>\n\tnumber of interations to use for inv ("+ntos(iterations)+")");
			message("--damping=<damping>\n\tmatrix damping value for inv ("+ntos(damping)+")");
			message("--smoothing=<smoothing>\n\tmatrix smoothing value for inv ("+ntos(smoothing)+")");
			message("--axis=<axis>\n\taxis of rotation (X)");
			message("--interpolation=<method>\n\twhich interpolation method to use (bilinear)");
			message("--intermediates=on/off\n\twrite PNG files of intermediate steps (off)");
			return 1;
		}
		else if ( fswitch[i].var("input_g") )
		{
			ip3dgfilename = fswitch[i].val();
		}
		else if ( fswitch[i].var("input_q") )
		{
			ip3dqfilename = fswitch[i].val();
		}
		else if ( fswitch[i].var("input_i") )
		{
			std::istringstream itmp(fswitch[i].val());
			itmp >> ip3dqi;
			if ( ip3dqi < 1 || 5 < ip3dqi )
			{
				message("\"input_i\" must be between 1 and 5.\nUse the --help option to learn more.");
				throw; return 1;
			}
		}
		else if ( fswitch[i].var("output") || fswitch[i].var("o") )
		{
			op3dfilename = fswitch[i].val();
		}
		else if ( fswitch[i].var("xresolution") || fswitch[i].var("xres") )
		{
			std::istringstream choice(fswitch[i].val());
			choice >> xresolution;
		}
		else if ( fswitch[i].var("yresolution") || fswitch[i].var("yres") )
		{
			std::istringstream choice(fswitch[i].val());
			choice >> yresolution;
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
		else if ( fswitch[i].var("axis") )
		{
			if ( fswitch[i].val() == std::string("X") || fswitch[i].val() == std::string("YZ") )
			{
				rot_axis = Angle::X;
			}
			else if ( fswitch[i].val() == std::string("Y") || fswitch[i].val() == std::string("ZX") )
			{
				rot_axis = Angle::Y;
			}
			//else if ( fswitch[i].val() == std::string("Z") || fswitch[i].val() == std::string("XY") )
			//  {
			//    rot_axis = Angle::Z;
			//  }
			else
			{
				message("\"axis\" must be one of X or Y.\nUse the --help option to learn more.");
				throw; return 1;
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
	// read grid
	message("reading the grid");
	grid_input mygrid;
	mygrid.type() = structured;
	mygrid.format() = Binary;
	mygrid.precision() = Single;
	mygrid.multidomain() = false;
	mygrid.blanking() = false;
	mygrid.load_grid() = true;
	mygrid.gridfile() = ip3dgfilename;
	mygrid.datafile() = ip3dqfilename;
	mygrid.qdata() = ip3dqi;
	grid<real> the_grid(mygrid);
	the_grid.read_data(mygrid);
	
	std::string gdn = op3dfilename;
	the_grid.give_dataname(op3dfilename);
	
	std::valarray<real> angles;
	
	// generate projection lines
	message("setting up rays");
	real scale_x = 1.0;
	real scale_y = 1.0;
	std::valarray<real> min2(2);
	std::valarray<real> max2(2);
	std::valarray<real> min3 = the_grid.min();
	std::valarray<real> max3 = the_grid.max();
	switch(rot_axis)
	{
		case Angle::X: case Angle::YZ:
			min2[0] = min3[1];
			min2[1] = min3[2];
			max2[0] = max3[1];
			max2[1] = max3[2];
			scale_x = std::abs(max3[0]-min3[0]) / real(xresolution + 1);
			scale_y = norm(&min2, &max2) / real(yresolution + 1);
			break;
		case Angle::Y: case Angle::ZX:
			min2[0] = min3[2];
			min2[1] = min3[0];
			max2[0] = max3[2];
			max2[1] = max3[0];
			scale_x = norm(&min2, &max2) / real(xresolution + 1);
			scale_y = std::abs(max3[1]-min3[1]) / real(yresolution + 1);
			break;
	}
	min2.resize(0);
	max2.resize(0);
	min3.resize(0);
	max3.resize(0);
	real scale = ( scale_x + scale_y ) / 2.0;
	std::vector< std::valarray<real> > ipoints(xresolution*yresolution, std::valarray<real>(3));
	size_t counter = 0;
	for (size_t j=0; j<yresolution; j++)
	{
		for (size_t i=0; i<xresolution; i++)
		{
			ipoints[counter][0] = ( real(i) - real(xresolution-1) / 2.0 ) * scale_x;
			ipoints[counter][1] = ( real(j) - real(yresolution-1) / 2.0 ) * scale_y;
			ipoints[counter][2] = 0.0;
			ipoints[counter] += the_grid.center();
			counter++;
		}
	}
	
	angles.resize(nangles);
	for (size_t i=0; i<angles.size(); i++)
	{
		angles[i] = ((real)(180*i))/((real)(nangles));
	}
	
	Rotation<real> rot(3);
	rot.set_origin(the_grid.center());
	std::valarray<real> islope(3);
	islope[0] = 0;
	islope[1] = 0;
	islope[2] = 1;
	std::vector< line<real> > rays(0);
	for (size_t i=0; i<angles.size(); i++)
	{
		rot.reset(angles[i], rot_axis);
		std::valarray<real> slope = rot.O(islope);
		for (size_t j=0; j<ipoints.size(); j++)
		{
			line<real> tline(slope, rot(ipoints[j]));
			rays.push_back(tline);
		}
	}
	
	// project rays
	message("projecting rays");
	SparseMatrix<real> A(0,the_grid.ncells());
	std::valarray<real> b;
	std::valarray<bool> blanks(true, rays.size());
	Tomography::projection(the_grid, rays, blanks, b, A, walkfast, true, inv_on, true);
	
	if ( write_intermediates ) // save projections
	{
		std::string o = op3dfilename + "_1-projections.png";
		message("saving '"+o+"'");
		std::valarray<real> otb(real(0), blanks.size());
		otb[blanks] = b;
		Tomography::pngwrite(o, yresolution, xresolution, otb, rot_axis, angles, scale, true, true);
	}
	
	if ( fft_on )
	{
		size_t resolution;
		size_t nslices;
		switch(rot_axis)
		{
			case Angle::X: case::Angle::YZ:
				resolution = yresolution;
				nslices = xresolution;
				break;
			case Angle::Y: case::Angle::ZX:
				resolution = xresolution;
				nslices = yresolution;
				break;
		}
		
		// create series of sinograms
		std::valarray<real> ftb(real(0), blanks.size());
		ftb[blanks] = b;
		std::vector< std::valarray<real> > sinogram(nslices, std::valarray<real>(nangles*resolution));
		for (size_t i=0; i<nslices; i++)
		{
			for (size_t j=0; j<nangles; j++)
			{
				switch(rot_axis)
				{
					case Angle::X: case::Angle::YZ:
						sinogram[i][std::slice(j*resolution,resolution,1)] =
						ftb[std::slice(j*resolution*nslices+i,resolution,nslices)];
						break;
					case Angle::Y: case::Angle::ZX:
						sinogram[i][std::slice(j*resolution,resolution,1)] =
						ftb[std::slice(j*resolution*nslices+i*resolution,resolution,1)];
						break;
				}
			}
		}
		ftb.resize(0);
		if ( write_intermediates ) // save sinograms
		{
			std::string o = op3dfilename + "_2-sinograms.png";
			message("saving '"+o+"'");
			std::valarray<real> osinogram(sinogram.size()*nangles*resolution);
			for (size_t i=0; i<sinogram.size(); i++)
			{
				osinogram[std::slice(i*nangles*resolution, nangles*resolution, 1)] = sinogram[i];
			}
			Tomography::pngwrite(o, nangles, resolution, osinogram, rot_axis, angles, scale, false, true);
		}
		
		// perform 1D Fourier Transform on each tomogram in each sinogram
		message("performing 1D FFT on tomograms");
		fftw_complex* one_d_in = (fftw_complex*)fftw_malloc(resolution * sizeof(fftw_complex));
		fftw_complex* one_d_out = (fftw_complex*)fftw_malloc(resolution * sizeof(fftw_complex));
		fftw_plan p1 = fftw_plan_dft_1d(resolution, one_d_in, one_d_out, FFTW_FORWARD, FFTW_MEASURE);
		std::vector< std::valarray<real> > fft_re(sinogram.size(), std::valarray<real>(nangles*resolution));
		std::vector< std::valarray<real> > fft_im(sinogram.size(), std::valarray<real>(nangles*resolution));
		size_t DC = ( !(resolution%2) ) ? (resolution-2)/2 : (resolution-1)/2;
		real sres = std::sqrt(real(resolution));
		for (size_t slice=0; slice<sinogram.size(); slice++)
		{
			for (size_t j=0; j<nangles; j++)
			{
				for (size_t i=0; i<resolution; i++)
				{
					// shift the data so that negative parts are tagged on the end, rather than the beginning
					one_d_in[i][0] = sinogram[slice][j*resolution+(i+DC)%resolution];
					one_d_in[i][1] = 0;
				}
				fftw_execute(p1);
				for (size_t i=0; i<resolution; i++)
				{
					// be sure to shift the DFT output so that negative frequencies come before the DC (ready for 2D interpolation)
					fft_re[slice][j*resolution+i] = one_d_out[(i+(resolution-DC))%resolution][0]/sres;
					fft_im[slice][j*resolution+i] = one_d_out[(i+(resolution-DC))%resolution][1]/sres;
				}
			}
		}
		fftw_destroy_plan(p1);
		fftw_free(one_d_in);
		fftw_free(one_d_out);
		
		if ( write_intermediates ) // saving 1d Fourier Transforms
		{
			std::string o1fft = op3dfilename + "_3-fft_tomograms.png";
			message("saving '"+o1fft+"'");
			std::valarray<real> fout(2*sinogram.size()*nangles*resolution);
			for (size_t i=0; i<sinogram.size(); i++)
			{
				fout[std::slice(2*i*nangles*resolution,nangles*resolution,1)] = fft_re[i];
				fout[std::slice(2*i*nangles*resolution+nangles*resolution,nangles*resolution,1)] = fft_im[i];
			}
			Tomography::pngwrite(o1fft, nangles, resolution, fout, rot_axis, angles, scale, false, true);
		}
		
		// interpolate in 2D Fourier space
		switch (iterp)
		{
			case NN:
				message("natural neighbour interpolating in 2D Fourier space");
				for (size_t i=0; i<sinogram.size(); i++)
				{
					nn_interpolate(fft_re[i], angles);
					nn_interpolate(fft_im[i], angles);
					switch(rot_axis)
					{
						case Angle::X: case::Angle::YZ:
							rotate_matrix(resolution, resolution, fft_re[i]);
							rotate_matrix(resolution, resolution, fft_im[i]);
							break;
					}
				}
				break;
			case BILINEAR:
				message("bilinear interpolating in 2D Fourier space");
				for (size_t i=0; i<sinogram.size(); i++)
				{
					bilinear_interpolate(fft_re[i], angles);
					bilinear_interpolate(fft_im[i], angles);
					switch(rot_axis)
					{
						case Angle::X: case::Angle::YZ:
							rotate_matrix(resolution, resolution, fft_re[i]);
							rotate_matrix(resolution, resolution, fft_im[i]);
							break;
					}
				}
				break;
		}
		
		if ( write_intermediates ) // saving 2D Fourier space image
		{
			std::string o2fft = op3dfilename + "_4-fft_2d.png";
			message("saving '"+o2fft+"'");
			std::valarray<real> fout(2*sinogram.size()*resolution*resolution);
			for (size_t i=0; i<sinogram.size(); i++)
			{
				fout[std::slice(2*i*resolution*resolution,resolution*resolution,1)] = fft_re[i];
				fout[std::slice(2*i*resolution*resolution+resolution*resolution,resolution*resolution,1)] = fft_im[i];
			}
			Tomography::pngwrite(o2fft, resolution, resolution, fout, rot_axis, angles, scale, false, true);
		}
		
		// inverse 2D Fourier Transform
		message("performing reverse 2D FFT on plane");
		fftw_complex* two_d_in = (fftw_complex*)fftw_malloc(resolution * resolution * sizeof(fftw_complex));
		fftw_complex* two_d_out = (fftw_complex*)fftw_malloc(resolution * resolution * sizeof(fftw_complex));
		fftw_plan p2 = fftw_plan_dft_2d(resolution, resolution, two_d_in, two_d_out, FFTW_BACKWARD, FFTW_MEASURE);
		for (size_t s=0; s<sinogram.size(); s++)
		{
			for (size_t i=0; i<resolution; i++) // NOTE: swap x-y format from Column-major to Row-major for FFTW
			{
				for (size_t j=0; j<resolution; j++)
				{
					// be sure to shift the x,y frequencies so that the negative frequencies are reversed
					two_d_in[j*resolution+i][0] = fft_re[s][((i+DC)%resolution)*resolution+(j+DC)%resolution];
					two_d_in[j*resolution+i][1] = fft_im[s][((i+DC)%resolution)*resolution+(j+DC)%resolution];
				}
			}
			fftw_execute(p2);
			for (size_t i=0; i<resolution; i++) // NOTE: swap x-y format from Row-major back to Column-major
			{
				for (size_t j=0; j<resolution; j++)
				{
					// swap the x-y components so that negative frequencies come before DC
					fft_re[s][j*resolution+i] = two_d_out[((i+(resolution-DC))%resolution)*resolution+(j+(resolution-DC))%resolution][0] / resolution;
					fft_im[s][j*resolution+i] = two_d_out[((i+(resolution-DC))%resolution)*resolution+(j+(resolution-DC))%resolution][1] / resolution;
				}
			}
		}
		fftw_destroy_plan(p2);
		fftw_free(two_d_in);
		fftw_free(two_d_out);
		
		if ( write_intermediates ) // saving reverse images
		{
			std::string orf = op3dfilename + "_5-2d.png";
			message("saving '"+orf+"'");
			std::valarray<real> fout(2*sinogram.size()*resolution*resolution);
			for (size_t i=0; i<sinogram.size(); i++)
			{
				fout[std::slice(2*i*resolution*resolution,resolution*resolution,1)] = fft_re[i];
				fout[std::slice(2*i*resolution*resolution+resolution*resolution,resolution*resolution,1)] = fft_im[i];
			}
			Tomography::pngwrite(orf, resolution, resolution, fout, rot_axis, angles, scale, false, true);
		}
		
		// write output
		message("saving '"+op3dfilename+".png'");
		std::valarray<real> sqrfft(sinogram.size()*resolution*resolution);
		for (size_t i=0; i<sinogram.size(); i++)
		{
			sqrfft[std::slice(i*resolution*resolution,resolution*resolution,1)] =
			std::sqrt(fft_re[i]*fft_re[i] + fft_im[i]*fft_im[i]);
		}
		Tomography::pngwrite(op3dfilename+".png", resolution, resolution, sqrfft, rot_axis, angles, scale, true, true);
	}
	
	
	if ( inv_on )
	{
		if ( write_intermediates ) // saving matrix
		{
			std::string ma = op3dfilename + "_6.matrix";
			message("saving '"+ma+"'");
			A.write(ma);
		}
		
		// add smoothing
		if ( smoothing )
		{
			message("adding smoothing equations");
			AddSmoothing(A, b, the_grid, smoothing);
		}
		the_grid.clear_aux();
		
		// doing matrix inversion
		std::valarray<real> x(real(0), A.cols());
		message("inverting matrix");
		lsqr_input<real> linput(damping,0,0,0,iterations);
		LSQR(A, x, b, linput);
		
		for (size_t i=0; i<x.size(); i++)
		{
			the_grid[i] = x[i];
		}
		
		// write output
		message("saving '"+op3dfilename+"'");
		the_grid.write(op3dfilename, SaveGrid, Binary);
		the_grid.write(op3dfilename, SaveData, Binary);
	}
	
	return 0;
}
