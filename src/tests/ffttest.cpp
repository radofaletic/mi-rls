/**
 for testing ray tracing AND fourier techniques
 */

#include <cmath>
#include <fftw3.h>
#include <string>
#include <valarray>
#include <vector>

#ifndef USE_MESSAGES
#define USE_MESSAGES
#endif /* USE_MESSAGES */

#include "angles.h"
#include "argv.h"
#include "conversions.h"
#include "extra_math.h"
#include "front-end.h"
#include "grid.h"
#include "imglib.h"
#include "line.h"
#include "nngridr.h"
#include "rotation.h"
#include "sparse_matrix.h"
#include "tomography.h"

#ifdef USE_MESSAGES
#include <fstream>
std::ofstream messages;
#endif /* USE_MESSAGES */

int main(int argc, char* argv[])
{
#ifdef USE_MESSAGES
	messages.open((std::string(argv[0])+std::string(".messages")).c_str());
#endif /* USE_MESSAGES */
	
	std::vector<args> fswitch = get_args(argc, argv);
	
	std::string the_pfilename;
	grid_input gridinputs;
	gridinputs.type() = structured;
	gridinputs.format() = Formatted;
	gridinputs.precision() = Single;
	gridinputs.multidomain() = false;
	gridinputs.blanking() = false;
	gridinputs.gridfile() = "data/phantoms/dorn_3d.PFG";
	gridinputs.datafile() = "data/phantoms/dorn_3d.PFS";
	gridinputs.qdata() = 1;
	int rotations = 32;
	std::string the_dataname = "dorn";
	for (std::size_t i=0; i<fswitch.size(); i++)
	{
		if ( fswitch[i].var("type") )
		{
			if ( fswitch[i].val("structured") || fswitch[i].val("Structured") )
			{
				gridinputs.type() = structured;
			}
		}
		else if ( fswitch[i].var("format") )
		{
			if ( fswitch[i].val("formatted") || fswitch[i].val("Formatted") )
			{
				gridinputs.format() = Formatted;
			}
			else if ( fswitch[i].val("unformatted") || fswitch[i].val("Unformatted") )
			{
				gridinputs.format() = Unformatted;
			}
			else if ( fswitch[i].val("binary") || fswitch[i].val("Binary") )
			{
				gridinputs.format() = Binary;
			}
		}
		else if ( fswitch[i].var("precision") )
		{
			if ( fswitch[i].val("single") || fswitch[i].val("Single") )
			{
				gridinputs.precision() = Single;
			}
			else if ( fswitch[i].val("double") || fswitch[i].val("Double") )
			{
				gridinputs.precision() = Double;
			}
		}
		else if ( fswitch[i].var("multidomain") )
		{
			if ( fswitch[i].val("yes") || fswitch[i].val("true") || fswitch[i].val("1") || fswitch[i].val("on") )
			{
				gridinputs.multidomain() = true;
			}
			else if ( fswitch[i].val("no") || fswitch[i].val("false") || fswitch[i].val("0") || fswitch[i].val("off") )
			{
				gridinputs.multidomain() = false;
			}
		}
		else if ( fswitch[i].var() == "blanking" )
		{
			if ( fswitch[i].val("yes") || fswitch[i].val("true")|| fswitch[i].val("1") || fswitch[i].val("on") )
			{
				gridinputs.blanking() = true;
			}
			else if ( fswitch[i].val("no") || fswitch[i].val("false") || fswitch[i].val("0") || fswitch[i].val("off") )
			{
				gridinputs.blanking() = false;
			}
		}
		else if ( fswitch[i].var("gridfile") )
		{
			gridinputs.gridfile() = fswitch[i].val();
		}
		else if ( fswitch[i].var("datafile") )
		{
			gridinputs.datafile() = fswitch[i].val();
		}
		else if ( fswitch[i].var("qdata") )
		{
			gridinputs.qdata() = std::stol(fswitch[i].val());
		}
		else if ( fswitch[i].var("pngfile") )
		{
			the_pfilename = fswitch[i].val();
			the_dataname = the_pfilename.substr(0,the_pfilename.size()-4);
		}
		else if ( fswitch[i].var("dataname") )
		{
			the_dataname = fswitch[i].val();
		}
		else if ( fswitch[i].var("rotations") )
		{
			rotations = std::stol(fswitch[i].val());
		}
	}
    std::size_t Nrows;
    std::size_t Ncols;
	std::valarray<double> b;
	if ( the_pfilename.size() == 0 )
	{
		/*
		 do test projections
		 */
		
		//
		// read in grid and a-priori data
		//
		
		grid<double> test_grid(gridinputs);
		test_grid.read_data(gridinputs);
		test_grid.give_dataname(the_dataname);
		std::string the_filename = test_grid.get_filename();
		
		test_grid.set_basis(unity);
		
		//
		// a blank projection pane
		//
		
		message("creating generic projection plane");
        std::size_t meshsize = 24;
		//std::size_t meshsize = 128;
		std::valarray<double>* ratt = new std::valarray<double>(test_grid.scale()/10.0f,test_grid.dim());
		std::valarray<double>* ratu = new std::valarray<double>(test_grid.min() + *ratt);
		std::valarray<double>* ratv = new std::valarray<double>(test_grid.max() - *ratt);
		std::vector< line<double> > rays = Tomography::synthetic_rays(rotations, meshsize, *ratu, *ratv, Nrows, Ncols);
		delete ratt;
		delete ratu;
		delete ratv;
		
		//
		// get stuck into the ray tracing
		//
		message("ray tracing");
		b.resize(rays.size());
		std::valarray<bool> blanks;
		SparseMatrix<double> A(0,0);
		Tomography::projection(test_grid, rays, blanks, b, A, walkfast, true, false, true);
		message("\t...done\n");
		
		the_pfilename = the_dataname + ".png";
		message("writing PNG file "+the_pfilename);
		Tomography::pngwrite(the_filename+"-projected.png", Nrows, Ncols, b);
		message("\t...done");
		
		test_grid.clear_data();
		test_grid.clear();
	}
	else
	{
		
	}
	double scale_x = 0;
	double scale_y = 0;
	Angle::axes axis = Angle::Y;
	std::valarray<double> angles;
	std::vector< std::valarray<double> > tpoints;
	message("reading PNG file "+the_pfilename);
	pngread(the_pfilename, the_dataname, Nrows, Ncols, scale_x, scale_y, b, axis, angles);
	std::vector< std::valarray<double> > data(b.size()/(Nrows*Ncols));
	for (std::size_t i=0; i<data.size(); i++)
	{
		data[i].resize(Nrows*Ncols);
		for (std::size_t j=0; j<data[i].size(); j++)
		{
			data[i] = b[std::slice(i*Nrows*Ncols,Nrows*Ncols,1)];
		}
	}
	message("\t...done");
	
	/*
	 Now we do the Fourier inversion:
	 - make the output data array
	 - loop over each slice
	 {
	 - get each strip for this slice, and put them into a new array
	 - perform 1D FFT on each strip
	 - put each of these FFT'd strips onto a 2D plane
	 - interpolate the 2D plane
	 - do a reverse 2D FFT on the plane
	 - save the final data to the output data array
	 }
	 */
	
	// set up Fourier inversion tools
	message("setting up Fourier tools");
	const std::size_t number_of_slices = ( axis == Angle::X ) ? Ncols : Nrows;
	const std::size_t output_dimension = ( axis == Angle::X ) ? Nrows : Ncols;
	const std::size_t strip_start_factor = ( axis == Angle::X ) ? 1 : Ncols;
	const std::size_t strip_stride = ( axis == Angle::X ) ? Ncols : 1;
	std::vector< std::valarray<double> > tpointso;
	std::vector< std::valarray<double> > output_re(number_of_slices);
	std::vector< std::valarray<double> > output_im(number_of_slices);
	std::vector< std::valarray<double> > output_fft_re(number_of_slices);
	std::vector< std::valarray<double> > output_fft_im(number_of_slices);
	for (std::size_t i=0; i<number_of_slices; i++)
	{
		output_re[i].resize(output_dimension*output_dimension);
		output_im[i].resize(output_dimension*output_dimension);
		output_fft_re[i].resize(output_dimension*output_dimension);
		output_fft_im[i].resize(output_dimension*output_dimension);
	}
	fftw_complex* one_d_in = (fftw_complex*)fftw_malloc(output_dimension * sizeof(fftw_complex));
	fftw_complex* one_d_out = (fftw_complex*)fftw_malloc(output_dimension * sizeof(fftw_complex));
	fftw_plan p1 = fftw_plan_dft_1d(output_dimension, one_d_in, one_d_out, FFTW_FORWARD, FFTW_MEASURE);
	fftw_complex* two_d_in = (fftw_complex*)fftw_malloc(output_dimension*output_dimension * sizeof(fftw_complex));
	fftw_complex* two_d_out = (fftw_complex*)fftw_malloc(output_dimension*output_dimension * sizeof(fftw_complex));
	fftw_plan p2 = fftw_plan_dft_2d(output_dimension, output_dimension, two_d_in, two_d_out, FFTW_BACKWARD, FFTW_MEASURE);
	std::vector< fftw_complex* > strips(data.size());
	for (std::size_t i=0; i<strips.size(); i++)
	{
		strips[i] = (fftw_complex*)fftw_malloc(output_dimension * sizeof(fftw_complex));
	}
	double re_max = 0;
	int re_i = 0;
	int re_j = 0;
	message("\t...done");
	
	// loop over each slice
	for (std::size_t slice=0; slice<number_of_slices; slice++)
	{
		message("computing Fourier inversion for slice "+std::to_string(slice+1)+" of "+std::to_string(number_of_slices));
		// get each strip for this slice, and put it into the "strips" array
		debug("setting up strips for slice");
		for (std::size_t i=0; i<strips.size(); i++)
		{
			for (std::size_t j=0; j<output_dimension; j++)
			{
				strips[i][j][0] = data[i][strip_start_factor*slice + strip_stride*j];
				strips[i][j][1] = 0;
			}
		}
		// do a 1D FFT for each strip in "strips"
		debug("computing 1D Fourier transform for each strip");
        std::size_t topz = ( output_dimension%2 == 0 ) ? output_dimension/2 : (output_dimension-1)/2;
        std::size_t DC = ( output_dimension%2 == 0 ) ? topz - 1 : topz;
		double sod = std::sqrt(double(output_dimension));
		for (std::size_t i=0; i<strips.size(); i++)
		{
			for (std::size_t j=0; j<output_dimension; j++)
			{
				one_d_in[j][0] = strips[i][j][0];
				one_d_in[j][1] = strips[i][j][1];
			}
			fftw_execute(p1);
			for (std::size_t j=0; j<output_dimension; j++)
			{ // also for negative DFT frequencies
				strips[i][j][0] = one_d_out[(j+topz+1)%output_dimension][0]/sod;
				strips[i][j][1] = one_d_out[(j+topz+1)%output_dimension][1]/sod;
			}
		}
		// put strips onto a plane
		debug("putting each strip onto a plane");
		sod = double(output_dimension);
		if ( tpoints.size() == 0 )
		{
			if ( angles.size() != 0 && angles.size() != strips.size() )
			{
				debug(std::string(argv[0]),"angles are inconsistent");
			}
			else if ( angles.size() == 0 )
			{
				angles.resize(strips.size());
				for (std::size_t i=0; i<angles.size(); i++)
				{
					angles[i] = i * 180 / angles.size();
				}
			}
			std::vector< std::valarray<double> > opoints(output_dimension);
			for (std::size_t i=0; i<opoints.size(); i++)
			{
				opoints[i].resize(2);
				opoints[i][0] = double(i) - double(DC);
				opoints[i][1] = 0;
			}
			tpoints.resize(opoints.size()*strips.size());
			for (std::size_t i=0; i<strips.size(); i++)
			{
				Rotation<double> prot(2);
				prot.set(angles[i]);
				for (std::size_t j=0; j<opoints.size(); j++)
				{
					tpoints[i*opoints.size()+j].resize(2);
					tpoints[i*opoints.size()+j] = prot.O(opoints[j]);
				}
			}
		}
		std::valarray<double> slice_plane_re(output_dimension*strips.size());
		std::valarray<double> slice_plane_im(output_dimension*strips.size());
		for (std::size_t i=0; i<strips.size(); i++)
		{
			for (std::size_t j=0; j<output_dimension; j++)
			{
				slice_plane_re[i*output_dimension+j] = strips[i][j][0];
				slice_plane_im[i*output_dimension+j] = strips[i][j][1];
			}
		}
		// interpolate the data on the plane (needs to be put into real data first)
		debug("interpolating the data on the plane");
		if ( tpointso.size() == 0 )
		{
			tpointso.resize(tpoints.size());
			for (std::size_t i=0; i<tpointso.size(); i++)
			{
				tpointso[i].resize(tpoints[i].size());
				tpointso[i] = tpoints[i];
			}
			nngridr(tpointso, slice_plane_re, output_dimension, output_dimension);
		}
		else
		{
			nngridrc(tpoints, slice_plane_re, output_dimension, output_dimension);
		}
		nngridrc(tpoints, slice_plane_im, output_dimension, output_dimension);
		for (std::size_t j=0; j<output_dimension; j++)
		{
			for (std::size_t i=0; i<output_dimension; i++)
			{
                std::size_t p = j*output_dimension+i;
				output_fft_re[slice][p] = two_d_in[p][0] = slice_plane_re[p];
				output_fft_im[slice][p] = two_d_in[p][1] = slice_plane_im[p];
				if ( output_fft_re[slice][p] > re_max )
				{
					re_max = output_fft_re[slice][p];
					re_i = i;
					re_j = j;
				}
			}
		}
		// do reverse 2D FFT on plane
		debug("computing the 2D reverse Fourier transform of the plane");
		fftw_execute(p2);
		// save the data to the final data array
		debug("put the data into the output array");
		for (std::size_t i=0; i<output_dimension*output_dimension; i++)
		{
			output_re[slice][i] = two_d_out[i][0]/sod;
			output_im[slice][i] = two_d_out[i][1]/sod;
		}
		message("\t...done");
	}
	for (std::size_t i=0; i<strips.size(); i++)
	{
		fftw_free(strips[i]);
	}
	strips.clear();
	fftw_destroy_plan(p2);
	fftw_free(two_d_in);
	fftw_free(two_d_out);
	fftw_destroy_plan(p1);
	fftw_free(one_d_in);
	fftw_free(one_d_out);
	
	// putting the 2D data back into order
	message("putting 2D data back into order");
	for (std::size_t slice=0; slice<number_of_slices; slice++)
	{
		fftw_complex itemp[output_dimension];
		fftw_complex otemp[output_dimension];
		for (std::size_t j=0; j<output_dimension; j++)
		{
			for (std::size_t i=0; i<output_dimension; i++)
			{
				itemp[i][0] = output_re[slice][j*output_dimension+i];
				itemp[i][1] = output_im[slice][j*output_dimension+i];
			}
			for (int i=0; i<output_dimension; i++)
			{
				otemp[i][0] = itemp[(i-re_i)%output_dimension][0];
				otemp[i][1] = itemp[(i-re_i)%output_dimension][1];
			}
			for (std::size_t i=0; i<output_dimension; i++)
			{
				output_re[slice][j*output_dimension+i] = otemp[i][0];
				output_im[slice][j*output_dimension+i] = otemp[i][1];
			}
		}
		for (std::size_t i=0; i<output_dimension; i++)
		{
			for (std::size_t j=0; j<output_dimension; j++)
			{
				itemp[j][0] = output_re[slice][j*output_dimension+i];
				itemp[j][1] = output_im[slice][j*output_dimension+i];
			}
			for (int j=0; j<output_dimension; j++)
			{
				otemp[j][0] = itemp[(j-re_j)%output_dimension][0];
				otemp[j][1] = itemp[(j-re_j)%output_dimension][1];
			}
			for (std::size_t j=0; j<output_dimension; j++)
			{
				output_re[slice][j*output_dimension+i] = otemp[j][0];
				output_im[slice][j*output_dimension+i] = otemp[j][1];
			}
		}
	}
	
	// write the output data
	message("writing output PNG files");
	b.resize(output_re.size()*output_dimension*output_dimension);
	// 2D FFT
	for (std::size_t i=0; i<output_fft_re.size(); i++)
	{
		for (std::size_t j=0; j<output_fft_re[i].size(); j++)
		{
			b[std::slice(i*output_dimension*output_dimension,output_dimension*output_dimension,1)] = output_fft_re[i];
		}
	}
	pngwrite(the_dataname+"-fft2_re.png", output_dimension, output_dimension, b, tpointso);
	for (std::size_t i=0; i<output_fft_im.size(); i++)
	{
		for (std::size_t j=0; j<output_fft_im[i].size(); j++)
		{
			b[std::slice(i*output_dimension*output_dimension,output_dimension*output_dimension,1)] = output_fft_im[i];
		}
	}
	pngwrite(the_dataname+"-fft2_im.png", output_dimension, output_dimension, b, tpointso);
	// 2D data
	for (std::size_t i=0; i<output_re.size(); i++)
	{
		for (std::size_t j=0; j<output_re[i].size(); j++)
		{
			b[std::slice(i*output_dimension*output_dimension,output_dimension*output_dimension,1)] = output_re[i];
		}
	}
	pngwrite(the_dataname+"-inverted_re.png", output_dimension, output_dimension, b, tpointso);
	for (std::size_t i=0; i<output_im.size(); i++)
	{
		for (std::size_t j=0; j<output_im[i].size(); j++)
		{
			b[std::slice(i*output_dimension*output_dimension,output_dimension*output_dimension,1)] = output_im[i];
		}
	}
	pngwrite(the_dataname+"-inverted_im.png", output_dimension, output_dimension, b, tpointso);
	for (std::size_t i=0; i<output_re.size(); i++)
	{
		for (std::size_t j=0; j<output_re[i].size(); j++)
		{
			output_re[i][j] = std::sqrt(output_re[i][j]*output_re[i][j] + output_im[i][j]*output_im[i][j]);
		}
	}
	for (std::size_t i=0; i<output_re.size(); i++)
	{
		for (std::size_t j=0; j<output_re[i].size(); j++)
		{
			b[std::slice(i*output_dimension*output_dimension,output_dimension*output_dimension,1)] = output_re[i];
		}
	}
	pngwrite(the_dataname+"-inverted.png", output_dimension, output_dimension, b, tpointso);
	message("\t...done");
	
	// finish
	message("FINISHED running " + std::string(argv[0]));
	message("please wait until the programme exits");
	return 0;
}
