/**
 A programme to perform a linear tomographic inversion on any greyscale PNG file
 
 Rado Faletic
 21st April 2005
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

int main(int argc, char* argv[])
{
	std::string ifilename = "input.png";
	std::string mfilename = "";
	std::string ip3dgfilename = "";
	std::string op3dfilename = "output";
    std::size_t int_res = 0;
	interpolation_method iterp = NN;
	bool ramp_filter = false;
	bool write_intermediates = false;
	bool fft_on = true;
	bool inv_on = true;
    std::size_t iterations = 2500;
	double damping = 0.0;
    double smoothing = 1.0;
    std::size_t nangles = 0;
	std::string stretch = "middle";
	bool all_positives = true;
	bool compression = true;
	
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
			message("--input=<inputfile>\n\tthe PNG data file ("+ifilename+")");
			message("--input_g=<inputfile>\n\tthe Plot3D grid file (none)");
			message("--matrix=<inputfile>\n\tread in a preexisting matrix ("+mfilename+")");
			message("--method=fft/inv\n\twhich inversion method to use (fft & inv)");
			message("--output=<outputfile>\n\tthe output Plot3D/PNG file, without file extension ("+op3dfilename+")");
			message("--iterations=<iterations>\n\tnumber of interations to use for inv ("+std::to_string(iterations)+")");
			message("--damping=<damping>\n\tmatrix damping value for inv ("+std::to_string(damping)+")");
			message("--smoothing=<smoothing>\n\tmatrix smoothing value for inv ("+std::to_string(smoothing)+")");
			message("--positive=yes/no\n\tshould the matrix solution be positive (yes)");
			message("--compression=yes/no\n\tcompress the matrix before solving? (yes)");
			message("--resolution=<n>\n\tresolution for FFT interpolation, determined by input ("+std::to_string(int_res)+")");
			message("--interpolation=<method>\n\twhich interpolation method to use (nn)");
			message("--filter=none/ramp\n\twhich filter to use (none)");
			message("--intermediates=on/off\n\twrite PNG files of intermediate steps (off)");
			message("--stretch=inner/middle/outer\n\tthe inversion domain to use ("+stretch+")");
			message("--nangles=<n>\n\tnumber of angles for an axi-symmetric case ("+std::to_string(nangles)+")");
			return 1;
		}
		else if ( fswitch[i].var("input") )
		{
			ifilename = fswitch[i].val();
		}
		else if ( fswitch[i].var("input_g") )
		{
			ip3dgfilename = fswitch[i].val();
		}
		else if ( fswitch[i].var("matrix") )
		{
			mfilename = fswitch[i].val();
		}
		else if ( fswitch[i].var("output") || fswitch[i].var("o") )
		{
			op3dfilename = fswitch[i].val();
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
		else if ( fswitch[i].var("positive")  )
		{
			if ( fswitch[i].val() == "yes" )
			{
				all_positives = true;
			}
			else
			{
				all_positives = false;
			}
		}
		else if ( fswitch[i].var("compression")  )
		{
			if ( fswitch[i].val() == "yes" )
			{
				compression = true;
			}
			else
			{
				compression = false;
			}
		}
		else if ( fswitch[i].var("resolution") || fswitch[i].var("res") )
		{
			std::istringstream choice(fswitch[i].val());
			choice >> int_res;
		}
		else if ( fswitch[i].var("interpolation") || fswitch[i].var("interpolate") )
		{
			if ( fswitch[i].val() == std::string("BILINEAR") || fswitch[i].val() == std::string("bilinear") )
			{
				iterp = BILINEAR;
			}
			else
			{
				iterp = NN;
			}
		}
		else if ( fswitch[i].var("filter")  )
		{
			ramp_filter = ( fswitch[i].val() == "ramp" ) ? true : false;
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
		else if ( fswitch[i].var("stretch") )
		{
			stretch = fswitch[i].val();
		}
		else if ( fswitch[i].var("nangles") )
		{
			std::istringstream choice(fswitch[i].val());
			choice >> nangles;
		}
		else
		{
			message("Unrecognised option '"+fswitch[i].var()+"'.\nUse the --help option to learn more.");
			throw;
            return 1;
		}
	}
	
	// read PNG image file
	message("reading '"+ifilename+"'");
    std::size_t xresolution;
    std::size_t yresolution;
	std::valarray<double> b;
	std::valarray<bool> blanks;
	Angle::axes rot_axis;
	std::valarray<double> angles;
    double scale;
	bool realdata;
	Tomography::pngread(ifilename, yresolution, xresolution, b, blanks, rot_axis, angles, scale, realdata);
	if ( nangles && xresolution*yresolution == b.size() ) // axi-symmetric case
	{
		angles.resize(nangles);
		std::valarray<double> single_b = b;
		std::valarray<bool> single_blanks = blanks;
		b.resize(nangles*single_b.size());
		blanks.resize(nangles*single_blanks.size());
		for (std::size_t i=0; i<nangles; i++)
		{
			angles[i] = ((double)(180*i))/((double)(nangles));
			b[std::slice(i*single_b.size(), single_b.size(), 1)] = single_b;
			blanks[std::slice(i*single_blanks.size(), single_blanks.size(), 1)] = single_blanks;
		}
	}
	else
	{
		if ( !angles.size() )
		{
			angles.resize(b.size()/(xresolution*yresolution));
			for (std::size_t i=0; i<angles.size(); i++)
			{
				angles[i] = (double(180*i)) / (double(angles.size()));
			}
		}
		nangles = angles.size();
	}
	
	// Fourier inversion
	if ( fft_on )
	{
        std::size_t resolution;
        std::size_t nslices;
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
            default:
                break;
		}
        std::size_t ires = ( int_res ) ? int_res : resolution;
		std::valarray<double>& fb = b;
		std::valarray<double> fangles = angles;
		
		// create series of sinograms
		std::vector< std::valarray<double> > sinogram(nslices, std::valarray<double>(nangles*resolution));
		for (std::size_t i=0; i<nslices; i++)
		{
			for (std::size_t j=0; j<nangles; j++)
			{
				switch(rot_axis)
				{
					case Angle::X: case::Angle::YZ:
						sinogram[i][std::slice(j*resolution,resolution,1)] =
						fb[std::slice(j*resolution*nslices+i,resolution,nslices)];
						break;
					case Angle::Y: case::Angle::ZX:
						sinogram[i][std::slice(j*resolution,resolution,1)] =
						fb[std::slice(j*resolution*nslices+i*resolution,resolution,1)];
						break;
                    default:
                        break;
				}
			}
		}
		if ( write_intermediates ) // save sinograms
		{
			std::string o = op3dfilename + "_1-sinograms.png";
			message("saving '"+o+"'");
			std::valarray<double> osinogram(sinogram.size()*nangles*resolution);
			for (std::size_t i=0; i<sinogram.size(); i++)
			{
				osinogram[std::slice(i*nangles*resolution, nangles*resolution, 1)] = sinogram[i];
			}
			Tomography::pngwrite(o, nangles, resolution, osinogram, rot_axis, fangles, scale, false, false);
		}
		
		// perform 1D Fourier Transform on each tomogram in each sinogram
		message("performing 1D FFT on tomograms");
		fftw_complex* one_d_in = (fftw_complex*)fftw_malloc(resolution * sizeof(fftw_complex));
		fftw_complex* one_d_out = (fftw_complex*)fftw_malloc(resolution * sizeof(fftw_complex));
		fftw_plan p1 = fftw_plan_dft_1d(resolution, one_d_in, one_d_out, FFTW_FORWARD, FFTW_MEASURE);
		std::vector< std::valarray<double> > fft_re(sinogram.size(), std::valarray<double>(nangles*resolution));
		std::vector< std::valarray<double> > fft_im(sinogram.size(), std::valarray<double>(nangles*resolution));
        std::size_t DC = ( !(resolution%2) ) ? (resolution-2)/2 : (resolution-1)/2;
		for (std::size_t slice=0; slice<sinogram.size(); slice++)
		{
			for (std::size_t j=0; j<nangles; j++)
			{
				for (std::size_t i=0; i<resolution; i++)
				{
					// shift the data so that negative parts are tagged on the end, rather than the beginning
					one_d_in[i][0] = sinogram[slice][j*resolution+(i+DC)%resolution];
					one_d_in[i][1] = 0;
				}
				fftw_execute(p1);
				for (std::size_t i=0; i<resolution; i++)
				{
					// apply the filter, if selected this is the standard ramp filter multiplied by a Butterworth low-pass filter
                    double filt = ( ramp_filter ) ? std::abs(double(i-DC)) / std::sqrt(1+std::pow((std::abs(double(i-DC))/double(resolution-DC)),10)) : double(1);
					// be sure to shift the DFT output so that negative frequencies come before the DC (ready for 2D interpolation)
					fft_re[slice][j*resolution+i] = filt * one_d_out[(i+(resolution-DC))%resolution][0];
					fft_im[slice][j*resolution+i] = filt * one_d_out[(i+(resolution-DC))%resolution][1];
				}
			}
		}
		fftw_destroy_plan(p1);
		fftw_free(one_d_in);
		fftw_free(one_d_out);
		
		if ( write_intermediates ) // saving 1d Fourier Transforms
		{
			std::string o1fft = op3dfilename + "_2-fft_tomograms.png";
			message("saving '"+o1fft+"'");
			std::valarray<double> fout(2*sinogram.size()*nangles*resolution);
			for (std::size_t i=0; i<sinogram.size(); i++)
			{
				fout[std::slice(2*i*nangles*resolution,nangles*resolution,1)] = fft_re[i];
				fout[std::slice(2*i*nangles*resolution+nangles*resolution,nangles*resolution,1)] = fft_im[i];
			}
			Tomography::pngwrite(o1fft, nangles, resolution, fout, rot_axis, fangles, scale, false, false);
		}
		
		// interpolate in 2D Fourier space
		switch (iterp)
		{
			case NN:
				message("natural neighbour interpolating in 2D Fourier space, "+std::to_string(sinogram.size())+" slices");
				for (std::size_t i=0; i<sinogram.size(); i++)
				{
					counter("slice ", sinogram.size(), i+1);
					nn_interpolate(fft_re[i], fangles, ires, ires);
					nn_interpolate(fft_im[i], fangles, ires, ires);
				}
				break;
			case BILINEAR:
				message("bilinear interpolating in 2D Fourier space, "+std::to_string(sinogram.size())+" slices");
				for (std::size_t i=0; i<sinogram.size(); i++)
				{
					counter("slice ", sinogram.size(), i+1);
					bilinear_interpolate(fft_re[i], fangles, ires, ires);
					bilinear_interpolate(fft_im[i], fangles, ires, ires);
				}
				break;
		}
		
		if ( write_intermediates ) // saving 2D Fourier space image
		{
			std::string o2fft = op3dfilename + "_3-fft_2d.png";
			message("saving '"+o2fft+"'");
			std::valarray<double> fout(2*sinogram.size()*ires*ires);
			for (std::size_t i=0; i<sinogram.size(); i++)
			{
				fout[std::slice(2*i*ires*ires,ires*ires,1)] = fft_re[i];
				fout[std::slice(2*i*ires*ires+ires*ires,ires*ires,1)] = fft_im[i];
			}
			Tomography::pngwrite(o2fft, ires, ires, fout, rot_axis, fangles, scale, false, false);
		}
		
		// inverse 2D Fourier Transform
		message("performing reverse 2D FFT on plane");
		fftw_complex* two_d_in = (fftw_complex*)fftw_malloc(ires * ires * sizeof(fftw_complex));
		fftw_complex* two_d_out = (fftw_complex*)fftw_malloc(ires * ires * sizeof(fftw_complex));
		fftw_plan p2 = fftw_plan_dft_2d(ires, ires, two_d_in, two_d_out, FFTW_BACKWARD, FFTW_MEASURE);
		DC = ( !(ires%2) ) ? (ires-2)/2 : (ires-1)/2;
		for (std::size_t s=0; s<sinogram.size(); s++)
		{
			counter("slice ", sinogram.size(), s+1);
			for (std::size_t i=0; i<ires; i++) // NOTE: swap x-y format from Column-major to Row-major for FFTW
			{
				for (std::size_t j=0; j<ires; j++)
				{
					// be sure to shift the x,y frequencies so that the negative frequencies are reversed
					two_d_in[i*ires+j][0] = fft_re[s][((j+DC)%ires)*ires+(i+DC)%ires];
					two_d_in[i*ires+j][1] = fft_im[s][((j+DC)%ires)*ires+(i+DC)%ires];
				}
			}
			fftw_execute(p2);
			for (std::size_t j=0; j<ires; j++) // NOTE: swap x-y format from Row-major back to Column-major
			{
				for (std::size_t i=0; i<ires; i++)
				{
					// swap the x-y components so that negative frequencies come before DC
					fft_re[s][j*ires+i] = two_d_out[((i+(ires-DC))%ires)*ires+(j+(ires-DC))%ires][0] / ( ires * ires );
					fft_im[s][j*ires+i] = two_d_out[((i+(ires-DC))%ires)*ires+(j+(ires-DC))%ires][1] / ( ires * ires );
				}
			}
		}
		fftw_destroy_plan(p2);
		fftw_free(two_d_in);
		fftw_free(two_d_out);
		if ( rot_axis == Angle::X || rot_axis == Angle::YZ )
		{
			for (std::size_t i=0; i<sinogram.size(); i++)
			{
				rotate_matrix(ires, ires, fft_re[i], false);
				rotate_matrix(ires, ires, fft_im[i], false);
				horizontal_flip(fft_re[i], ires);
				horizontal_flip(fft_im[i], ires);
			}
		}
		
		if ( write_intermediates ) // saving reverse images
		{
			std::string orf = op3dfilename + "_4-2d.png";
			message("saving '"+orf+"'");
			std::valarray<double> fout(2*sinogram.size()*ires*ires);
			for (std::size_t i=0; i<sinogram.size(); i++)
			{
				fout[std::slice(2*i*ires*ires,ires*ires,1)] = fft_re[i];
				fout[std::slice(2*i*ires*ires+ires*ires,ires*ires,1)] = fft_im[i];
			}
			Tomography::pngwrite(orf, ires, ires, fout, rot_axis, fangles, scale, false, false);
		}
		
		std::valarray<double> sqrfft(sinogram.size()*ires*ires);
		for (std::size_t i=0; i<sinogram.size(); i++)
		{
			sqrfft[std::slice(i*ires*ires,ires*ires,1)] = std::sqrt(fft_re[i]*fft_re[i] + fft_im[i]*fft_im[i]);
			//for (std::size_t j=0; j<ires*ires; j++)
			//  {
			//    sqrfft[i*ires*ires+j] = std::sqrt(fft_re[i][j]*fft_re[i][j]+fft_im[i][j]*fft_im[i][j]);
			//    if ( fft_re[i][j] < double(0) )
			//      {
			//       sqrfft[i*ires*ires+j] = double(0);
			//      }
			//  }
		}
		fft_re.resize(0);
		fft_im.resize(0);
		
		if ( write_intermediates ) // saving reverse images
		{
			message("saving '"+op3dfilename+"_5-fft.png'");
			Tomography::pngwrite(op3dfilename+"_5-fft.png", ires, ires, sqrfft, rot_axis, fangles, scale, false, false);
		}
		
        std::size_t fres = ires;
        std::size_t ssize = sinogram.size();
		DC = 0;
		if ( stretch == "inner" || stretch == "middle" )
		{
			fres = ires;
		}
		else if ( stretch == "outer" )
		{
			fres = std::size_t(double(ires)*std::sqrt(double(0.5))+double(0.5));
			ssize = fres;
            std::size_t tmp = sinogram.size() - ssize;
			DC = ( !(tmp%2) ) ? (tmp-2)/2 : (tmp-1)/2;
			DC = tmp - DC;
			DC *= ires*ires;
		}
		std::valarray<double> cfft(double(0), ssize*fres*fres);
		for (std::size_t s=0; s<ssize; s++)
		{
			for (std::size_t j=0; j<fres; j++)
			{
				cfft[std::slice(s*fres*fres+j*fres,fres,1)] =
				sqrfft[std::slice(DC+s*ires*ires+ires*((ires-fres)/2)+((ires-fres)/2)+j*ires,fres,1)];
			}
		}
		sqrfft.resize(0);
		
		// write output
		message("saving '"+op3dfilename+"_fft.png'");
		Tomography::pngwrite(op3dfilename+"_fft.png", fres, fres, cfft, rot_axis, fangles, scale, true, true);
	}
	
	// matrix inversion
	if ( inv_on )
	{
		// grid
		message("creating grid");
		grid_input mygrid;
        double scale_d = scale;
		if ( ip3dgfilename.size() )
		{
			mygrid.type() = structured;
			mygrid.format() = Binary;
			mygrid.precision() = Single;
			mygrid.multidomain() = false;
			mygrid.blanking() = false;
			mygrid.load_grid() = true;
			mygrid.gridfile() = ip3dgfilename;
		}
		else
		{
			mygrid.type() = structured;
			mygrid.load_grid() = false;
			mygrid.g_nX() = xresolution;
			mygrid.g_nY() = yresolution;
			mygrid.g_scale() = scale;
			if ( int_res )
			{
				if ( stretch == "inner" || stretch == "middle" )
				{
					switch(rot_axis)
					{
						case Angle::X: case::Angle::YZ:
							mygrid.g_nX() = std::size_t(xresolution * double(int_res) / double(yresolution));
							mygrid.g_nY() = int_res;
							scale_d *= double(yresolution) / int_res;
							break;
						case Angle::Y: case::Angle::ZX:
							mygrid.g_nX() = int_res;
							mygrid.g_nY() = std::size_t(yresolution * double(int_res) / double(xresolution));
							scale_d *= double(xresolution) / int_res;
							break;
                        default:
                            break;
					}
				}
				else if ( stretch == "outer" )
				{
					switch(rot_axis)
					{
						case Angle::X: case::Angle::YZ:
							mygrid.g_nY() = std::size_t(std::sqrt(0.5) * int_res);
							mygrid.g_nX() = std::size_t(xresolution * double(mygrid.g_nY()) / double(yresolution));
							scale_d *= std::sqrt(double(2)) * double(yresolution) / int_res;
							break;
						case Angle::Y: case::Angle::ZX:
							mygrid.g_nX() = std::size_t(std::sqrt(0.5) * int_res);
							mygrid.g_nY() = std::size_t(yresolution * double(mygrid.g_nX()) / double(xresolution));
							scale_d *= std::sqrt(double(2)) * double(xresolution) / int_res;
							break;
                        default: break;
					}
				}
			}
			switch(rot_axis)
			{
				case Angle::X: case::Angle::YZ:
					mygrid.g_nZ() = mygrid.g_nY();
					break;
				case Angle::Y: case::Angle::ZX:
					mygrid.g_nZ() = mygrid.g_nX();
					break;
                default:
                    break;
			}
		}
		grid<double> the_grid(mygrid, false);
		the_grid.give_dataname(op3dfilename);
		
		SparseMatrix<double> A(0,the_grid.ncells());
		std::valarray<double> inv_b(0);
		
		if ( mfilename.size() )
		{
			message("reading matrix '"+mfilename+"'");
			A.read(mfilename, inv_b);
			if ( inv_b.size() == b.size() ) inv_b = b;
		}
		else
		{
			message("generating grid properties");
			the_grid.auxs();
			// generate projection lines
			message("setting up "+std::to_string(nangles*xresolution*yresolution)+" rays:");
			std::vector< std::valarray<double> > ipoints(0);
			for (std::size_t j=0; j<yresolution; j++)
			{
				for (std::size_t i=0; i<xresolution; i++)
				{
					std::valarray<double> tmp = the_grid.center();
					tmp[0] += ( double(i) - double(xresolution-1) / 2 ) * scale_d;
					tmp[1] += ( double(j) - double(yresolution-1) / 2 ) * scale_d;
					ipoints.push_back(tmp);
				}
			}
			Rotation<double> rot(3);
			rot.set_origin(the_grid.center());
			std::valarray<double> islope(3);
			islope[0] = 0;
			islope[1] = 0;
			islope[2] = 1;
			std::vector< line<double> > rays(0);
			for (std::size_t i=0; i<angles.size(); i++)
			{
				rot.reset(angles[i], rot_axis);
				std::valarray<double> slope = rot.O(islope);
				for (std::size_t j=0; j<ipoints.size(); j++)
				{
					//counter("ray", rays.size()+1);
					line<double> tline(slope, rot(ipoints[j]));
					rays.push_back(tline);
				}
			}
			
			// project rays
			message("projecting "+std::to_string(rays.size())+" rays: ");
			Tomography::projection(the_grid, rays, blanks, inv_b, A, walkfast, true, inv_on, true);
			
			if ( write_intermediates ) // saving projections
			{
				std::string o = op3dfilename + "_5-projections.png";
				message("saving '"+o+"'");
				std::valarray<double> otb(double(0), rays.size());
				otb[blanks] = inv_b;
				Tomography::pngwrite(o, yresolution, xresolution, otb, rot_axis, angles, scale, false, true);
			}
			rays.clear();
			
            std::size_t bbs = 0;
			for (std::size_t i=0; i<blanks.size(); i++)
			{
				if ( blanks[i] )
				{
					bbs++;
				}
			}
			inv_b.resize(bbs);
			inv_b = b[blanks];
			
			message("saving '"+op3dfilename+".matrix'");
			A.write(op3dfilename+"_inv.matrix", inv_b);
		}
		b.resize(0);
        std::size_t Nx = the_grid.nex();
        std::size_t Ny = the_grid.ney();
        std::size_t Nz = the_grid.nez();
		Nx = ( Nx <= 1 ) ? 1 : Nx - 1;
		Ny = ( Ny <= 1 ) ? 1 : Ny - 1;
		Nz = ( Nz <= 1 ) ? 1 : Nz - 1;
		
		// adding damping and smoothing
		if ( smoothing )
		{
			message("adding smoothing equations");
			if ( mfilename.size() )
			{
				the_grid.auxs();
			}
			AddSmoothing(A, inv_b, the_grid, smoothing);
		}
		the_grid.clear();
		
		// doing matrix inversion
		std::valarray<double> x(double(0), A.cols());
        std::size_t cs = A.cols();
		bool write_with_trans = false;
		if ( compression )
		{
			message("compressing matrix");
			A.Compress(x);
		}
		message("inverting matrix");
		lsqr_input<double> linput(damping,0,0,0,iterations,all_positives);
		lsqr_output<double> Aoutput = LSQR(A, x, inv_b, linput);
		if ( compression )
		{
			message("uncompressing matrix from " + std::to_string(A.cols()) + " to " + std::to_string(cs));
			if ( A.cols() != cs )
			{
				write_with_trans = true;
			}
            double filler = ( double(0) <= x.min() ) ? double(0) : x.min()-(x.max()-x.min())/9;
			A.Uncompress(x, filler);
		}
		std::valarray<bool> Arefs = *(A.Referenced());
		A.clear();
		Aoutput.write(op3dfilename+"_inv.lsqr");
		
		if ( !compression )
		{
			// if any of the cells are unreferenced, that means that data does not exist for them,
			// so set those cells to a nominal value of the scale of the solution divided by -9
			// (so that it is 10% of the total scale)
			for (std::size_t i=0; i<Arefs.size(); i++)
			{
				if ( !Arefs[i] )
				{
					write_with_trans = true;
					break;
				}
			}
			if ( write_with_trans )
			{
				x[!Arefs] = x.min() - (x.max()-x.min()/9);
			}
			Arefs.resize(0);
		}
		
		// write output
		std::valarray<double> xo(double(0), x.size());
        std::size_t rdim = 0;
        std::size_t cdim = 0;
		switch(rot_axis)
		{
			case Angle::X: case Angle::YZ:
				for (std::size_t k=0; k<Nx; k++)
				{
					for (std::size_t j=0; j<Ny; j++)
					{
						xo[std::slice(k*Ny*Nz+j*Nz, Nz, 1)] = x[std::slice(j*Nx+k, Nz, Nx*Ny)];
					}
				}
				rdim = Ny;
				cdim = Nz;
				break;
			case Angle::Y: case Angle::ZX:
				for (std::size_t k=0; k<Ny; k++)
				{
					for (std::size_t j=0; j<Nz; j++)
					{
						xo[std::slice(k*Nz*Nx+j*Nx, Nx, 1)] = x[std::slice(k*Nx+j*Nx*Ny, Nx, 1)];
					}
				}
				rdim = Nz;
				cdim = Nx;
				break;
            default:
                break;
		}
		
		op3dfilename += "_inv";
		message("saving '"+op3dfilename+".png'");
		Tomography::pngwrite(op3dfilename+".png", rdim, cdim, xo, rot_axis, angles, scale, true, write_with_trans);
		xo.resize(0);
	}
	
	return 0;
}
