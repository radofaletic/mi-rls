/**
 A programme to perform an FFT tomographic inversion on any 3D Plot3D file
 
 Rado Faletic
 7th July 2004
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
#include "sparse_matrix.h"
#include "tomography.h"

typedef double real;

int main(int argc, char* argv[])
{
	std::string ip3dgfilename = "input.PBG";
	std::string ip3dqfilename = "input.PBS";
	short ip3dqi = 1;
	std::string opngfilename = "output.png";
	std::size_t resolution = 32;
    std::size_t nslices = 32;
    std::size_t nangles = 12;
	Angle::axes rot_axis = Angle::X;
	interpolation_method iterp = BILINEAR;
	bool write_intermediates = false;
	
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
			message("--input_g=<inputfile>\n\tthe Plot3D grid file ("+ip3dgfilename+")");
			message("--input_q=<inputfile>\n\tthe Plot3D Q file ("+ip3dqfilename+")");
			message("--input_i=<n>\n\tthe Plot3D Q variable number between 1 and 5 ("+std::to_string(ip3dqi)+")");
			message("--output=<outputfile>\n\tthe output PNG file ("+opngfilename+")");
			message("--resolution=<res>\n\tresolution of each slice ("+std::to_string(resolution)+")");
			message("--slices=<n>\n\tthe number of slices ("+std::to_string(nslices)+")");
			message("--angles=<nangles>\n\tthe number of angles, or projections ("+std::to_string(nangles)+")");
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
				throw;
                return 1;
			}
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
		else if ( fswitch[i].var("slices") )
		{
			std::istringstream choice(fswitch[i].val());
			choice >> nslices;
		}
		else if ( fswitch[i].var("angles") )
		{
			std::istringstream choice(fswitch[i].val());
			choice >> nangles;
			if ( !nangles )
			{
				message("\"angles\" must be non-zero.\nUse the --help option to learn more.");
				throw;
                return 1;
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
				throw;
                return 1;
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
			throw;
            return 1;
		}
	}
	if ( opngfilename.substr(opngfilename.size()-4,4) != ".png" && opngfilename.substr(opngfilename.size()-4,4) != ".PNG" )
	{
		opngfilename += ".png";
	}
	
	// read PNG image file
	std::valarray<double> angles;
	
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
	grid<double> the_grid(mygrid);
	the_grid.read_data(mygrid);
	
	std::string gdn = opngfilename.substr(0,opngfilename.size()-4);
	the_grid.give_dataname(gdn);
	
	// generate projection lines
	message("setting up rays");
    double dlength = 1;
    double dlength_s = 1;
	std::valarray<double> min2(2);
	std::valarray<double> max2(2);
	std::valarray<double> min3 = the_grid.min();
	std::valarray<double> max3 = the_grid.max();
	switch(rot_axis)
	{
		case Angle::X: case Angle::YZ:
			min2[0] = min3[1];
			min2[1] = min3[2];
			max2[0] = max3[1];
			max2[1] = max3[2];
			dlength = norm(&min2, &max2);
			dlength_s = std::abs(max3[0]-min3[0]);
			break;
		case Angle::Y: case Angle::ZX:
			min2[0] = min3[2];
			min2[1] = min3[0];
			max2[0] = max3[2];
			max2[1] = max3[0];
			dlength = norm(&min2, &max2);
			dlength_s = std::abs(max3[1]-min3[1]);
			break;
        default:
            break;
	}
    double scale = dlength / (resolution + 1);
    double scale_s = dlength_s / (nslices + 1);
	dlength /= 2;
	dlength_s /= 2;
	std::vector< std::valarray<double> > ipoints(resolution*nslices, std::valarray<double>(3));
    std::size_t counter = 0;
	switch(rot_axis)
	{
		case Angle::X: case::Angle::YZ:
			for (std::size_t j=1; j<resolution+1; j++)
			{
				for (std::size_t i=1; i<nslices+1; i++)
				{
					ipoints[counter][0] = i * scale_s - dlength_s;
					ipoints[counter][1] = j * scale - dlength;
					ipoints[counter][2] = 0;
					ipoints[counter] += the_grid.center();
					counter++;
				}
			}
			break;
		case Angle::Y: case::Angle::ZX:
			for (std::size_t j=1; j<nslices+1; j++)
			{
				for (std::size_t i=1; i<resolution+1; i++)
				{
					ipoints[counter][0] = i * scale - dlength;
					ipoints[counter][1] = j * scale_s - dlength_s;
					ipoints[counter][2] = 0;
					ipoints[counter] += the_grid.center();
					counter++;
				}
			}
			break;
        default:
            break;
	}
	
	angles.resize(nangles);
	for (std::size_t i=0; i<angles.size(); i++)
	{
		angles[i] = ((double)(180*i))/((double)(nangles));
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
			line<double> tline(slope, rot(ipoints[j]));
			rays.push_back(tline);
		}
	}
	
	// project rays
	message("projecting rays");
	SparseMatrix<double> A(0,the_grid.ncells());
	std::valarray<double> tproj;
	std::valarray<bool> blanks(true, rays.size());
	Tomography::projection(the_grid, rays, blanks, tproj, A, walkfast, true, false, true);
	std::valarray<double> projected(double(0), blanks.size());
	projected[blanks] = tproj;
	tproj.resize(0);
	
	if ( write_intermediates ) // save projections
	{
		std::string o = opngfilename.substr(0,opngfilename.size()-4);
		o += "_1-projections.png";
		message("saving '"+o+"'");
		switch(rot_axis)
		{
			case Angle::X: case::Angle::YZ:
				Tomography::pngwrite(o, resolution, nslices, projected, rot_axis, angles, scale, false, true);
				break;
			case Angle::Y: case::Angle::ZX:
				Tomography::pngwrite(o, nslices, resolution, projected, rot_axis, angles, scale, false, true);
				break;
            default:
                break;
		}
	}
	
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
					projected[std::slice(j*resolution*nslices+i,resolution,nslices)];
					break;
				case Angle::Y: case::Angle::ZX:
					sinogram[i][std::slice(j*resolution,resolution,1)] =
					projected[std::slice(j*resolution*nslices+i*resolution,resolution,1)];
					break;
                default:
                    break;
			}
		}
	}
	if ( write_intermediates ) // save sinograms
	{
		std::string o = opngfilename.substr(0,opngfilename.size()-4);
		o += "_2-sinograms.png";
		message("saving '"+o+"'");
		std::valarray<double> osinogram(sinogram.size()*nangles*resolution);
		for (std::size_t i=0; i<sinogram.size(); i++)
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
	std::vector< std::valarray<double> > fft_re(sinogram.size(), std::valarray<double>(nangles*resolution));
	std::vector< std::valarray<double> > fft_im(sinogram.size(), std::valarray<double>(nangles*resolution));
    std::size_t DC = ( !(resolution%2) ) ? (resolution-2)/2 : (resolution-1)/2;
    double sres = std::sqrt(double(resolution));
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
		std::string o1fft = opngfilename.substr(0,opngfilename.size()-4);
		o1fft += "_3-fft_tomograms.png";
		message("saving '"+o1fft+"'");
		std::valarray<double> fout(2*sinogram.size()*nangles*resolution);
		for (std::size_t i=0; i<sinogram.size(); i++)
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
			for (std::size_t i=0; i<sinogram.size(); i++)
			{
				nn_interpolate(fft_re[i], angles);
				nn_interpolate(fft_im[i], angles);
				switch(rot_axis)
				{
					case Angle::X: case::Angle::YZ:
						rotate_matrix(resolution, resolution, fft_re[i]);
						rotate_matrix(resolution, resolution, fft_im[i]);
						break;
                    default:
                        break;
				}
			}
			break;
		case BILINEAR:
			message("bilinear interpolating in 2D Fourier space");
			for (std::size_t i=0; i<sinogram.size(); i++)
			{
				bilinear_interpolate(fft_re[i], angles);
				bilinear_interpolate(fft_im[i], angles);
				switch(rot_axis)
				{
					case Angle::X: case::Angle::YZ:
						rotate_matrix(resolution, resolution, fft_re[i]);
						rotate_matrix(resolution, resolution, fft_im[i]);
						break;
                    default:
                        break;
				}
			}
			break;
        default:
            break;
	}
	
	if ( write_intermediates ) // saving 2D Fourier space image
	{
		std::string o2fft = opngfilename.substr(0,opngfilename.size()-4);
		o2fft += "_4-fft_2d.png";
		message("saving '"+o2fft+"'");
		std::valarray<double> fout(2*sinogram.size()*resolution*resolution);
		for (std::size_t i=0; i<sinogram.size(); i++)
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
	for (std::size_t s=0; s<sinogram.size(); s++)
	{
		for (std::size_t i=0; i<resolution; i++) // NOTE: swap x-y format from Column-major to Row-major for FFTW
		{
			for (std::size_t j=0; j<resolution; j++)
			{
				// be sure to shift the x,y frequencies so that the negative frequencies are reversed
				two_d_in[j*resolution+i][0] = fft_re[s][((i+DC)%resolution)*resolution+(j+DC)%resolution];
				two_d_in[j*resolution+i][1] = fft_im[s][((i+DC)%resolution)*resolution+(j+DC)%resolution];
			}
		}
		fftw_execute(p2);
		for (std::size_t i=0; i<resolution; i++) // NOTE: swap x-y format from Row-major back to Column-major
		{
			for (std::size_t j=0; j<resolution; j++)
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
		std::string orf = opngfilename.substr(0,opngfilename.size()-4);
		orf += "_5-2d.png";
		message("saving '"+orf+"'");
		std::valarray<double> fout(2*sinogram.size()*resolution*resolution);
		for (std::size_t i=0; i<sinogram.size(); i++)
		{
			fout[std::slice(2*i*resolution*resolution,resolution*resolution,1)] = fft_re[i];
			fout[std::slice(2*i*resolution*resolution+resolution*resolution,resolution*resolution,1)] = fft_im[i];
		}
		Tomography::pngwrite(orf, resolution, resolution, fout, rot_axis, angles, scale, false, true);
	}
	
	// write output
	message("saving '"+opngfilename+"'");
	std::valarray<double> sqrfft(sinogram.size()*resolution*resolution);
	for (std::size_t i=0; i<sinogram.size(); i++)
	{
		sqrfft[std::slice(i*resolution*resolution,resolution*resolution,1)] =
		std::sqrt(fft_re[i]*fft_re[i] + fft_im[i]*fft_im[i]);
	}
	Tomography::pngwrite(opngfilename, resolution, resolution, sqrfft, rot_axis, angles, scale, true, true);
	return 0;
}
