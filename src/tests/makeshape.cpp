
#include <fftw3.h>
#include <valarray>

#include "angles.h"
#include "tomography.h"

int main(int argc, char* argv[])
{
    const std::size_t sos = 100;
	std::valarray<double> shape(0.0, sos*sos);
	for (std::size_t j=0; j<sos; j++)
	{
		for (std::size_t i=0; i<sos; i++)
		{
			if ( ( sos/4 <= i+1 && i+1 < sos/2 && sos/10 <= j+1 && j+1 <= sos*9/10 ) ||
				( sos/2 <= i+1 && i+1 <= sos*3/4 && sos/10 <= j+1 && j+2 < sos/2 ) )
			{
				shape[j*sos+i] = 1.0;
			}
		}
	}
	
	std::valarray<double> angles(0);
	double scale = 1.0;
	Angle::axes sino_axis = Angle::Y;
	Tomography::pngwrite("shape_orig.png", sos, sos, shape, sino_axis, angles, scale, true, false);
	
	fftw_complex* two_d_in = (fftw_complex*)fftw_malloc(sos * sos * sizeof(fftw_complex));
	fftw_complex* two_d_out = (fftw_complex*)fftw_malloc(sos * sos * sizeof(fftw_complex));
	fftw_plan p2 = fftw_plan_dft_2d(sos, sos, two_d_in, two_d_out, FFTW_FORWARD, FFTW_MEASURE);
    std::size_t DC = ( !(sos%2) ) ? (sos-2)/2 : (sos-1)/2;
	for (std::size_t i=0; i<sos; i++) // NOTE: swap x-y format from Column-major to Row-major for FFTW
	{
		for (std::size_t j=0; j<sos; j++)
		{
			// be sure to shift the x,y frequencies so that the negative frequencies are reversed
			two_d_in[i*sos+j][0] = shape[((j+DC)%sos)*sos+(i+DC)%sos];
			two_d_in[i*sos+j][1] = 0.0;
		}
	}
	fftw_execute(p2);
	DC = sos - DC;
	std::valarray<double> fft_re(0.0, sos*sos);
	std::valarray<double> fft_im(0.0, sos*sos);
	for (std::size_t j=0; j<sos; j++)  // NOTE: swap x-y format from Row-major back to Column-major
	{
		for (std::size_t i=0; i<sos; i++)
		{
			// swap the x-y components so that negative frequencies come before DC
			fft_re[j*sos+i] = two_d_out[((i+DC)%sos)*sos+(j+DC)%sos][0] / sos;
			fft_im[j*sos+i] = two_d_out[((i+DC)%sos)*sos+(j+DC)%sos][1] / sos;
		}
	}
	fftw_destroy_plan(p2);
	fftw_free(two_d_in);
	fftw_free(two_d_out);
	
	std::valarray<double> fout(fft_re.size()+fft_im.size());
	fout[std::slice(0,fft_re.size(),1)] = fft_re;
	fout[std::slice(fft_re.size(),fft_im.size(),1)] = fft_im;
	Tomography::pngwrite("shape_fft.png", sos, sos, fout, sino_axis, angles, scale, true, false);
}
