/**
 for testing Rado's png facilities
 */

#include <cmath>
#include <cstring>
#include <iostream>
#include <png.h>
#include <string>
#include <valarray>
#include <vector>

#define PNG_DEBUG 2

#include "angles.h"
#include "imglib.h"

int main(int argc, char* argv[])
{
	std::string filename = "pngtest";
	Angle::axes axis = Angle::Y;
	std::valarray<double> center;
	std::vector< std::valarray<double> > tpoints;
	std::valarray<double> angles(0);
	std::vector< std::valarray<double> > data_1(15);
    std::size_t Nrows = 193;
    std::size_t Ncols = 278;
	for (std::size_t i=0; i<data_1.size(); i++)
	{
		data_1[i].resize(Nrows*Ncols);
	}
	for (int i=0; i<Nrows; i++)
	{
		for (int j=0; j<Ncols; j++)
		{
			data_1[0][i*Ncols+j] = (i*i+j*j < 12500) ? 0 : std::sqrt(double(i*i + j*j - 12500)) + 1.0;
		}
	}
	for (std::size_t i=1; i<data_1.size(); i++)
	{
		data_1[i] = data_1[0];
	}
	std::valarray<double> data_1x(data_1.size() * Nrows * Ncols);
	for (std::size_t s=0; s<data_1.size(); s++) {
		for (int i=0; i<Nrows; i++) {
			for (int j=0; j<Ncols; j++) {
				data_1x[s * Nrows * Ncols + i * Ncols + j] = data_1[s][i * Ncols + j];
			}
		}
	}
	bool success1 = pngwrite(filename+"_1.png", Nrows, Ncols, data_1x, tpoints, axis, angles);
	if ( success1 )
	{
		std::cout << filename << "_1.png written successfully.\n";
	}
	else
	{
		std::cout << filename << "_1.png write FAILED.\n";
	}
	std::string dataname = "";
	double scale_x = 0.0;
	double scale_y = 0.0;
	std::valarray<double> data_2x;
	bool success2 = pngread(filename+"_1.png", dataname, Nrows, Ncols, scale_x, scale_y, data_2x, axis, angles);
	if ( success2 )
	{
		std::cout << "successfully read " << dataname << " from pngtest.png.\n";
	}
	else
	{
		std::cout << "pngtest.png read FAILED.\n";
	}
	bool success3 = pngwrite(filename+"_2.png", Nrows, Ncols, data_2x, tpoints, axis, angles);
	if ( success3 )
	{
		std::cout << filename << "_2.png written successfully.\n";
	}
	else
	{
		std::cout << filename << "_2.png write FAILED.\n";
	}
	
	if ( success1 && success2 && success3 )
	{
		std::vector< std::valarray<double> > data_2(data_1.size());
		for (std::size_t s=0; s<data_1.size(); s++) {
			data_2[s].resize(Nrows * Ncols);
			for (std::size_t i=0; i<Nrows * Ncols; i++) {
				data_2[s][i] = data_2x[s * Nrows * Ncols + i];
			}
		}
		
        std::size_t counter_ = 0;
		std::cout << std::endl;
		for (std::size_t i=0; i<data_1.size(); i++)
		{
			for (std::size_t j=0; j<data_1[0].size(); j++)
			{
				double pd = 100*std::abs(data_2[i][j]-data_1[i][j])/std::max(data_2[i][j],data_1[i][j]);
				if ( pd > 1.0 )
				{
					//std::cout << "data_1[" << i << "][" << j << "] = " << data_1[i][j] << "\tdata_2[" << i << "][" << j << "] = " << data_2[i][j] << "\t% = " << pd << std::endl;
					counter_++;
				}
			}
		}
		std::cout << "ERRORS: " << counter_ << " out of " << data_1.size() * data_1[0].size() << ", or " << 100.0*counter_/(data_1.size()*data_1[0].size()) << "%" << std::endl;
	}
}
