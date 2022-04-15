#if DEBUG >= 5
#ifndef USE_MESSAGES
#define USE_MESSAGES
#endif /* USE_MESSAGES */
#endif /* DEBUG */

#include <cstring>
#include <iostream>
#include <string>
#include <valarray>
#include <vector>
#include "imglib.h"
#include "nngridr.h"

#ifdef USE_MESSAGES
#include <fstream>
std::ofstream messages;
#endif /* USE_MESSAGES */

int main(int argc, char* argv[])
{
#ifdef USE_MESSAGES
	messages.open((std::string(argv[0])+std::string(".messages")).c_str());
#endif /* USE_MESSAGES */
	
	const int Nrows = 22;
	const int Ncols = 22;
	const std::string argcv = ( argc > 1 ) ? argv[1] : "";
	
	std::ifstream filep(argcv.c_str());
	if ( !filep )
	{
		debug(std::string(argv[0]),"Unable to open data file "+argcv);
		return 1;
	}
	std::vector< std::valarray<double> > points;
	std::vector< std::valarray<double> > data;
	data.resize(1);
	while ( !filep.eof() )
	{
		std::valarray<double> tempp(2);
		filep >> tempp[0];
		if ( filep.eof() )
		{
			continue;
		}
		filep >> tempp[1];
		if ( filep.eof() )
		{
			continue;
		}
		double tempd;
		filep >> tempd;
		if ( filep.eof() )
		{
			continue;
		}
		points.push_back(tempp);
		std::valarray<double> temp = data[0];
		data[0].resize(temp.size()+1, tempd);
		data[0][std::slice(0,temp.size(),1)] = temp;
	}
	filep.close();
	
	nngridr(points, data[0], Nrows, Ncols);
	
	const std::string png_file = argcv + ".png";
	std::valarray<double> png_angles(double(0), 1);
	std::valarray<double> dataX(data.size() * data[0].size());
	for (size_t s=0; s<data.size(); s++) {
		for (size_t i=0; i<data[s].size(); i++) {
			dataX[s * data[s].size()] = data[s][i];
		}
	}
	pngwrite(png_file, Nrows, Ncols, dataX, points, Angle::X, png_angles);
}
