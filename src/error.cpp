/**
 A programme to calculate the error between two PNG files
 
 Rado Faletic
 26th September 2004
 22nd April 2022
 */

#include <algorithm>
#include <cmath>
#include <string>
#include <valarray>

#include "angles.h"
#include "argv.h"
#include "extra_math.h"
#include "file.h"
#include "front-end.h"
#include "tomography.h"

int main(int argc, char* argv[])
{
	std::string sourcefile = "source.png";
	std::vector< std::string > datafile(1, "datafile.png");
	std::vector< std::string > output(1, "datafile.difference.png");
	bool L2error = true;
	
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
			message("--source=<inputfile>\n\tthe original image (" + sourcefile + ")");
			message("--data=<inputfile>\n\tthe reconstructed image(s) (" + datafile[0] + ")");
			message("--type=L2\n\tthe type of error to calculate (L2)");
			return 1;
		}
		else if ( fswitch[i].var("source") )
		{
			sourcefile = fswitch[i].val();
		}
		else if ( fswitch[i].var("data") )
		{
			std::string tmp = fswitch[i].val();
            std::size_t ndf = count(tmp.begin(), tmp.end(), ',');
			if ( ndf )
			{
				datafile.clear();
				std::string::size_type oget_c = 0;
				std::string::size_type get_c = tmp.find(",");
				while ( get_c != tmp.npos )
				{
					datafile.push_back(tmp.substr(oget_c, get_c - oget_c));
					oget_c = get_c+1;
					get_c = tmp.find(",", oget_c);
				}
				datafile.push_back(tmp.substr(oget_c, tmp.size() - oget_c));
			}
			else
			{
				datafile[0] = tmp;
			}
		}
		else if ( fswitch[i].var("type") )
		{
			if ( fswitch[i].val() == "L2" )
			{
				L2error = true;
			}
		}
		else
		{
			message("Unrecognised option '" + fswitch[i].var() + "'.\nUse the --help option to learn more.");
			throw;
            return 1;
		}
	}
	if ( sourcefile.substr(sourcefile.size() - 4, 4) != ".png" && sourcefile.substr(sourcefile.size() - 4, 4) != ".PNG" )
	{
		sourcefile += ".png";
	}
	output.resize(datafile.size());
	for (std::size_t i=0; i<datafile.size(); i++)
	{
		if ( datafile[i].substr(datafile[i].size() - 4, 4) != ".png" && datafile[i].substr(datafile[i].size() - 4, 4) != ".PNG" )
		{
			datafile[i] += ".png";
		}
		output[i] = datafile[i].substr(0, datafile[i].size() - 3) + "difference.png";
	}
	
	// read original PNG image
	message("reading '" + sourcefile + "'");
    std::size_t NrowsS, NcolsS;
	std::valarray<double> source;
	Angle::axes sino_axis;
	std::valarray<double> angles;
    double scale;
	bool realdata;
	std::valarray<bool> blanks;
	Tomography::pngread(sourcefile, NrowsS, NcolsS, source, blanks, sino_axis, angles, scale, realdata);
	
	// read PNG image file(s)
    std::size_t NrowsD, NcolsD;
	std::vector< std::valarray<double> > data(datafile.size());
	Angle::axes sino_axisD;
	std::valarray<double> anglesD;
    double scaleD;
	for (std::size_t i=0; i<data.size(); i++)
	{
		message("reading '"+datafile[i]+"'");
		Tomography::pngread(datafile[i], NrowsD, NcolsD, data[i], blanks, sino_axisD, anglesD, scaleD, realdata);
		if ( source.size() != data[i].size() || NrowsD != NrowsS || NcolsD != NcolsS )
		{
			message("files '" + sourcefile + "' and '" + datafile[i] + "' have different dimensions.");
			throw;
            return 1;
		}
	}
	
	std::vector< std::valarray<double> > difference = data;
	std::for_each(difference.begin(), difference.end(), SubtractVector<double>(source));
	std::for_each(difference.begin(), difference.end(), Abs<double>);
	
	std::valarray<double> error(double(0), difference.size());
	
	if ( L2error )
	{ // use the L2 error (from Toft, pg 151)
		for (std::size_t i=0; i<error.size(); i++)
		{
			error[i] = norm(difference[i]) / norm(source);
		}
	}
	
	// write the un-normalised difference maps
	for (std::size_t i=0; i<difference.size(); i++)
	{
		message("writing '"+output[i]+"'");
		Tomography::pngwrite(output[i], NrowsS, NcolsS, difference[i], sino_axis, angles, scale, true, false);
		difference[i].resize(0);
	}
	difference.clear();
	
	std::vector<std::string> labels(datafile.size() + 1);
	labels[0] = sourcefile.substr(0, sourcefile.size() - 4);
	for (std::size_t i=0; i<datafile.size(); i++)
	{
		labels[i+1] = datafile[i].substr(0, datafile[i].size() - 4) + " " + std::to_string(error[i]);
	}
	
	std::vector< std::valarray<double> > cv(0);
	if ( NrowsS*NcolsS != source.size() ) // we have a series of slices from a 3d image
	{
        std::size_t slices = source.size() / ( NrowsS * NcolsS );
		cv.push_back(source[std::slice(0, slices, (NrowsS * NcolsS + NcolsS + 1) % source.size())]);
		for (std::size_t i=0; i<data.size(); i++)
		{
			cv.push_back(data[i][std::slice(0, slices, (NrowsD * NcolsD + NcolsD + 1) % data[i].size())]);
		}
	}
	else // a normal image
	{
		cv.push_back(source[std::slice(0, NrowsS, NcolsS + 1)]);
		for (std::size_t i=0; i<data.size(); i++)
		{
			cv.push_back(data[i][std::slice(0, NrowsD, NcolsD + 1)]);
		}
	}
	source.resize(0);
	data.clear();
	
	sourcefile = datafile[0].substr(0, datafile[0].size() - 8) + ".csv";
	message("writing '" + sourcefile + "'");
	write_CSV(sourcefile, cv, labels);
	
	for (std::size_t i=0; i<datafile.size(); i++)
	{
		message("L2 error for " + datafile[i] + ": " + std::to_string(error[i]));
	}
	
	return 0;
}
