/**
 A programme to calculate the error between two PNG files
 
 Rado Faletic
 12th November 2004
 22nd April 2022, updated to c++20
 */

#include <sstream>
#include <string>
#include <valarray>
#include <vector>

#include "argv.h"
#include "file.h"
#include "front-end.h"

int main(int argc, char* argv[])
{
	std::string ifile = "input.csv";
	std::string ofile = "";
    std::size_t factor = 2;
	
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
			message("--input=<inputfile>\n\tthe original CSV ("+ifile+")");
			message("--output=<outputfile>\n\tthe output CSV (determined by input)");
			message("--factor=<n>\n\texpansion factor ("+std::to_string(factor)+")");
			return 1;
		}
		else if ( fswitch[i].var("input") )
		{
			ifile = fswitch[i].val();
		}
		else if ( fswitch[i].var("output") )
		{
			ofile = fswitch[i].val();
		}
		else if ( fswitch[i].var("factor") )
		{
			std::istringstream choice(fswitch[i].val());
			choice >> factor;
		}
		else
		{
			message("Unrecognised option '"+fswitch[i].var()+"'.\nUse the --help option to learn more.");
			throw;
            return 1;
		}
	}
	if ( ifile.substr(ifile.size()-4,4) != ".csv" && ifile.substr(ifile.size()-4,4) != ".CSV" )
	{
		ifile += ".csv";
	}
	if ( !ofile.size() )
	{
		ofile = ifile.substr(0,ifile.size()-4)+"_expand.csv";
	}
	if ( ofile.substr(ofile.size()-4,4) != ".csv" && ofile.substr(ofile.size()-4,4) != ".CSV" )
	{
		ofile += ".csv";
	}
	
	std::vector<double> tmpv(0);
	std::ifstream input(ifile.c_str());
	std::string tmps;
    double filler = double(0);
	while (std::getline(input, tmps))
	{
        double tmp = std::stod(tmps);
		tmpv.push_back(tmp);
		if ( tmp <= filler ) filler = tmp;
	}
	std::vector< std::valarray<double> > cv(1, std::valarray<double>(filler,factor*tmpv.size()));
	for (std::size_t i=0; i<tmpv.size(); i++)
	{
		cv[0][factor*i] = tmpv[i];
		if ( i )
		{
			for (std::size_t j=1; j<=factor-1; j++)
			{
				cv[0][factor*(i-1)+j] = cv[0][factor*(i-1)] + (double(j)/double(factor))*(cv[0][factor*i]-cv[0][factor*(i-1)]);
			}
		}
	}
	std::vector<std::string> label(1,"");
	
	write_CSV(ofile, cv, label);
	
	return 0;
}
