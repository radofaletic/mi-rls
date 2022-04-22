/**
 A programme to calculate the error between two Plot3D Q-files
 
 Rado Faletic
 11th August 2004
 22nd April 2022
 */

#include <cmath>
#include <string>
#include <valarray>

#include "argv.h"
#include "extra_math.h"
#include "file.h"
#include "front-end.h"
#include "plot3d.h"

int main(int argc, char* argv[])
{
	std::string pfile = "source.PBG";
	std::string qfile = "source.PBS";
	std::string datafile1 = "datafile.PBS";
	std::string datafile2 = "";
	std::string output1 = "datafile.error.PBS";
	std::string output2 = "";
	
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
			message("--source=<inputfile>\n\tthe original Plot3D file (" + pfile + "," + qfile + ")");
			message("--data=<inputfile>\n\tthe reconstructed Plot3D images (" + datafile1 + ")");
			return 1;
		}
		else if ( fswitch[i].var("source") )
		{
			std::string tmp = fswitch[i].val();
			std::string::size_type get_c = tmp.find(",");
			if ( get_c != tmp.npos )
			{
				pfile = tmp.substr(0, get_c);
				qfile = tmp.substr(get_c + 1, tmp.size() - pfile.size() - 1);
			}
			else
			{
				pfile = "";
				qfile = tmp;
			}
		}
		else if ( fswitch[i].var("data") )
		{
			std::string tmp = fswitch[i].val();
			std::string::size_type get_c = tmp.find(",");
			if ( get_c != tmp.npos )
			{
				datafile1 = tmp.substr(0, get_c);
				datafile2 = tmp.substr(get_c + 1, tmp.size() - datafile1.size() - 1);
			}
			else
			{
				datafile1 = tmp;
			}
		}
		else
		{
			message("Unrecognised option '" + fswitch[i].var() + "'.\nUse the --help option to learn more.");
			throw;
            return 1;
		}
	}
	output1 = ( datafile1[datafile1.size() - 4] == '.' ) ? datafile1.substr(0, datafile1.size() - 3) : datafile1 + ".";
	output1 += "error.PBS";
	if ( datafile2.size() )
	{
		output2 = ( datafile2[datafile2.size() - 4] == '.' ) ? datafile2.substr(0, datafile2.size() - 3) : datafile2 + ".";
		output2 += "error.PBS";
	}
	
	// read original PNG image
	message("reading '" + pfile + "," + qfile + "'");
	std::valarray<double> X, Y, Z;
	std::valarray<bool> B;
	std::valarray<std::size_t> Nx, Ny, Nz;
	Plot3D::read(X, Y, Z, B, Nx, Ny, Nz, pfile, Binary, Single, false, false);
	std::valarray<double> source;
	std::vector< std::valarray<double> > zero(4);
	Plot3D::read_data(source, zero[0], zero[1], zero[2], zero[3], qfile, Binary, Single, false);
	
	// read PNG image file 1
	message("reading '"+datafile1+"'");
	std::valarray<double> data1;
	Plot3D::read_data(data1, zero[0], zero[1], zero[2], zero[3], datafile1, Binary, Single, false);
	if ( source.size() != data1.size() )
	{
		message("files '" + pfile + "," + qfile + "' and '" + datafile1 + "' have different dimensions.");
		throw;
        return 1;
	}
	
	std::valarray<double> data2;
	if ( datafile2.size() )
	{
		Plot3D::read_data(data2, zero[0], zero[1], zero[2], zero[3], datafile2, Binary, Single, false);
		if ( source.size() != data2.size() )
		{
			message("files '" + pfile + "," + qfile + "' and '" + datafile2 + "' have different dimensions.");
			throw;
            return 1;
		}
	}
	
	// use the L2 error (from Toft, pg 151)
	std::valarray<double> difference1 = std::abs(data1-source);
	data1.resize(0);
    double L2error1 = norm(difference1) / norm(source);
	std::valarray<double> difference2(0);
    double L2error2 = double(0);;
	if ( datafile2.size() )
	{
		difference2.resize(data2.size());
		difference2 = std::abs(data2-source);
		data2.resize(0);
		L2error2 = norm(difference2) / norm(source);
	}
	
	for (std::size_t i=0; i<zero.size(); i++)
	{
		zero[i] = std::valarray<double>(double(0), zero[i].size());
	}
	
	message("writing '" + output1 + "'");
	Plot3D::write_data(difference1, zero[0], zero[1], zero[2], zero[3], Nx[0], Ny[0], Nz[0], output1, Binary);
	difference1.resize(0);
	if ( datafile2.size() )
	{
		message("writing '" + output2 + "'");
		Plot3D::write_data(difference2, zero[0], zero[1], zero[2], zero[3], Nx[0], Ny[0], Nz[0], output2, Binary);
		difference2.resize(0);
	}
	
	message("L2 error for " + datafile1 + ": " + std::to_string(L2error1));
	if ( datafile2.size() )
	{
		message("L2 error for " + datafile2 + ": " + std::to_string(L2error2));
	}
	
	return 0;
}
