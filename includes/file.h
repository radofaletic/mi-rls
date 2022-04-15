/*
 generic defines for file attributes, especially Fortran file formats
 
 Rado Faletic
 Department of Physics
 Faculty of Science
 Australian National University  ACT  0200
 Australia
 
 Rado.Faletic@anu.edu.au
 6th August 2003
 */


#ifndef _FILE_
#define _FILE_


#include <algorithm>
#include <fstream>
#include <string>
#include <valarray>
#include <vector>


enum dataformat
{
	Formatted,
	Unformatted,
	Binary
};

enum dataprecision
{
	Single,
	Double
};

template<class T> void write_CSV(const std::string&, const std::vector< std::valarray<T> >&, const std::vector<std::string>&);


template<class T> void
write_CSV(const std::string& filename,
		  const std::vector< std::valarray<T> >& csvs,
		  const std::vector<std::string>& labls)
{
	std::string csvfilename = filename;
	if ( csvfilename.substr(csvfilename.size()-4, 4) != ".csv" && csvfilename.substr(csvfilename.size()-4, 4) != ".CSV" )
	{
		csvfilename += ".csv";
	}
	std::vector< std::valarray<T> > csv = csvs;
	std::vector<std::string> labels = labls;
	
	size_t maxlength = 0;
	for (unsigned short i=0; i<csv.size(); i++)
	{
		maxlength = std::max(maxlength, csv[i].size());
	}
	
	for (unsigned short i=0; i<csv.size(); i++)
	{
		if ( csv[i].size() != maxlength )
		{
			std::valarray<T> tmp = csv[i];
			csv[i].resize(maxlength, T(0));
			csv[i][std::slice(0,tmp.size(),1)] = tmp;
		}
	}
	
	std::ofstream csvfile(csvfilename.c_str());
	if ( labels.size() )
	{
		if ( labels.size() != csv.size() )
		{
			labels.resize(csv.size(), "data");
		}
		for (unsigned short i=0; i<labels.size(); i++)
		{
			if ( i )
			{
				csvfile << ",";
			}
			csvfile << labels[i];
		}
		csvfile << std::endl;
	}
	for (size_t j=0; j<maxlength; j++)
	{
		for (unsigned short i=0; i<csv.size(); i++)
		{
			if ( i )
			{
				csvfile << ",";
			}
			csvfile << csv[i][j];
		}
		csvfile << std::endl;
	}
	csvfile.close();
}

#endif /* _FILE_ */
