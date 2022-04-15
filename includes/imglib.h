/*
 imglib
 
 for creating equi-spaced rectangular planar grids in 3D
 and read/write functions for text and PNG files
 
 Rado Faletic
 Department of Physics, Faculty of Science
 Australian National University  ACT  0200
 Australia
 
 Rado.Faletic@anu.edu.au
 23rd May 2004
 */


#ifndef _IMGLIB_
#define _IMGLIB_


/* ---------- header files ---------- */
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <ctime>
#include <fstream>
#include <iostream>
#include <png.h>
#include <string>
#include <valarray>
#include <vector>
/* ---------- user header files ---------- */
#include "angles.h"
#include "conversions.h"
#include "extra_math.h"
#include "file.h"
#include "fortran_io.h"
#include "front-end.h"
/* ---------------------------------- */


/* ------------------------------------------- */
/* ---------- function declarations ---------- */
/* ------------------------------------------- */


// The function "blank_2dgrid" generates a structured
// 2D mesh in 3D. Note that Nrows is the number
// of points in the y-direction and Ncols is the
// number of points in the x-direction
template<class T> std::vector< std::valarray<T> >
blank_2dgrid(const unsigned short& = 2,                // dimension
			 const size_t& = 128, const size_t& = 128, // Nrows, Ncols
			 const T& = T(1),                          // spacing
			 const Angle::axes& = Angle::XY);          // axes

template<class T> std::vector< std::valarray<T> >
blank_2dgrid(const std::valarray<T>&,                  // shift vector
			 const size_t& =128, const size_t& =128,   // Nrows, Ncols
			 const T& = T(1),                          // spacing
			 const Angle::axes& = Angle::XY);          // axes

template<class T> std::vector< std::valarray<T> >
blank_2dgrid(size_t& Nrows, size_t& Ncols,             // Nrows, Ncols
			 const std::valarray<T>&,                  // lower left point
			 const std::valarray<T>&,                  // upper right point
			 const size_t& = 128,                      // Nrowscols
			 const Angle::axes& = Angle::XY);          // axes

// each mesh point has data associated with it,
// and here we save that data for multiple
// meshes based on the one initial blank_2dgrid mesh
//  b3d file format:
//                  number_of_projections
//                  axes
//                  Nrows Ncols
//                  rotation_origin
//                  array_of_initial_3D_points

// write initial mesh data to the file
template<class T> bool
write_b3d_begin(const std::string&,                      // filename
				const size_t&,                           // number of projections
				const Angle::axes&,                      // axes
				const std::valarray<T>&,                 // origin of rotation
				const std::vector< std::valarray<T> >&,  // points
				const size_t& Nrows, const size_t& Ncols,// Nrows, Ncols
				const dataformat& = Unformatted);

// append projection data to the file
template<class T> bool
write_b3d_append(const std::string&,                     // filename
				 const T&,                               // angle
				 const std::valarray<T>&,                // data
				 const dataformat& = Unformatted);

// read projection data from the file
template<class T> bool
read_b3d(const std::string&,                             // filename
		 Angle::axes&,                                   // axes
		 std::valarray<T>&,                              // origin
		 std::vector< std::valarray<T> >&,               // points
		 size_t&, size_t&,                               // Nrows, Ncols
		 std::valarray<T>&,                              // angles
		 std::vector< std::valarray<T> >&,               // data
		 const dataformat& = Unformatted);

// write projection data to a PNG file
template<class T> bool
pngwrite(const std::string&,                             // filename
		 const size_t&, const size_t&,                   // Nrows, Ncols
		 const std::valarray<T>&,                        // data
		 const std::vector< std::valarray<T> >&,         // grid
		 const Angle::axes& = Angle::Y,                  // rotation axis
		 const std::valarray<T>& = std::valarray<T>(0)); // angles

// read projection data from a PNG file
template<class T> bool
pngread(const std::string&,                              // filename
		std::string&,                                    // dataname
		size_t&, size_t&,                                // Nrows, Ncols
		T&, T&,                                          // sCALx, sCALy
		std::valarray<T>&,                               // data
		Angle::axes&,                                    // rotation axis
		std::valarray<T>&);                              // angles


/* ------------------------------------------ */
/* ---------- function definitions ---------- */
/* ------------------------------------------ */


/* ---------- blank_2dgrid ---------- */
/* create a regular lattice with
 Nrows nodes in the y-direction
 and Ncols nodes in the
 x-direction                   */
template<class T> std::vector< std::valarray<T> >
blank_2dgrid(const unsigned short& dim,
			 const size_t& Nrows,
			 const size_t& Ncols,
			 const T& spacing,
			 const Angle::axes& axis)
{
	std::vector< std::valarray<T> > image(Nrows*Ncols);
	for (size_t i=0; i<image.size(); i++)
	{
		image[i].resize(dim, T(0));
	}
	
	for (size_t i=0; i<Nrows; i++)
	{
		T Y = T(i) * spacing;
		for (size_t j=0; j<Ncols; j++)
		{
			if ( dim > 2 )
			{
				switch(axis)
				{
					case Angle::XY:
						image[i*Ncols+j][0] = T(j) * spacing;
						image[i*Ncols+j][1] = Y;
						break;
					case Angle::YZ:
						image[i*Ncols+j][1] = T(j) * spacing;
						image[i*Ncols+j][2] = Y;
						break;
					case Angle::ZX:
						image[i*Ncols+j][2] = T(j) * spacing;
						image[i*Ncols+j][0] = Y;
						break;
				}
			}
			else
			{
				image[i*Ncols+j][0] = T(j) * spacing;
				image[i*Ncols+j][1] = Y;
			}
		}
	}
	// grid values are held in Cartesian format,
	// but our 3D-image must be in row-col format
	for (size_t i=0; i<size_t(std::floor(T(Nrows)/2)); i++)
	{
		std::swap(image[i*Nrows], image[image.size()-1-i]);
	}
	return image;
}
/* ----------------------------- */

/* ---------- blank_2dgrid ---------- */
/* create a regular lattice with
 Nrows nodes in the X direction
 and Ncols nodes in the Y
 direction and then shift the
 whole thing in the shifter
 direction                     */
template<class T> std::vector< std::valarray<T> >
blank_2dgrid(const std::valarray<T>& shifter,
			 const size_t& Nrows,
			 const size_t& Ncols,
			 const T& spacing,
			 const Angle::axes& axis)
{
	std::vector< std::valarray<T> > image = blank_2dgrid(shifter.size(), Nrows, Ncols, spacing, axis);
	std::for_each(image.begin(), image.end(), AddVector<T>(shifter));
	return image;
}
/* ----------------------------- */

/* ---------- blank_2dgrid ---------- */
/* create a regular lattice with
 lower left corner at vmin,
 upper right corner at vmax, and
 with nn nodes along the longer
 of the two directions         */
template<class T> std::vector< std::valarray<T> >
blank_2dgrid(size_t& Nrows, size_t& Ncols,
			 const std::valarray<T>& vmin,
			 const std::valarray<T>& vmax,
			 const size_t& nn,
			 const Angle::axes& axis)
{
	std::valarray<T> themin(T(0),vmin.size());
	
	// default to the xy_image configuration
	unsigned short d1 = 0;
	unsigned short d2 = 1;
	if ( vmin.size() > 2 )
	{
		switch(axis)
		{
			case Angle::YZ:
				d1 = 1;
				d2 = 2;
				break;
			case Angle::ZX:
				d1 = 2;
				d2 = 0;
				break;
		}
	}
	
	if ( vmin.size() == 2 )
	{
		themin = vmin;
	}
	else
	{
		themin[d1] = vmin[d1];
		themin[d2] = vmin[d2];
	}
	
	T xdist = std::abs(vmax[d1] - vmin[d1]);
	T ydist = std::abs(vmax[d2] - vmin[d2]);
	
	T spacing = T(0);
	if ( xdist > ydist )
	{
		Nrows = size_t( T(nn) * ydist / xdist );
		Ncols = nn;
		spacing = xdist / T(Ncols);
	}
	else
	{
		Nrows = nn;
		Ncols = size_t( T(nn) * xdist / ydist );
		spacing = ydist / T(Nrows);
	}
	
	return blank_2dgrid(themin, Nrows, Ncols, spacing, axis);
}
/* ----------------------------- */

/* ---------- write_b3d_begin ---------- */
/* write the header for a series of 2D
 meshes in 3D                          */
template<class T> bool
write_b3d_begin(const std::string& filename,
				const size_t& no_proj,
				const Angle::axes& axes,
				const std::valarray<T>& rorigin,
				const std::vector< std::valarray<T> >& tpoints,
				const size_t& Nrows, const size_t& Ncols,
				const dataformat& format)
{
	std::ofstream file;
	switch(format)
	{
		case Unformatted: case Binary:
			file.open(filename.c_str(), std::ios_base::binary);
			break;
		default:
			file.open(filename.c_str());
			break;
	}
	bool byte_swapping = false;
	dataprecision precision = Double;
	if ( sizeof(T) <= 4 )
	{
		precision = Single;
	}
	
	// write the number of projections
	Fortran::iwrite_line(file, no_proj, byte_swapping, format);
	// write the axes of rotation
	size_t ax = 0;
	switch(axes)
	{
		case Angle::YZ:
			ax = 1;
			break;
		case Angle::ZX:
			ax = 2;
			break;
	}
	Fortran::iwrite_line(file, ax, byte_swapping, format);
	// write Nrows, Ncols
	std::valarray<size_t> dim(2);
	dim[0] = Nrows;
	dim[1] = Ncols;
	Fortran::iwrite_vector(file, dim, byte_swapping, format);
	// write the rotation origin
	Fortran::fwrite_vector(file, rorigin, byte_swapping, format, precision);
	// write the initial mesh points
	for (size_t i=0; i<tpoints.size(); i++)
	{
		Fortran::fwrite_vector(file, tpoints[i], byte_swapping, format, precision);
	}
	
	file.close();
	return true;
}
/* ------------------------------------- */

/* ---------- write_b3d_append ---------- */
/* append data for a series of 2D meshes in
 3D                                     */
template<class T> bool
write_b3d_append(const std::string& filename,
				 const T& angle,
				 const std::valarray<T>& data,
				 const dataformat& format)
{
	std::ofstream file;
	switch(format)
	{
		case Unformatted: case Binary:
			file.open(filename.c_str(), std::ios_base::app|std::ios_base::binary);
			break;
		default:
			file.open(filename.c_str(), std::ios_base::app);
			break;
	}
	bool byte_swapping = false;
	dataprecision precision = Double;
	if ( sizeof(T) <= 4 )
	{
		precision = Single;
	}
	
	// write the number of projections
	Fortran::fwrite_line(file, angle, byte_swapping, format, precision);
	// write the rotation origin
	Fortran::fwrite_vector(file, data, byte_swapping, format, precision);
	
	file.close();
	return true;
}
/* -------------------------------------- */

/* ---------- read_b3d ---------- */
/* read data for a series of 2D
 image meshes in 3D             */
template<class T> bool
read_b3d(const std::string& filename,
		 Angle::axes& axes,
		 std::valarray<T>& center,
		 std::vector< std::valarray<T> >& tpoints,
		 size_t& Nrows, size_t& Ncols,
		 std::valarray<T>& angles,
		 std::vector< std::valarray<T> >& data,
		 const dataformat& format)
{
	std::ifstream file;
	switch(format)
	{
		case Unformatted: case Binary:
			file.open(filename.c_str(), std::ios_base::binary);
			break;
		default:
			file.open(filename.c_str());
			break;
	}
	
	if ( !file )
	{
		debug("read_b3d","no file to read from");
		throw; return false;
	}
	bool byte_swapping = false;
	dataprecision precision = Double;
	if ( sizeof(T) <= 4 )
	{
		precision = Single;
	}
	
	// read the number of projections
	size_t Nrots;
	Fortran::iread_line(file, Nrots, byte_swapping, format);
	// read the axes
	size_t ax = 0;
	Fortran::iread_line(file, ax, byte_swapping, format);
	switch(ax)
	{
		case 1:
			axes = Angle::YZ;
			break;
		case 2:
			axes = Angle::ZX;
			break;
		default:
			axes = Angle::XY;
			break;
	}
	
	// read Nrows, Ncols
	std::valarray<size_t> dim;
	Fortran::iread_vector(file, 2, dim, byte_swapping, format);
	// read the origin
	Fortran::fread_vector(file, 3, center, byte_swapping, format, precision);
	// read the initial mesh points
	Nrows = dim[0];
	Ncols = dim[1];
	tpoints.resize(Nrows*Ncols);
	for (size_t i=0; i<tpoints.size(); i++)
	{
		Fortran::fread_vector(file, 3, tpoints[i], byte_swapping, format, precision);
	}
	
	// read the data
	angles.resize(Nrots);
	data.resize(Nrots);
	for (size_t i=0; i<Nrots; i++)
	{
		Fortran::fread_line(file, angles[i], byte_swapping, format, precision);
		Fortran::fread_vector(file, Nrows*Ncols, data[i], byte_swapping, format, precision);
	}
	
	file.close();
	return true;
}
/* ------------------------------ */

/* ---------- pngwrite ---------- */
template<class T> bool
pngwrite(const std::string& myfilename,
		 const size_t& iNrows, const size_t& iNcols,
		 const std::valarray<T>& input_data,
		 const std::vector< std::valarray<T> >& tpoints,
		 const Angle::axes& raxis,
		 const std::valarray<T>& my_angles)
{
	if ( tpoints.size() && ( tpoints.size() != iNrows*iNcols ) )
	{
		debug("pngwrite","inconsistent PNG image dimensions");
		throw; return false;
	}
	std::vector< std::valarray<T> > data;
	for (size_t i=0; i<input_data.size()/(iNrows*iNcols); i++)
	{
		std::valarray<T> tdata = input_data[std::slice(i*iNrows*iNcols,iNrows*iNcols,1)];
		data.push_back(tdata);
	}
	
	// the configuration for the collage of images
	png_uint_32 c_cols = data.size();
	png_uint_32 c_rows = 1;
	size_t wh_sqr = (c_cols*iNcols+c_cols)*(c_cols*iNcols+c_cols) +
	(c_rows*iNrows+c_rows)*(c_rows*iNrows+c_rows);
	
	for (size_t i=1; i<data.size()+1; i++)
	{
		for (size_t j=1; j<data.size()+1; j++)
		{
			if ( i * j < data.size() )
			{
				continue;
			}
			size_t wh_width = i*iNcols+i;
			size_t wh_height = j*iNrows+j;
			if ( wh_width < wh_height )
			{
				continue;
			}
			size_t wh_sqr_n = wh_width*wh_width + wh_height*wh_height;
			if ( wh_sqr_n < wh_sqr )
			{
				c_cols = i;
				c_rows = j;
				wh_sqr = wh_sqr_n;
			}
			else if ( wh_sqr_n == wh_sqr )
			{
				if ( i < c_cols )
				{
					c_cols = i;
				}
				if ( j < c_rows )
				{
					c_rows = j;
				}
			}
		}
	}
	
	// the size of the image collage
	png_uint_32 Nrows = ( c_rows * iNrows ) + ( c_rows - 1 );
	png_uint_32 Ncols = ( c_cols * iNcols ) + ( c_cols - 1 );
	
	// set up angle facilities
	std::valarray<T> angles = my_angles;
	bool use_angles = true;
	if ( !angles.size() || angles.size() != data.size() )
	{
		use_angles = false;
		angles.resize(data.size());
		for (size_t i=0; i<angles.size(); i++)
		{
			angles[i] = i * 180 / angles.size();
		}
	}
	
	// set up the PNG file and structures
	std::string filename = myfilename;
	std::string filetitle = myfilename;
	if ( filename.substr(filename.size()-4,4) != std::string(".png") )
	{
		filename += ".png";
	}
	else
	{
		filetitle = filename.substr(0,filename.size()-4);
	}
	png_FILE_p fp = std::fopen(filename.c_str(), "wb");
	if ( !fp )
	{
		std::fclose(fp);
		debug("pngwrite","could not open PNG file '"+filename+"'");
		throw; return false;
	}
	png_structp png_ptr =
	png_create_write_struct(PNG_LIBPNG_VER_STRING, png_voidp_NULL, png_error_ptr_NULL, png_error_ptr_NULL);
	if ( !png_ptr )
	{
		std::fclose(fp);
		debug("pngwrite","error creating PNG pointer");
		throw; return false;
	}
	png_infop info_ptr = png_create_info_struct(png_ptr);
	if ( !info_ptr )
	{
		std::fclose(fp);
		png_destroy_write_struct(&png_ptr, png_infopp_NULL);
		debug("pngwrite","error creating PNG info_pointer");
		throw; return false;
	}
	if ( setjmp(png_jmpbuf(png_ptr)) )
	{
		std::fclose(fp);
		png_destroy_write_struct(&png_ptr, &info_ptr);
		debug("pngwrite","error creating PNG buffer");
		throw; return false;
	}
	png_init_io(png_ptr, fp);
	
	// IHDR
	int bit_depth = 16; // 16bit grayscale
	png_set_IHDR(png_ptr, info_ptr, Ncols, Nrows, bit_depth, PNG_COLOR_TYPE_GRAY, PNG_INTERLACE_NONE,
				 PNG_COMPRESSION_TYPE_DEFAULT, PNG_FILTER_TYPE_DEFAULT);
	
	// tRNS
	T dmax = data[0].max();
	T dmin = data[0].min();
	for (size_t i=0; i<data.size(); i++)
	{
		dmax = std::max(dmax, data[i].max());
		dmin = std::min(dmin, data[i].min());
	}
	if ( is_zero(dmin) && !eq(dmin, dmax) )
	{
		dmin = dmax;
		for (size_t i=0; i<data.size(); i++)
		{
			for (png_uint_32 j=0; j<data[i].size(); j++)
			{
				dmin = ( data[i][j] < dmin && !is_zero(data[i][j]) ) ? data[i][j] : dmin;
			}
		}
	}
	png_bytep trans = NULL;
	png_uint_16 num_trans = 1;
	png_color_16 trans_values[1];
	trans_values[0].gray = 0;
	bool has_tRNS = false;
	if ( dmin > T(0) || dmax < T(0) )
	{
		has_tRNS = true;
	}
	png_set_tRNS(png_ptr, info_ptr, trans, num_trans, trans_values);
	
	// pCAL
	png_uint_32 no_colors = png_uint_32(ipow(2,bit_depth));
	png_int_32 X0 = 1;
	png_int_32 X1 = no_colors - 1;
	int nparams = 2; // linear_data = p0 + p1 * sample / ( X1 - X0 );
	png_charp units = "";
	T p0 = dmin - ( dmax - dmin ) * X0 / ( X1 - X0 );
	T p1 = dmax - dmin;
	if ( eq(dmax, dmin) )
	{
		p0 = T(0);
		p1 = dmax * ( X1 - X0 ) / X1;
	}
	std::string p0s = ntoIEEE(p0);
	std::string p1s = ntoIEEE(p1);
	png_charp params[2] = { (char*)p0s.c_str(), (char*)p1s.c_str() };
	png_set_pCAL(png_ptr, info_ptr, "projected", X0, X1, PNG_EQUATION_LINEAR, nparams, units, params);
	
	// sCAL
	double sCAL_width = 0.0;
	double sCAL_height = 0.0;
	if ( tpoints.size() )
	{
		size_t sCAL_wcounter = 0;
		size_t sCAL_hcounter = 0;
		for (size_t i=0; i<tpoints.size(); i++)
		{
			if ( i%Ncols )
			{
				std::valarray<T> sCALdiff = tpoints[i] - tpoints[i-1];
				sCAL_width += norm(sCALdiff);
				sCAL_wcounter++;
			}
			if ( Ncols <= i )
			{
				std::valarray<T> sCALdiff = tpoints[i] - tpoints[i-Ncols];
				sCAL_height += norm(sCALdiff);
				sCAL_hcounter++;
			}
		}
		sCAL_width /= sCAL_wcounter;
		sCAL_height /= sCAL_hcounter;
		png_set_sCAL(png_ptr, info_ptr, PNG_SCALE_METER, sCAL_width, sCAL_height);
	}
	
	// tIME
	png_time mod_time;
	png_convert_from_time_t(&mod_time, std::time(0));
	png_set_tIME(png_ptr, info_ptr, &mod_time);
	
	// tEXT
	int num_text = 10+data.size();
	png_text text_ptr[num_text];
	text_ptr[0].key = "Title";
	text_ptr[0].text = (char*)filetitle.c_str();
	text_ptr[0].text_length = filetitle.size();
	text_ptr[0].compression = PNG_TEXT_COMPRESSION_NONE;
	text_ptr[1].key = "Author";
	text_ptr[1].text = "Rado Faletic";
	text_ptr[1].text_length = 12;
	text_ptr[1].compression = PNG_TEXT_COMPRESSION_NONE;
	text_ptr[2].key = "Description";
	text_ptr[2].text = "tomographic projection";
	text_ptr[2].text_length = 22;
	text_ptr[2].compression = PNG_TEXT_COMPRESSION_NONE;
	std::string png_cyear = "(c) 2002, ";
	png_cyear += ntos(mod_time.year);
	png_cyear += " Rado Faletic";
	text_ptr[3].key = "Copyright";
	text_ptr[3].text = (char*)png_cyear.c_str();
	text_ptr[3].text_length = png_cyear.size();
	text_ptr[3].compression = PNG_TEXT_COMPRESSION_NONE;
	char* ctime = png_convert_to_rfc1123(png_ptr, &mod_time);
	text_ptr[4].key = "Creation Time";
	text_ptr[4].text = ctime;
	text_ptr[4].text_length = std::strlen(ctime);
	text_ptr[4].compression = PNG_TEXT_COMPRESSION_NONE;
	text_ptr[5].key = "Software";
	text_ptr[5].text = "imglib::pngwrite";
	text_ptr[5].text_length = 16;
	text_ptr[5].compression = PNG_TEXT_COMPRESSION_NONE;
	text_ptr[6].key = "Source";
	text_ptr[6].text = "data";
	text_ptr[6].text_length = 4;
	text_ptr[6].compression = PNG_TEXT_COMPRESSION_NONE;
	std::string png_collage_txt = ntos(data.size());
	png_collage_txt += ":";
	png_collage_txt += ntos(c_cols);
	png_collage_txt += ",";
	png_collage_txt += ntos(c_rows);
	text_ptr[7].key = "collage";
	text_ptr[7].text = (char*)png_collage_txt.c_str();
	text_ptr[7].text_length = png_collage_txt.size();
	text_ptr[7].compression = PNG_TEXT_COMPRESSION_NONE;
	text_ptr[8].key = "raxis";
	text_ptr[8].text = ( raxis == Angle::X ) ? (char*)"X" : (char*)"Y";
	text_ptr[8].text_length = 1;
	text_ptr[8].compression = PNG_TEXT_COMPRESSION_NONE;
	std::string png_angles_txt = "default";
	if ( use_angles )
	{
		png_angles_txt = vtoIEEE(angles);
	}
	text_ptr[9].key = "angles";
	text_ptr[9].text = (char*)png_angles_txt.c_str();
	text_ptr[9].text_length = png_angles_txt.size();
	text_ptr[9].compression = PNG_TEXT_COMPRESSION_NONE;
	std::vector<std::string> ds(data.size());
	std::vector<std::string> dd(data.size());
	for (size_t i=0; i<data.size(); i++)
	{
		ds[i] = "data_"+ntos(i+1);
		dd[i] = vtoIEEE(data[i]);
		text_ptr[10+i].key = (char*)ds[i].c_str();
		text_ptr[10+i].text = (char*)dd[i].c_str();
		text_ptr[10+i].text_length = dd[i].size();
		//text_ptr[10+i].compression = PNG_TEXT_COMPRESSION_zTXt;
		text_ptr[10+i].compression = PNG_TEXT_COMPRESSION_NONE;
	}
	png_set_text(png_ptr, info_ptr, text_ptr, num_text);
	
	// convert from T to integer
	std::valarray<png_uint_16> idata(Nrows*Ncols);
	for (size_t i=0; i<Nrows*Ncols; i++)
	{
		idata[i] = 0;
	}
	T m = p1 / ( X1 - X0 );
	for (png_uint_32 i=0; i<c_rows; i++)
	{
		for (png_uint_32 j=0; j<c_cols; j++)
		{
			size_t bd = i*c_cols+j;
			if ( bd < data.size() )
			{
				png_uint_32 bcol = j * iNcols + j;
				png_uint_32 brow = i * iNrows + i;
				for (png_uint_32 ii=0; ii<iNrows; ii++)
				{
					for (png_uint_32 jj=0; jj<iNcols; jj++)
					{
						T tmp = (data[bd][ii*iNcols+jj] - p0)/m;
						if ( tmp < T(0) )
						{
							tmp = T(0);
						}
						else if ( T(0) < tmp && tmp < T(X0) )
						{
							tmp = T(X0);
						}
						else if ( T(X1) < tmp )
						{
							tmp = T(X1);
						}
						if ( has_tRNS && is_zero(data[bd][ii*iNcols+jj]) )
						{
							tmp = T(0);
						}
						idata[(brow+ii) * Ncols + (bcol+jj)] = png_uint_16(std::floor(tmp));
					}
				}
			}
		}
	}
	// create row pointers
	std::valarray<png_bytep> row_pointers(Nrows);
	for (png_uint_32 i=0; i<Nrows; i++)
	{
		row_pointers[i] = (png_bytep)(&(idata[i*Ncols]));
	}
	// png_set_rows
	png_set_rows(png_ptr, info_ptr, &(row_pointers[0]));
	
	// IEND
	if ( is_nbo() )
	{
		png_write_png(png_ptr, info_ptr, PNG_TRANSFORM_IDENTITY, png_voidp_NULL);
	}
	else
	{
		png_write_png(png_ptr, info_ptr, PNG_TRANSFORM_SWAP_ENDIAN, png_voidp_NULL);
	}
	
	png_destroy_write_struct(&png_ptr, &info_ptr);
	std::fclose(fp);
	return true;
}
/* ------------------------------ */

/* ---------- pngread ---------- */
template<class T> bool
pngread(const std::string& myfilename,
		std::string& dataname,
		size_t& iNrows, size_t& iNcols,
		T& sCAL_width, T& sCAL_height,
		std::valarray<T>& output_data,
		Angle::axes& raxis,
		std::valarray<T>& angles)
{
	// set up the PNG file and structures
	std::string filename = myfilename;
	if ( filename.substr(filename.size()-4,4) != std::string(".png") &&
		filename.substr(filename.size()-4,4) != std::string(".PNG") )
	{
		filename += ".png";
	}
	png_FILE_p fp = std::fopen(filename.c_str(), "rb");
	if ( !fp )
	{
		std::fclose(fp);
		debug("pngwrite","error opening PNG file '"+filename+"'");
		throw; return false;
	}
	png_byte sig2check[8];
	std::fread(sig2check, 1, 8, fp);
	if ( png_sig_cmp(sig2check, 0, 8) )
	{
		std::fclose(fp);
		debug("pngwrite","error in signature of PNG file '"+filename+"'");
		throw; return false;
	}
	
	png_structp png_ptr =
	png_create_read_struct(PNG_LIBPNG_VER_STRING, png_voidp_NULL, png_error_ptr_NULL, png_error_ptr_NULL);
	if ( !png_ptr )
	{
		std::fclose(fp);
		debug("pngwrite","error creating PNG pointer");
		throw; return false;
	}
	
	png_infop info_ptr = png_create_info_struct(png_ptr);
	if ( !info_ptr )
	{
		std::fclose(fp);
		png_destroy_read_struct(&png_ptr, png_infopp_NULL, png_infopp_NULL);
		debug("pngwrite","error creating PNG info_pointer");
		throw; return false;
	}
	
	png_infop end_info = png_create_info_struct(png_ptr);
	if ( !end_info )
	{
		std::fclose(fp);
		png_destroy_read_struct(&png_ptr, &info_ptr, png_infopp_NULL);
		debug("pngwrite","error creating PNG end_info");
		throw; return false;
	}
	
	if ( setjmp(png_jmpbuf(png_ptr)) )
	{
		std::fclose(fp);
		png_destroy_read_struct(&png_ptr, &info_ptr, &end_info);
		debug("pngwrite","error creating PNG buffer");
		throw; return false;
	}
	png_init_io(png_ptr, fp);
	png_set_sig_bytes(png_ptr, 8);
	
	// read the PNG file
	if ( is_nbo() )
	{
		png_read_png(png_ptr, info_ptr, PNG_TRANSFORM_IDENTITY, png_voidp_NULL);
	}
	else
	{
		png_read_png(png_ptr, info_ptr, PNG_TRANSFORM_SWAP_ENDIAN, png_voidp_NULL);
	}
	
	// IHDR
	png_uint_32 Ncols, Nrows;
	int bit_depth;
	int color_type;
	int interlace_method, compression_method, filter_method;
	if ( !png_get_IHDR(png_ptr, info_ptr, &Ncols, &Nrows, &bit_depth, &color_type, &interlace_method, &compression_method,
					   &filter_method) )
	{
		std::fclose(fp);
		png_destroy_read_struct(&png_ptr, &info_ptr, &end_info);
		debug("pngwrite","error reading header in PNG file '"+filename+"'");
		throw; return false;
	}
	png_uint_32 no_colors = png_uint_32(ipow(2,bit_depth));
	
	png_bytep tRNS_trans;
	int tRNS_num_trans;
	png_color_16p tRNS_trans_values;
	bool has_tRNS = false;
	png_uint_16 trans_value;
	if ( png_get_tRNS(png_ptr, info_ptr, &tRNS_trans, &tRNS_num_trans, &tRNS_trans_values) )
	{
		if ( tRNS_num_trans > 0 )
		{
			has_tRNS = true;
			trans_value = tRNS_trans_values[0].gray;
		}
	}
	
	// pCAL
	png_charp pCAL_purpose;
	png_int_32 X0 = 0;
	png_int_32 X1 = no_colors - 1;
	int pCAL_type = PNG_EQUATION_LINEAR;
	int pCAL_nparams = 2;
	png_charp pCAL_units;
	png_charpp pCAL_params;
	T p0 = T(0);
	T p1 = T(0);
	std::string pCAL_param;
	if ( png_get_pCAL(png_ptr, info_ptr, &pCAL_purpose, &X0, &X1, &pCAL_type, &pCAL_nparams, &pCAL_units, &pCAL_params) )
	{
		switch(pCAL_type)
		{
			case PNG_EQUATION_LINEAR:
				p0 = IEEEton(std::string(pCAL_params[0]));
				p1 = IEEEton(std::string(pCAL_params[1]));
				break;
			default:
				break;
		}
	}
	
	// sCAL
	int sCAL_unit;
	double sCAL_wt = 0.0;
	double sCAL_ht = 0.0;
	png_get_sCAL(png_ptr, info_ptr, &sCAL_unit, &sCAL_wt, &sCAL_ht);
	sCAL_width = sCAL_wt;
	sCAL_height = sCAL_ht;
	
	// tEXT
	int num_text;
	png_textp text_ptr;
	std::vector< std::valarray<T> > data(0);
	angles.resize(0);
	png_uint_32 c_cols = 1;
	png_uint_32 c_rows = 1;
	bool IEEEdata = false;
	if ( png_get_text(png_ptr, info_ptr, &text_ptr, &num_text) )
	{
		for (size_t i=0; i<num_text; i++)
		{
			std::string key = text_ptr[i].key;
			std::string text = text_ptr[i].text;
			if ( key == std::string("Title") )
			{
				dataname = text;
			}
			else if ( key == std::string("collage") )
			{
				size_t ds = text.find(":");
				long int ndc = stoi(text.substr(0,ds));
				data.resize(ndc);
				ds++;
				text = text.substr(ds, text.size()-ds);
				ds = text.find(",");
				c_cols = stoi(text.substr(0,ds));
				ds++;
				text = text.substr(ds, text.size()-ds);
				c_rows = stoi(text);
			}
			else if ( key == std::string("raxis") )
			{
				if ( text == std::string("X") )
				{
					raxis = Angle::X;
				}
				else
				{
					raxis = Angle::Y;
				}
			}
		}
		iNcols = ( Ncols + 1 - c_cols ) / c_cols;
		iNrows = ( Nrows + 1 - c_rows ) / c_rows;
		if ( data.size() )
		{
			for (size_t i=0; i<num_text; i++)
			{
				std::string key = text_ptr[i].key;
				std::string text = text_ptr[i].text;
				if ( key == std::string("angles") )
				{
					if ( text == std::string("default") )
					{
						continue;
					}
					angles.resize(data.size());
					angles = IEEEtov<T>(text);
				}
				if ( key.substr(0,5) == std::string("data_") )
				{
					size_t dj = stoi(key.substr(5,key.size()-5));
					data[dj-1].resize(iNrows*iNcols);
					data[dj-1] = IEEEtov<T>(text);
					IEEEdata = true;
				}
			}
		}
	}
	
	if ( !IEEEdata )
	{
		// get the row_pointers
		png_bytepp row_pointers = png_get_rows(png_ptr, info_ptr);
		
		// initialize bit data
		std::valarray<png_uint_16> idata(Nrows*Ncols);
		for (size_t i=0; i<Nrows*Ncols; i++)
		{
			idata[i] = 0;
		}
		
		// create png_bytep pointers into idata
		std::vector<png_bytep> cidata(Nrows*Ncols*sizeof(png_uint_16));
		cidata[0] = (png_bytep)(&(idata[0]));
		for (size_t i=1; i<cidata.size(); i++)
		{
			cidata[i] = (cidata[i-1])+1;
		}
		for (size_t i=0; i<Nrows; i++)
		{
			for (size_t j=0; j<Ncols*sizeof(png_uint_16); j++)
			{
				*(cidata[i*Ncols*sizeof(png_uint_16)+j]) = row_pointers[i][j];
			}
		}
		
		// convert from integer to T
		for (size_t i=0; i<data.size(); i++)
		{
			data[i].resize(iNcols*iNrows);
		}
		for (png_uint_32 i=0; i<c_rows; i++)
		{
			for (png_uint_32 j=0; j<c_cols; j++)
			{
				size_t bd = i*c_cols+j;
				if ( bd < data.size() )
				{
					png_uint_32 bcol = j * iNcols + j;
					png_uint_32 brow = i * iNrows + i;
					for (png_uint_32 ii=0; ii<iNrows; ii++)
					{
						for (png_uint_32 jj=0; jj<iNcols; jj++)
						{
							size_t idx = (brow+ii) * Ncols + (bcol+jj);
							if ( has_tRNS && idata[idx] == trans_value )
							{
								data[bd][ii*iNcols+jj] = T(0);
							}
							else
							{
								data[bd][ii*iNcols+jj] = p0 + p1 * idata[idx] / ( X1 - X0 );
							}
						}
					}
				}
			}
		}
	}
	
	png_destroy_read_struct(&png_ptr, &info_ptr, &end_info);
	std::fclose(fp);
	
	output_data.resize(data.size()*iNrows*iNcols);
	for (size_t i=0; i<data.size(); i++)
	{
		output_data[std::slice(i*iNrows*iNcols, iNrows*iNcols, 1)] = data[i];
	}
	return true;
}
/* ----------------------------- */


#endif /* _IMGLIB_ */
