/**
 tomography
 
 A collection of C++ routines to handle tomographic functions.
 
 Rado Faletic
 Department of Physics
 Faculty of Science
 Australian National University  ACT  0200
 Australia
 
 Rado.Faletic@anu.edu.au
 
 21st April 2005
 22nd April 2022, updated to C++20 and updated libpng requirements
 */





#ifndef _TOMOGRAPHY_
#define _TOMOGRAPHY_





/* ---------- standard header files ---------- */
#include <bit>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <ctime>
#include <limits>
#include <png.h>
#include <string>
#include <valarray>
#include <vector>





/* ---------- user header files ---------- */
#include "angles.h"
#include "conversions.h"
#include "extra_math.h"
#include "front-end.h"
#include "grid.h"
#include "line.h"
#include "rotation.h"
#include "sparse_matrix.h"





/* ---------- Tomography namespace ---------- */

namespace Tomography
{

/// the engine for performing tomographic projections
template<class T> void projection(grid<T>&,                        // grid
                                  const std::vector< line<T> >&,   // rays
                                  std::valarray<bool>&,            // blanks
                                  std::valarray<T>&,               // projections
                                  SparseMatrix<T>&,                // Frechet derivative
                                  const grid_proj_method& = brute, // projection method
                                  const bool& = true,              // return projections
                                  const bool& = true,              // return derivatives
                                  const bool& = true);             // save projections

/// write a series of tomographic projections to a PNG image file
template<class T> bool pngwrite(const std::string&,                            // filename
                                const std::size_t&, const std::size_t&,        // Nrows, Ncols
                                const std::valarray<T>&,                       // data
                                const Angle::axes& = Angle::ZX,                // rotation axis
                                const std::valarray<T>& = std::valarray<T>(0), // angles
                                const T& = T(1),                               // scale
                                const bool& = true,                            // write real data
                                const bool& = false);                          // transparency

/// read a series of tomographic projections from a PNG image file
template<class T> bool pngread(const std::string&,         // filename
                               std::size_t&, std::size_t&, // Nrows, Ncols
                               std::valarray<T>&,          // data
                               std::valarray<bool>&,       // blanks
                               Angle::axes&,               // rotation axis
                               std::valarray<T>&,          // angles
                               T&,                         // scale
                               bool&);                     // real values?

/// create rays from a real set of tomograms
template<class T> std::vector< line<T> > create_rays(const std::size_t&, const std::size_t&,        // Nrows, Ncols
                                                     const std::size_t&,                            // number of datum
                                                     const Angle::axes& = Angle::XY,                // rotation axis
                                                     const std::valarray<T>& = std::valarray<T>(0), // angles
                                                     const T& = T(1));                              // scale

/// generate a synthetic set of rays
template<class T> std::vector< line<T> > synthetic_rays(const std::size_t&,          // number of angles
                                                        const std::size_t&,          // mesh size
                                                        const std::valarray<T>&,     // min
                                                        const std::valarray<T>&,     // max
                                                        std::size_t&, std::size_t&); // Nrows, Ncols
}





/* ---------- function definitions ---------- */





/* ---------- projection ---------- */
template<class T> void Tomography::projection(grid<T>& g,
                                              const std::vector< line<T> >& ray,
                                              std::valarray<bool>& blanks,
                                              std::valarray<T>& out_b,
                                              SparseMatrix<T>& A,
                                              const grid_proj_method& the_method,
                                              const bool& residuals_on,
                                              const bool& matrix_on,
                                              const bool& grid_samples_on)
{
    std::vector<T> b(0);
    if ( blanks.size() != ray.size() )
    {
        std::valarray<bool> obl = blanks;
        blanks.resize(ray.size(), true);
        if ( obl.size() <= blanks.size() )
        {
            blanks[std::slice(0,obl.size(),1)] = obl;
        }
        else
        {
            blanks = obl[std::slice(0,blanks.size(),1)];
        }
    }
    if ( !matrix_on )
    {
        A.clear();
    }
    
    for (std::size_t i=0; i<ray.size(); i++)
    {
        counter("projection", ray.size(), i + 1);
        if ( !blanks[i] )
        {
            continue;
        }
        
        std::vector< two_numbers<T> > projections(0);
        
        T bi = g.raytrace(ray[i], projections, the_method);
        if ( !projections.size() )
        {
            blanks[i] = false;
            continue;
        }
        else
        {
            if ( residuals_on )
            {
                b.push_back(bi);
            }
            if ( matrix_on )
            {
                A.AddRow(projections);
            }
            if ( grid_samples_on )
            {
                for (std::size_t j=0; j<projections.size(); j++)
                {
                    g[projections[j].n_integer] += projections[j].n_real;
                }
            }
        }
    }
    out_b.resize(b.size());
    std::copy(b.begin(), b.end(), &(out_b[0]));
}





/* ---------- pngwrite ---------- */
/// NOTE: input data should be in standard numerical format, ie starting at the bottom left. This function automatically converts this to image format, ie starting at the top left.
template<class T> bool Tomography::pngwrite(const std::string& myfilename,
                                            const std::size_t& iNrows, const std::size_t& iNcols,
                                            const std::valarray<T>& input_data_in,
                                            const Angle::axes& raxis,
                                            const std::valarray<T>& my_angles,
                                            const T& scale,
                                            const bool& IEEEdata,
                                            const bool& tpcy)
{
    int bit_depth = 16; // 16bit grayscale
    png_uint_32 no_colors = png_uint_32(ipow(2,bit_depth));
    png_int_32 X0 = ( tpcy ) ? 0 : 1;
    png_int_32 X1 = no_colors - 1;
    std::valarray<T> input_data = input_data_in;
    if ( !IEEEdata ) // stretch to use the full dynamic range in the integer values
    {
        T min = input_data.min();
        T max = input_data.max();
        T m = T( X1 - X0 ) / ( max - min );
        T c = T(X1) - m * max;
        input_data *= m;
        input_data += c;
    }
    std::vector< std::valarray<T> > data(0);
    
    // the configuration for the collage of images
    png_uint_32 c_cols = 1;
    png_uint_32 c_rows = 1;
    
    if ( input_data.size() == iNrows * iNcols )
    {
        data.push_back(input_data);
    }
    else
    {
        for (std::size_t i=0; i<input_data.size()/(iNrows*iNcols); i++)
        {
            std::valarray<T> tdata = input_data[std::slice(i*iNrows*iNcols,iNrows*iNcols,1)];
            data.push_back(tdata);
        }
        c_cols = data.size();
        std::size_t wh_sqr = (c_cols*iNcols+c_cols-1)*(c_cols*iNcols+c_cols-1) + (c_rows*iNrows+c_rows-1)*(c_rows*iNrows+c_rows-1);
        for (std::size_t i=1; i<data.size()+1; i++)
        {
            for (std::size_t j=1; j<data.size()+1; j++)
            {
                if ( i * j < data.size() )
                {
                    continue;
                }
                std::size_t wh_width = i*iNcols+i-1;
                std::size_t wh_height = j*iNrows+j-1;
                if ( wh_width < wh_height )
                {
                    continue;
                }
                std::size_t wh_sqr_n = wh_width*wh_width + wh_height*wh_height;
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
    }
    
    // convert from numerical format to image format
    for (std::size_t d=0; d<data.size(); d++)
    {
        for (std::size_t j=0; j<iNrows/2; j++)
        {
            for (std::size_t i=0; i<iNcols; i++)
            {
                std::swap(data[d][j*iNcols+i], data[d][(iNrows-1-j)*iNcols+i]);
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
    }
    if ( angles.size() != data.size() && angles.size() == iNcols )
    {
        use_angles = true;
    }
    
    // set up the PNG file and structures
    std::string filename = myfilename;
    std::string filetitle = myfilename;
    if ( filename.substr(filename.size()-4,4) == std::string(".PNG") )
    {
        filename = filename.substr(0,filename.size()-3) + "png";
        filetitle = filename.substr(0,filename.size()-4);
    }
    else if ( filename.substr(filename.size()-4,4) != std::string(".png") )
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
        return false;
    }
    png_structp png_ptr =
    png_create_write_struct(PNG_LIBPNG_VER_STRING, nullptr, (png_error_ptr)nullptr, (png_error_ptr)nullptr);
    if ( !png_ptr )
    {
        std::fclose(fp);
        return false;
    }
    png_infop info_ptr = png_create_info_struct(png_ptr);
    if ( !info_ptr )
    {
        std::fclose(fp);
        png_destroy_write_struct(&png_ptr, (png_infopp)nullptr);
        return false;
    }
    if ( setjmp(png_jmpbuf(png_ptr)) )
    {
        std::fclose(fp);
        png_destroy_write_struct(&png_ptr, &info_ptr);
        return false;
    }
    png_init_io(png_ptr, fp);
    
    // IHDR
    png_set_IHDR(png_ptr, info_ptr, Ncols, Nrows, bit_depth,
                 PNG_COLOR_TYPE_GRAY, PNG_INTERLACE_NONE, PNG_COMPRESSION_TYPE_DEFAULT, PNG_FILTER_TYPE_DEFAULT);
    
    // tRNS
    png_bytep trans = NULL;
    png_uint_16 num_trans = 1;
    png_color_16 trans_values[1];
    trans_values[0].gray = 0;
    png_set_tRNS(png_ptr, info_ptr, trans, num_trans, trans_values);
    
    // pCAL
    int nparams = 2; // real_data = p0 + p1 * pixel_value / ( X1 - X0 );
    png_charp units = (char*)"px";
    T dmax = input_data.max();
    T dmin = input_data.min();
    T p0 = dmin - ( dmax - dmin ) * T(X0) / T( X1 - X0 );
    T p1 = dmax - dmin;
    std::string p0s = ntoIEEE(p0);
    std::string p1s = ntoIEEE(p1);
    png_charp params[2] = { (char*)p0s.c_str(), (char*)p1s.c_str() };
    png_set_pCAL(png_ptr, info_ptr, "projection", X0, X1, PNG_EQUATION_LINEAR, nparams, units, params);
    
    // sCAL
    if ( !is_zero(scale) )
    {
        png_set_sCAL(png_ptr, info_ptr, PNG_SCALE_METER, scale, scale);
    }
    
    // tIME
    png_time mod_time;
    png_convert_from_time_t(&mod_time, std::time(0));
    png_set_tIME(png_ptr, info_ptr, &mod_time);
    
    // necessities for IEEEdata
    std::vector<std::string> IEEEs(0);
    std::vector<std::string> IEEEsd(0);
    std::string sdn = "";
    if ( IEEEdata )
    {
        std::size_t nodier = 1600; // number of data in each row
        
        std::valarray<T> IEdata(data.size()*iNrows*iNcols);
        for (std::size_t i=0; i<data.size(); i++)
        {
            IEdata[std::slice(i*iNrows*iNcols, iNrows*iNcols, 1)] = data[i];
        }
        for (std::size_t i=0; i<IEdata.size(); i+=nodier)
        {
            std::size_t dvs = i + nodier;
            dvs = ( IEdata.size() < dvs ) ? IEdata.size() - i : nodier;
            std::valarray<T> data_v = IEdata[std::slice(i, dvs, 1)];
            IEEEs.push_back(vtoIEEE(data_v));
        }
        IEdata.resize(0);
        
        IEEEsd.resize(IEEEs.size());
        for (std::size_t i=0; i<IEEEsd.size(); i++)
        {
            IEEEsd[i] = "IEEEdata" + std::to_string(i);
        }
        
        sdn = std::to_string(IEEEs.size());
    }
    
    // tEXT
    int num_text = ( IEEEdata ) ? 11 + IEEEs.size() : 10;
    png_text text_ptr[num_text];
    text_ptr[0].key = (char*)"Title";
    text_ptr[0].text = (char*)filetitle.c_str();
    text_ptr[0].text_length = filetitle.size();
    text_ptr[0].compression = PNG_TEXT_COMPRESSION_NONE;
    text_ptr[1].key = (char*)"Author";
    text_ptr[1].text = (char*)"Rado Faletic";
    text_ptr[1].text_length = 12;
    text_ptr[1].compression = PNG_TEXT_COMPRESSION_NONE;
    text_ptr[2].key = (char*)"Description";
    text_ptr[2].text = (char*)"tomographic projection";
    text_ptr[2].text_length = 22;
    text_ptr[2].compression = PNG_TEXT_COMPRESSION_NONE;
    std::string png_cyear = "(c) 2002, ";
    png_cyear += std::to_string(mod_time.year);
    png_cyear += " Rado Faletic";
    text_ptr[3].key = (char*)"Copyright";
    text_ptr[3].text = (char*)png_cyear.c_str();
    text_ptr[3].text_length = png_cyear.size();
    text_ptr[3].compression = PNG_TEXT_COMPRESSION_NONE;
    text_ptr[4].key = (char*)"Creation Time";
    text_ptr[4].text = (char*)"                             ";
    png_convert_to_rfc1123_buffer(text_ptr[4].text, &mod_time);
    text_ptr[4].text_length = std::strlen(text_ptr[4].text);
    text_ptr[4].compression = PNG_TEXT_COMPRESSION_NONE;
    text_ptr[5].key = (char*)"Software";
    text_ptr[5].text = (char*)"tomography.h";
    text_ptr[5].text_length = 12;
    text_ptr[5].compression = PNG_TEXT_COMPRESSION_NONE;
    text_ptr[6].key = (char*)"Source";
    text_ptr[6].text = (char*)"data";
    text_ptr[6].text_length = 4;
    text_ptr[6].compression = PNG_TEXT_COMPRESSION_NONE;
    std::string png_collage_txt = std::to_string(data.size());
    png_collage_txt += ":";
    png_collage_txt += std::to_string(c_cols);
    png_collage_txt += ",";
    png_collage_txt += std::to_string(c_rows);
    text_ptr[7].key = (char*)"collage";
    text_ptr[7].text = (char*)png_collage_txt.c_str();
    text_ptr[7].text_length = png_collage_txt.size();
    text_ptr[7].compression = PNG_TEXT_COMPRESSION_NONE;
    text_ptr[8].key = (char*)"raxis";
    text_ptr[8].text = ( raxis == Angle::Y || raxis == Angle::ZX ) ? (char*)"Y" : (char*)"X";
    text_ptr[8].text_length = 1;
    text_ptr[8].compression = PNG_TEXT_COMPRESSION_NONE;
    std::string png_angles_txt = "default";
    if ( use_angles )
    {
        png_angles_txt = vtoIEEE(angles);
    }
    text_ptr[9].key = (char*)"angles";
    text_ptr[9].text = (char*)png_angles_txt.c_str();
    text_ptr[9].text_length = png_angles_txt.size();
    text_ptr[9].compression = PNG_TEXT_COMPRESSION_NONE;
    if ( IEEEdata )
    {
        text_ptr[10].key = (char*)"ndatarows"; // stores how many data arrays there are
        text_ptr[10].text = (char*)sdn.c_str();
        text_ptr[10].text_length = sdn.size();
        text_ptr[10].compression = PNG_TEXT_COMPRESSION_NONE;
        for (std::size_t i=0; i<IEEEs.size(); i++)
        {
            text_ptr[11+i].key = (char*)IEEEsd[i].c_str();
            text_ptr[11+i].text = (char*)IEEEs[i].c_str();
            text_ptr[11+i].text_length = IEEEs[i].size();
            text_ptr[11+i].compression = PNG_TEXT_COMPRESSION_zTXt;
        }
    }
    png_set_text(png_ptr, info_ptr, text_ptr, num_text);
    
    // convert from T to integer
    std::valarray<png_uint_16> idata(png_uint_16(0),Nrows*Ncols);
    T m = p1 / T( X1 - X0 );
    for (png_uint_32 i=0; i<c_rows; i++)
    {
        for (png_uint_32 j=0; j<c_cols; j++)
        {
            std::size_t bd = i*c_cols+j;
            if ( bd < data.size() )
            {
                png_uint_32 bcol = j * iNcols + j;
                png_uint_32 brow = i * iNrows + i;
                for (png_uint_32 ii=0; ii<iNrows; ii++)
                {
                    for (png_uint_32 jj=0; jj<iNcols; jj++)
                    {
                        T tmp = (data[bd][ii*iNcols+jj] - p0) / m;
                        png_uint_16 itmp = ( tmp + std::numeric_limits<T>::round_error() > std::ceil(tmp) )
                        ? png_uint_16(std::ceil(tmp)) : png_uint_16(std::floor(tmp));
                        if ( tmp < T(0) )
                        {
                            itmp = 0;
                        }
                        else if ( 0 < itmp && itmp < X0 )
                        {
                            itmp = X0;
                        }
                        else if ( X1 < itmp )
                        {
                            itmp = X1;
                        }
                        idata[(brow+ii) * Ncols + (bcol+jj)] = itmp;
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
    if (std::endian::native == std::endian::big)
    {
        png_write_png(png_ptr, info_ptr, PNG_TRANSFORM_IDENTITY, nullptr);
    }
    else
    {
        png_write_png(png_ptr, info_ptr, PNG_TRANSFORM_SWAP_ENDIAN, nullptr);
    }
    
    png_destroy_write_struct(&png_ptr, &info_ptr);
    std::fclose(fp);
    return true;
}





/* ---------- pngread ---------- */
/// NOTE: this function converts from image format, ie starting at the top left, into standard numerical format, ie starting at the bottom left.
template<class T> bool Tomography::pngread(const std::string& myfilename,
                                           std::size_t& iNrows, std::size_t& iNcols,
                                           std::valarray<T>& output_data,
                                           std::valarray<bool>& blanks,
                                           Angle::axes& raxis,
                                           std::valarray<T>& angles,
                                           T& scale,
                                           bool& IEEEdata)
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
        return false;
    }
    png_byte sig2check[8];
    std::fread(sig2check, 1, 8, fp);
    if ( png_sig_cmp(sig2check, 0, 8) )
    {
        std::fclose(fp);
        return false;
    }
    
    png_structp png_ptr =
    png_create_read_struct(PNG_LIBPNG_VER_STRING, nullptr, (png_error_ptr)nullptr, (png_error_ptr)nullptr);
    if ( !png_ptr )
    {
        std::fclose(fp);
        return false;
    }
    
    png_infop info_ptr = png_create_info_struct(png_ptr);
    if ( !info_ptr )
    {
        std::fclose(fp);
        png_destroy_read_struct(&png_ptr, (png_infopp)nullptr, (png_infopp)nullptr);
        return false;
    }
    
    png_infop end_info = png_create_info_struct(png_ptr);
    if ( !end_info )
    {
        std::fclose(fp);
        png_destroy_read_struct(&png_ptr, &info_ptr, (png_infopp)nullptr);
        return false;
    }
    
    if ( setjmp(png_jmpbuf(png_ptr)) )
    {
        std::fclose(fp);
        png_destroy_read_struct(&png_ptr, &info_ptr, &end_info);
        return false;
    }
    png_init_io(png_ptr, fp);
    png_set_sig_bytes(png_ptr, 8);
    
    // read the PNG file
    if (std::endian::native == std::endian::big)
    {
        png_read_png(png_ptr, info_ptr, PNG_TRANSFORM_IDENTITY, nullptr);
    }
    else
    {
        png_read_png(png_ptr, info_ptr, PNG_TRANSFORM_SWAP_ENDIAN, nullptr);
    }
    
    // IHDR
    png_uint_32 Ncols, Nrows;
    int bit_depth;
    int color_type;
    int interlace_method, compression_method, filter_method;
    if ( !png_get_IHDR(png_ptr, info_ptr, &Ncols, &Nrows, &bit_depth,
                       &color_type, &interlace_method,
                       &compression_method, &filter_method) )
    {
        std::fclose(fp);
        png_destroy_read_struct(&png_ptr, &info_ptr, &end_info);
        return false;
    }
    png_uint_32 no_colors = png_uint_32(ipow(2,bit_depth));
    
    // tRNS
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
    T p1 = T(X1);
    if ( png_get_pCAL(png_ptr, info_ptr, &pCAL_purpose, &X0, &X1, &pCAL_type, &pCAL_nparams, &pCAL_units, &pCAL_params) )
    {
        switch(pCAL_type)
        {
            case PNG_EQUATION_LINEAR:
                p0 = std::stod(std::string(pCAL_params[0]));
                p1 = std::stod(std::string(pCAL_params[1]));
                break;
            default:
                break;
        }
    }
    
    // sCAL
    int sCAL_unit;
    double sCAL_width = 1.0;
    double sCAL_height = 1.0;
    png_get_sCAL(png_ptr, info_ptr, &sCAL_unit, &sCAL_width, &sCAL_height);
    scale = 0.5 * (sCAL_width + sCAL_height);
    if ( is_zero(scale) )
    {
        scale = T(1);
    }
    
    // tEXT
    int num_text;
    png_textp text_ptr;
    std::vector< std::valarray<T> > data(0);
    angles.resize(0);
    png_uint_32 c_cols = 1;
    png_uint_32 c_rows = 1;
    IEEEdata = false;
    if ( png_get_text(png_ptr, info_ptr, &text_ptr, &num_text) )
    {
        for (std::size_t i=0; i<num_text; i++)
        {
            std::string key = text_ptr[i].key;
            if ( key == std::string("collage") )
            {
                std::string text = text_ptr[i].text;
                std::size_t ds = text.find(":");
                long int ndc = std::stol(text.substr(0,ds));
                data.resize(ndc);
                ds++;
                text = text.substr(ds, text.size()-ds);
                ds = text.find(",");
                c_cols = std::stol(text.substr(0,ds));
                ds++;
                text = text.substr(ds, text.size()-ds);
                c_rows = std::stol(text);
            }
            else if ( key == std::string("raxis") )
            {
                std::string text = text_ptr[i].text;
                if ( text == std::string("Y") )
                {
                    raxis = Angle::Y;
                }
                else
                {
                    raxis = Angle::X;
                }
            }
        }
        iNcols = ( Ncols + 1 - c_cols ) / c_cols;
        iNrows = ( Nrows + 1 - c_rows ) / c_rows;
        if ( !data.size() )
        {
            data.resize(1);
        }
        std::vector<std::string> IEEEs(0);
        for (std::size_t i=0; i<num_text; i++)
        {
            std::string key = text_ptr[i].key;
            if ( key == std::string("angles") )
            {
                std::string text = text_ptr[i].text;
                if ( text == std::string("default") )
                {
                    continue;
                }
                angles.resize(data.size());
                angles = IEEEtov<T>(text);
            }
            if ( key == std::string("ndatarows") )
            {
                std::string text = text_ptr[i].text;
                IEEEs.resize(std::stol(text));
                IEEEdata = true;
            }
        }
        if ( IEEEdata )
        {
            output_data.resize(data.size()*iNrows*iNcols);
            for (std::size_t i=0; i<num_text; i++)
            {
                std::string key = text_ptr[i].key;
                if ( key.substr(0,8) == "IEEEdata" )
                {
                    std::size_t data_n = std::stol(key.substr(8,key.size()-8));
                    IEEEs[data_n] = std::string(text_ptr[i].text);
                }
            }
            std::size_t data_i = 0;
            for (std::size_t i=0; i<IEEEs.size(); i++)
            {
                std::valarray<T> data_v = IEEEtov<T>(IEEEs[i]);
                output_data[std::slice(data_i,data_v.size(),1)] = data_v;
                data_i += data_v.size();
            }
        }
    }
    if ( c_cols == 1 && c_rows == 1 )
    {
        iNcols = Ncols;
        iNrows = Nrows;
    }
    if ( !data.size() )
    {
        data.resize(1);
    }
    
    blanks.resize(data.size()*iNcols*iNrows, true);
    if ( !IEEEdata || has_tRNS )
    {
        // bit data
        std::valarray<png_uint_16> idata(png_uint_16(0), Nrows*Ncols);
        // the row pointers
        png_bytepp row_pointers = png_get_rows(png_ptr, info_ptr);
        std::valarray<png_uint_16*> idatap(Nrows*Ncols);
        for (std::size_t i=0; i<Nrows; i++)
        {
            for(std::size_t j=0; j<Ncols; j++)
            {
                if ( bit_depth == 16 )
                {
                    idatap[i*Ncols+j] = (png_uint_16*)(row_pointers[i]+2*j);
                    idata[i*Ncols+j] = *(idatap[i*Ncols+j]);
                }
                else if ( bit_depth == 8 )
                {
                    idata[i*Ncols+j] = *(row_pointers[i]+j);
                }
            }
        }
        // convert from integer to T
        for (std::size_t i=0; i<data.size(); i++)
        {
            data[i].resize(iNcols*iNrows);
        }
        T m = p1 / T( X1 - X0 );
        for (png_uint_32 i=0; i<c_rows; i++)
        {
            for (png_uint_32 j=0; j<c_cols; j++)
            {
                std::size_t bd = i*c_cols+j;
                if ( bd >= data.size() )
                {
                    continue;
                }
                png_uint_32 bcol = j * iNcols + j;
                png_uint_32 brow = i * iNrows + i;
                for (png_uint_32 ii=0; ii<iNrows; ii++)
                {
                    for (png_uint_32 jj=0; jj<iNcols; jj++)
                    {
                        std::size_t idx = (brow+ii) * Ncols + (bcol+jj);
                        if ( has_tRNS && idata[idx] == trans_value )
                        {
                            blanks[bd*iNrows*iNcols+ii*iNcols+jj] = false;
                            data[bd][ii*iNcols+jj] = T(0);
                        }
                        else
                        {
                            data[bd][ii*iNcols+jj] = p0 + m * T(idata[idx]);
                        }
                    }
                }
            }
        }
        if ( !IEEEdata )
        {
            output_data.resize(data.size()*iNrows*iNcols);
            for (std::size_t i=0; i<data.size(); i++)
            {
                output_data[std::slice(i*iNrows*iNcols, data[i].size(), 1)] = data[i];
            }
        }
    }
    
    png_destroy_read_struct(&png_ptr, &info_ptr, &end_info);
    std::fclose(fp);
    
    // convert from image format to numerical format
    for (std::size_t d=0; d<output_data.size(); d+=(iNrows*iNcols))
    {
        for (std::size_t j=0; j<iNrows/2; j++)
        {
            for (std::size_t i=0; i<iNcols; i++)
            {
                std::swap(output_data[d+j*iNcols+i], output_data[d+(iNrows-1-j)*iNcols+i]);
                std::swap(blanks[d+j*iNcols+i], blanks[d+(iNrows-1-j)*iNcols+i]);
            }
        }
    }
    return true;
}





/* ---------- create_rays ---------- */
template<class T> std::vector< line<T> > Tomography::create_rays(const std::size_t& Nrows, const std::size_t& Ncols,
                                                                 const std::size_t& data_size,
                                                                 const Angle::axes& axis,
                                                                 const std::valarray<T>& input_angles,
                                                                 const T& scale)
{
    unsigned short dim = ( Nrows == 1 || Ncols == 1 ) ? 2 : 3;
    
    std::valarray<T> angles = input_angles;
    if ( !angles.size() )
    {
        angles.resize(data_size / ( Nrows * Ncols ));
        for (std::size_t i=0; i<angles.size(); i++)
        {
            angles[i] = T(i*180) / T(angles.size());
        }
    }
    
    std::vector< std::valarray<T> > ipoints(Nrows*Ncols, std::valarray<T>(T(0),dim));
    T max_height = ( T( Nrows - 1 ) / T(2) ) * scale;
    for (std::size_t i=0; i<Nrows; i++)
    {
        T y_value = max_height - i * scale;
        for (std::size_t j=0; j<Ncols; j++)
        {
            ipoints[i*Ncols+j][0] = ( T(j) - T(Ncols-1)/T(2) ) * scale;
            ipoints[i*Ncols+j][1] = y_value;
        }
    }
    std::valarray<T> islope(T(0), dim);
    islope[islope.size()-1] = T(1);
    
    std::vector< line<T> > rays(data_size);
    
    Rotation<T> rot(dim);
    std::valarray<T> slope = islope;
    std::size_t count = 0;
    std::size_t total = angles.size() * ipoints.size();
    debugn("creating " + std::to_string(total) + " rays: ");
    for (std::size_t i=0; i<angles.size(); i++)
    {
        rot.reset(T(angles[i]), axis);
        slope = rot.O(islope);
        for (std::size_t j=0; j<ipoints.size(); j++)
        {
            count++;
            counter("ray", total, count);
            line<T> temp_line(slope, rot.O(ipoints[j]));
            rays.push_back(temp_line);
        }
    }
    debugn("\n");
    
    return rays;
}





/* ---------- synthetic_rays ---------- */
template<class T> std::vector< line<T> > Tomography::synthetic_rays(const std::size_t& nangles,
                                                                    const std::size_t& meshsize,
                                                                    const std::valarray<T>& vmin,
                                                                    const std::valarray<T>& vmax,
                                                                    std::size_t& Nrows, std::size_t& Ncols)
{
    unsigned short dim = vmin.size();
    
    T xdist = std::abs(vmax[0] - vmin[0]);
    T ydist = std::abs(vmax[1] - vmin[1]);
    T spacing = T(1);
    
    std::valarray<T> origin = ( vmin + vmax ) / T(2);
    std::vector< std::valarray<T> > ipoints;
    std::valarray<T> islope(T(0), dim);
    islope[islope.size()-1] = T(1);
    std::vector< line<T> > rays(0);
    
    T y_value = origin[1];
    switch (dim)
    {
        case 2:
            Nrows = 1;
            Ncols = meshsize;
            spacing = xdist / T(meshsize);
            ipoints.resize(Nrows*Ncols, std::valarray<T>(T(y_value),dim));
            for (std::size_t i=0; i<Nrows; i++)
            {
                for (std::size_t j=0; j<Ncols; j++)
                {
                    ipoints[j][0] = vmin[0] + j * spacing;
                }
            }
            break;
        case 3:
            if ( xdist > ydist )
            {
                Nrows = std::size_t( T(meshsize) * ydist / xdist );
                Ncols = meshsize;
                spacing = xdist / T(Ncols);
            }
            else
            {
                Nrows = meshsize;
                Ncols = std::size_t( T(meshsize) * xdist / ydist );
                spacing = ydist / T(Nrows);
            }
            ipoints.resize(Nrows*Ncols, std::valarray<T>(T(origin[2]),dim));
            for (std::size_t i=0; i<Nrows; i++)
            {
                y_value = vmax[1] - i * spacing;
                for (std::size_t j=0; j<Ncols; j++)
                {
                    ipoints[i*Ncols+j][0] = vmin[0] + j * spacing;
                    ipoints[i*Ncols+j][1] = y_value;
                }
            }
            break;
        default:
            Nrows = 0;
            Ncols = 0;
            return rays;
            break;
    }
    
    Rotation<T> rot(dim);
    rot.set_origin(origin);
    std::valarray<T> slope = islope;
    std::size_t count = 0;
    std::size_t total = nangles * ipoints.size();
    debugn("creating " + std::to_string(total) + " rays: ");
    for (std::size_t i=0; i<nangles; i++)
    {
        rot.reset(T(i*180)/T(nangles), Angle::ZX);
        slope = rot.O(islope);
        for (std::size_t j=0; j<ipoints.size(); j++)
        {
            count++;
            counter("ray", total, count);
            line<T> temp_line(slope, rot(ipoints[j]));
            rays.push_back(temp_line);
        }
    }
    debugn("\n");
    return rays;
}





#endif /* _TOMOGRAPHY_ */
