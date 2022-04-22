/**
 plot3d
 
 routines specific to plot3d grids.
 
 Rado Faletic
 Department of Physics
 Faculty of Science
 Australian National University  ACT  0200
 Australia
 
 Rado.Faletic@anu.edu.au
 14th July 2004
 19th April 2022, updated to C++20 and updated libpng requirements 
 */

/**
 NOTE: see `Plot3D.txt' or `Plot3D.html' for more information about this data/file format
 
 Plot3D formats referenced from:
 CFD-GEOM
 Interactive Geometric Modeling and Grid Generation Software
 Appendix C: Grid File Formats
 Version 4.0, Feburary 1998
 CFD Research Corporation
 */





#ifndef _PLOT3D_
#define _PLOT3D_





/* ---------- standard header files ---------- */
#include <algorithm>
#include <bit>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <ctime>
#include <fstream>
#include <png.h>
#include <string>
#include <valarray>
#include <vector>





/* ---------- user header files ---------- */
#include "conversions.h"
#include "extra_math.h"
#include "file.h"
#include "fortran_io.h"
#include "front-end.h"
#include "polygon.h"
#include "polyhedron.h"
#include "shape.h"
#include "segment.h"





/* ---------- class & function declaration ---------- */

namespace Plot3D
{

template<class T> bool read_interface(std::vector< std::valarray<T> >&,
                                      std::valarray<std::size_t>&,
                                      std::valarray<std::size_t>&,
                                      std::valarray<std::size_t>&,
                                      std::string&,
                                      dataformat&,
                                      dataprecision&,
                                      bool&);

template<class T> bool read(std::valarray<T>&,
                            std::valarray<T>&,
                            std::valarray<T>&,
                            std::valarray<bool>&,
                            std::valarray<std::size_t>&,
                            std::valarray<std::size_t>&,
                            std::valarray<std::size_t>&,
                            const std::string&,
                            const dataformat& = Unformatted,
                            const dataprecision& = Single,
                            const bool& = true,
                            const bool& = true);

template<class T> void read_data_interface(std::valarray<T>&,
                                           std::string&,
                                           std::valarray<std::size_t>&,
                                           std::valarray<std::size_t>&,
                                           std::valarray<std::size_t>&,
                                           std::string&,
                                           dataformat&,
                                           dataprecision&,
                                           bool&);

template<class T> bool read_data(std::valarray<T>&,
                                 std::valarray<T>&,
                                 std::valarray<T>&,
                                 std::valarray<T>&,
                                 std::valarray<T>&,
                                 const std::string&,
                                 const dataformat& = Unformatted,
                                 const dataprecision& = Single,
                                 const bool& = true);

template<class T> bool extract_data(std::valarray<T>&,
                                    const std::string&,
                                    const dataformat& = Unformatted,
                                    const dataprecision& = Single,
                                    const bool& = true,
                                    const unsigned short& = 1);

template<class T> bool write_interface(const std::vector< std::valarray<T> >&,
                                       const std::valarray<std::size_t>&,
                                       const std::valarray<std::size_t>&,
                                       const std::valarray<std::size_t>&,
                                       const std::string&);

template<class T> bool write(const std::valarray<T>&,
                             const std::valarray<T>&,
                             const std::valarray<T>&,
                             const std::valarray<bool>&,
                             const std::size_t&,
                             const std::size_t&,
                             const std::size_t&,
                             const std::string&,
                             const dataformat& = Unformatted);

template<class T> bool write(const std::valarray<T>&,
                             const std::valarray<T>&,
                             const std::valarray<T>&,
                             const std::valarray<bool>&,
                             const std::valarray<std::size_t>&,
                             const std::valarray<std::size_t>&,
                             const std::valarray<std::size_t>&,
                             const std::string&,
                             const dataformat& = Unformatted,
                             const dataprecision& = Single,
                             const bool& = true,
                             const bool& = true);

template<class T> bool write_data_interface(const std::valarray<T>&,
                                            const std::string&,
                                            const std::valarray<std::size_t>&,
                                            const std::valarray<std::size_t>&,
                                            const std::valarray<std::size_t>&,
                                            const std::string&);

template<class T> bool write_data(const std::valarray<T>&,
                                  const std::valarray<T>&,
                                  const std::valarray<T>&,
                                  const std::valarray<T>&,
                                  const std::valarray<T>&,
                                  const std::size_t&,
                                  const std::size_t&,
                                  const std::size_t&,
                                  const std::string&,
                                  const dataformat& = Unformatted);

template<class T> bool write_data(const std::valarray<T>&,
                                  const std::valarray<T>&,
                                  const std::valarray<T>&,
                                  const std::valarray<T>&,
                                  const std::valarray<T>&,
                                  const std::valarray<std::size_t>&,
                                  const std::valarray<std::size_t>&,
                                  const std::valarray<std::size_t>&,
                                  const std::string&,
                                  const dataformat& = Unformatted,
                                  const dataprecision& = Single,
                                  const bool& = true);

bool write_var(const std::string&,
               const std::string& = "Q1",
               const std::string& = "Q2",
               const std::string& = "Q3",
               const std::string& = "Q4",
               const std::string& = "Q5");

void get_vertices(const unsigned short&,
                  const std::valarray<std::size_t>&,
                  const std::valarray<std::size_t>&,
                  const std::valarray<std::size_t>&,
                  std::valarray<std::size_t>&,
                  const std::size_t&);

void get_surrounds(const unsigned short&,
                   const std::valarray<std::size_t>&,
                   const std::valarray<std::size_t>&,
                   const std::valarray<std::size_t>&,
                   std::valarray<std::size_t>&,
                   const std::size_t&);

template<class T> void n2c(const std::vector< std::valarray<T> >&,
                           const std::valarray<std::size_t>&,
                           const std::valarray<std::size_t>&,
                           const std::valarray<std::size_t>&,
                           const std::valarray<T>&,
                           std::valarray<T>&);

template<class T> void c2n(const std::vector< std::valarray<T> >&,
                           const std::valarray<std::size_t>&,
                           const std::valarray<std::size_t>&,
                           const std::valarray<std::size_t>&,
                           const std::valarray<T>&,
                           std::valarray<T>&);

template<class T> bool write_neighbours(const std::vector< std::valarray<T> >&,
                                        const std::valarray<std::size_t>&,
                                        const std::valarray<std::size_t>&,
                                        const std::valarray<std::size_t>&,
                                        const std::vector< std::valarray<std::size_t> >&,
                                        const std::string&,
                                        const dataformat& = Unformatted);

enum values { NODES, CENTROIDS };

template<class T> bool to_PNG(const std::valarray<T>&,
                              const std::valarray<T>&,
                              const std::valarray<T>&,
                              const std::valarray<std::size_t>&,
                              const std::valarray<std::size_t>&,
                              const std::valarray<std::size_t>&,
                              const std::valarray<T>&,
                              const std::string&,
                              const Plot3D::values& = NODES);

}





/* ---------- function definitions ---------- */





/* ---------- read_interface ---------- */
template<class T> bool Plot3D::read_interface(std::vector< std::valarray<T> >& nodes,
                                              std::valarray<std::size_t>& Nx,
                                              std::valarray<std::size_t>& Ny,
                                              std::valarray<std::size_t>& Nz,
                                              std::string& gfilename,
                                              dataformat& format,
                                              dataprecision& fsize,
                                              bool& multidomain)
{
    bool blanking;
    
    std::string filename;
    bool open = false;
    
    std::valarray<T> X, Y, Z;
    
    // get file details
    while ( !open )
    {
        std::string stuff = "";
        filename = "";
        while ( filename == std::string("") ||
               filename == std::string(" ") ||
               filename == std::string("\n") )
        {
            filename = question("what is the name of your Plot3D file (without the extension)", std::string(""));
        }
        
        stuff = question("is your Plot3D file Unformatted, Formatted or Binary", std::string("Unformatted"));
        if ( stuff == std::string("Formatted") ||
            stuff == std::string("formatted") ||
            stuff == std::string("FORMATTED") ||
            stuff == std::string("F") ||
            stuff == std::string("f") )
        {
            format = Formatted;
        }
        else if ( stuff == std::string("Binary") ||
                 stuff == std::string("binary") ||
                 stuff == std::string("BINARY") ||
                 stuff == std::string("B") ||
                 stuff == std::string("b") )
        {
            format = Binary;
        }
        else
        {
            format = Unformatted;
        }
        
        stuff = "";
        
        stuff = question("is your Plot3D file single or double precision", std::string("Single"));
        if ( stuff == std::string("double") ||
            stuff == std::string("Double") ||
            stuff == std::string("DOUBLE") ||
            stuff == std::string("d") ||
            stuff == std::string("D") )
        {
            fsize = Double;
        }
        else
        {
            fsize = Single;
        }
        
        stuff = "";
        
        multidomain = yesno("is your Plot3D file multidomain", true);
        
        blanking = yesno("does your Plot3D file use blanking", true);
        
        
        stuff = ".P";
        switch(format)
        {
            case Formatted:
                stuff += 'F';
                break;
            case Unformatted:
                stuff += 'U';
                break;
            case Binary:
                stuff += 'B';
                break;
        }
        switch(fsize)
        {
            case Double:
                stuff += "GD";
                break;
            default:
                stuff += "G";
                break;
        }
        std::valarray<bool> B;
        if ( Plot3D::read(X, Y, Z, B, Nx, Ny, Nz, filename+stuff, format, fsize, multidomain, blanking) )
        {
            open = true;
        }
        else
        {
            message("Could not open the file named `" + filename + stuff + "', please re-enter your information...");
            open = false;
        }
    }
    
    nodes.resize(X.size(), std::valarray<T>(3));
    for (std::size_t loop=0; loop<nodes.size(); loop++)
    {
        nodes[loop][0] = X[loop];
        nodes[loop][1] = Y[loop];
        nodes[loop][2] = Z[loop];
    }
    
    gfilename = filename;
    return true;
}





/* ---------- read ---------- */
template<class T> bool Plot3D::read(std::valarray<T>& X,
                                    std::valarray<T>& Y,
                                    std::valarray<T>& Z,
                                    std::valarray<bool>& B,
                                    std::valarray<std::size_t>& Nx,
                                    std::valarray<std::size_t>& Ny,
                                    std::valarray<std::size_t>& Nz,
                                    const std::string& filename,
                                    const dataformat& format,
                                    const dataprecision& fsize,
                                    const bool& multidomain,
                                    const bool& blanking)
{
    X.resize(0);
    Y.resize(0);
    Z.resize(0);
    B.resize(0);
    Nx.resize(0);
    Ny.resize(0);
    Nz.resize(0);
    
    std::ifstream gridfile;
    switch(format)
    {
        case Unformatted: case Binary:
            gridfile.open(filename.c_str(),std::ios_base::binary);
            break;
        default:
            gridfile.open(filename.c_str());
            break;
    }
    if ( !gridfile )
    {
        debug("Plot3D::read", "can't open Plot3D grid file `" + filename + "'");
        throw;
        return false;
    }
    
    // byte swapping facility
    bool byte_swapping = (std::endian::native == std::endian::little);
    
    // get information about number of domains
    std::size_t domains = 1;
    if ( multidomain )
    {
        Fortran::iread_line(gridfile, domains, byte_swapping, format);
        if ( domains < 1 )
        {
            debug("Plot3D::read", "bad number of domains");
            gridfile.close();
            throw;
            return false;
        }
    }
    
    // get information about number of nodes
    std::size_t isize;
    Fortran::read_byte_info(gridfile, isize, byte_swapping, format);
    if ( format == Unformatted && domains * 3 * 4 != isize )
    {
        debug("Plot3D::read", "error reading Plot3D file `" + filename + "'");
        gridfile.close();
        throw;
        return false;
    }
    Nx.resize(domains);
    Ny.resize(domains);
    Nz.resize(domains);
    std::size_t numnodes = 0;
    for (std::size_t loop=0; loop<domains; loop++)
    {
        Fortran::iread_item(gridfile, Nx[loop], byte_swapping, format);
        Fortran::iread_item(gridfile, Ny[loop], byte_swapping, format);
        Fortran::iread_item(gridfile, Nz[loop], byte_swapping, format);
        if ( Nx[loop] < 1 )
        {
            Nx[loop] = 1;
        }
        if ( Ny[loop] < 1 )
        {
            Ny[loop] = 1;
        }
        if ( Nz[loop] < 1 )
        {
            Nx[loop] = 1;
        }
        numnodes += Nx[loop] * Ny[loop] * Nz[loop];
    }
    std::size_t itemp;
    Fortran::read_byte_info(gridfile, itemp, byte_swapping, format);
    if ( itemp != isize )
    {
        debug("Plot3D::read", "inconsistent Plot3D file `" + filename + "'");
        Nx.resize(0);
        Ny.resize(0);
        Nz.resize(0);
        gridfile.close();
        throw;
        return false;
    }
    X.resize(numnodes);
    Y.resize(numnodes);
    Z.resize(numnodes);
    B.resize(numnodes);
    
    // read actual grid nodes
    std::size_t zero = 0;
    for (std::size_t loop=0; loop<domains; loop++)
    {
        Fortran::read_byte_info(gridfile, isize, byte_swapping, format);
        // I values
        for (std::size_t k=0; k<Nz[loop]; k++)
        {
            for (std::size_t j=0; j<Ny[loop]; j++)
            {
                for (std::size_t i=0; i<Nx[loop]; i++)
                {
                    Fortran::fread_item(gridfile, X[zero + k * Ny[loop] * Nx[loop] + j * Nx[loop] + i], byte_swapping, format, fsize);
                }
            }
        }
        // J values
        for (std::size_t k=0; k<Nz[loop]; k++)
        {
            for (std::size_t j=0; j<Ny[loop]; j++)
            {
                for (std::size_t i=0; i<Nx[loop]; i++)
                {
                    Fortran::fread_item(gridfile, Y[zero + k * Ny[loop] * Nx[loop] + j * Nx[loop] + i], byte_swapping, format, fsize);
                }
            }
        }
        // K values
        for (std::size_t k=0; k<Nz[loop]; k++)
        {
            for (std::size_t j=0; j<Ny[loop]; j++)
            {
                for (std::size_t i=0; i<Nx[loop]; i++)
                {
                    Fortran::fread_item(gridfile, Z[zero + k * Ny[loop] * Nx[loop] + j * Nx[loop] + i], byte_swapping, format,fsize);
                }
            }
        }
        // blanking
        if ( blanking )
        {
            for (std::size_t k=0; k<Nz[loop]; k++)
            {
                for (std::size_t j=0; j<Ny[loop]; j++)
                {
                    for (std::size_t i=0; i<Nx[loop]; i++)
                    {
                        Fortran::bread_item(gridfile, B[zero + k * Ny[loop] * Nx[loop] + j * Nx[loop] + i], byte_swapping, format);
                    }
                }
            }
        }
        Fortran::read_byte_info(gridfile, itemp, byte_swapping, format);
        if ( itemp != isize )
        {
            debug("Plot3D::read", "inconsistent Plot3D file `" + filename + "'");
            X.resize(0);
            Y.resize(0);
            Z.resize(0);
            B.resize(0);
            Nx.resize(0);
            Ny.resize(0);
            Nz.resize(0);
            gridfile.close();
            throw;
            return false;
        }
        zero += Nx[loop] * Ny[loop] * Nz[loop];
    }
    
    gridfile.close();
    return true;
}





/* ---------- read_data_interface ---------- */
/// reads Plot3D data (or Q-file)
template<class T> void Plot3D::read_data_interface(std::valarray<T>& data,
                                                   std::string& dataname,
                                                   std::valarray<std::size_t>& Nx,
                                                   std::valarray<std::size_t>& Ny,
                                                   std::valarray<std::size_t>& Nz,
                                                   std::string& gfilename,
                                                   dataformat& format,
                                                   dataprecision& fsize,
                                                   bool& multidomain)
{
    std::string qfilename, vfilename;
    std::ifstream qfile;
    std::ifstream vfile;
    unsigned short variable = 1;
    std::vector<std::string> vars(10,"Q");
    std::string variablename;
    std::size_t numvar;
    std::size_t ft = 10;
    bool open = false;
    bool qscalar = true;
    bool qvector = true;
    
    std::size_t i = 0;
    
    // initialise the `vars' vector, for holding names of variables
    for (i=0; i<vars.size(); i++)
    {
        switch(i)
        {
            case 0:
                vars[i] += '1';
                break;
            case 1:
                vars[i] += '2';
                break;
            case 2:
                vars[i] += '3';
                break;
            case 3:
                vars[i] += '4';
                break;
            case 4:
                vars[i] += '5';
                break;
            case 5:
                vars[i] += '6';
                break;
            case 6:
                vars[i] += '7';
                break;
            case 7:
                vars[i] += '8';
                break;
            case 8:
                vars[i] += '9';
                break;
            case 9:
                vars[i] += "10";
                break;
        }
    }
    
    // see if user has already opened the grid file, and if so try to open the corresponding Q-file
    if ( gfilename == std::string("") || gfilename == std::string(" ") )
    {
        open = false;
    }
    else
    {
        message("Existing grid called `" + gfilename + "', will attempt to use this filename");
        qfilename = gfilename;
        vfilename = gfilename + ".VAR";
        qfilename += ".P";
        switch(format)
        {
            case Formatted:
                qfilename += 'F';
                break;
            case Unformatted:
                qfilename += 'U';
                break;
            case Binary:
                qfilename += 'B';
                break;
        }
        qfilename += "S";
        if ( fsize == Double )
        {
            qfilename += "D";
        }
        switch(format)
        {
            case Unformatted: case Binary:
                qfile.open(qfilename.c_str(),std::ios_base::binary);
                break;
            default:
                qfile.open(qfilename.c_str());
                break;
        }
        if ( !qfile )
        {
            ft -= 5;
            qscalar = false;
            message("could not find scalar Q-file `" + qfilename + "'");
        }
        else
        {
            qscalar = true;
        }
        qfile.close();
        switch(fsize)
        {
            case Double:
                qfilename.resize(qfilename.size()-2);
                break;
            default:
                qfilename.resize(qfilename.size()-1);
                break;
        }
        qfilename += "V";
        if ( fsize == Double )
        {
            qfilename += "D";
        }
        switch(format)
        {
            case Unformatted: case Binary:
                qfile.open(qfilename.c_str(),std::ios_base::binary);
                break;
            default:
                qfile.open(qfilename.c_str());
                break;
        }
        if ( !qfile )
        {
            ft -= 5;
            qvector = false;
            message("could not find vector Q-file `" + qfilename + "'");
        }
        else
        {
            qvector = true;
        }
        qfile.close();
        vfile.open(vfilename.c_str());
        if ( vfile && ft != 0 )
        {
            numvar = 0;
            std::getline(vfile,vars[numvar]);
            while ( !vfile.eof() )
            {
                numvar++;
                if ( vars[numvar-1][0] == ' ' )
                {
                    for (i=0; i<vars[numvar-1].size()-1; i++)
                    {
                        vars[numvar-1][i] = vars[numvar-1][i+1];
                    }
                    vars[numvar-1].resize(vars[numvar-1].size()-1);
                }
                
                if ( vars[numvar-1] == std::string("rho") ||
                    vars[numvar-1] == std::string("density") ||
                    vars[numvar-1] == std::string("Rho") ||
                    vars[numvar-1] == std::string("Density") ||
                    vars[numvar-1] == std::string("RHO") ||
                    vars[numvar-1] == std::string("DENSITY") )
                {
                    variable = numvar;
                }
                if ( numvar + 1 > vars.size() )
                {
                    vars.resize(vars.size()+1);
                }
                std::getline(vfile, vars[numvar]);
            }
        }
        else if ( ft > 0 )
        {
            numvar = ft;
        }
        else
        {
            open = false;
        }
        vfile.close();
        if ( numvar != vars.size() )
        {
            vars.resize(numvar);
        }
        if ( ft > 0 )
        {
            std::string number_word;
            number_word = question("Please select a variable to extract", vars, vars[variable - 1]);
            i = variable;
            variable = 0;
            if ( number_word == vars[0] )
            {
                variable = 1;
            }
            else if ( number_word == vars[1] )
            {
                variable = 2;
            }
            else if ( number_word == vars[2] )
            {
                variable = 3;
            }
            else if ( number_word == vars[3] )
            {
                variable = 4;
            }
            else if ( number_word == vars[4] )
            {
                variable = 5;
            }
            else if ( ft > 5 )
            {
                if ( number_word == vars[5] )
                {
                    variable = 6;
                }
                else if ( number_word == vars[6] )
                {
                    variable = 7;
                }
                else if ( number_word == vars[7] )
                {
                    variable = 8;
                }
                else if ( number_word == vars[8] )
                {
                    variable = 9;
                }
                else if ( number_word == vars[9] )
                {
                    variable = 10;
                }
            }
            if ( variable == 0 )
            {
                variable = i;
            }
            if ( variable <= 5 && qscalar )
            {
                switch(fsize)
                {
                    case Double:
                        qfilename.resize(qfilename.size()-2);
                        break;
                    default:
                        qfilename.resize(qfilename.size()-1);
                        break;
                }
                qfilename += "S";
                if ( fsize == Double )
                {
                    qfilename += "D";
                }
            }
            message("will read variable " + std::to_string(variable) + " (" + vars[variable - 1] + ") from the Q-file `" + qfilename + "'");
            if ( variable > 5 )
            {
                variable -= 5;
            }
            if ( Plot3D::extract_data(data, qfilename, format, fsize, multidomain, variable) )
            {
                open = true;
            }
            else
            {
                open = false;
                message("Could not automatically read the Plot3D Q-file named `" + qfilename + "', please enter your information...");
            }
        }
        else
        {
            open = false;
            message("Could not automatically read the Plot3D Q-file, please enter your information...");
        }
    }
    
    // get file details if not already gotten above
    while ( !open )
    {
        ft = 10;
        fsize = Single;
        multidomain = true;
        std::string stuff = "";
        qfilename = "";
        
        while ( qfilename == std::string("") ||
               qfilename == std::string(" ") ||
               qfilename == std::string("\n") )
        {
            qfilename = question("what is the name of your Plot3D Q-file (without the extension)", std::string(""));
        }
        vfilename = qfilename + ".VAR";
        qfilename += ".P";
        
        stuff = question("is your Plot3D Q-file Unformatted, Formatted or Binary", std::string("Unformatted"));
        if ( stuff == std::string("Formatted") ||
            stuff == std::string("formatted") ||
            stuff == std::string("FORMATTED") ||
            stuff == std::string("F") ||
            stuff == std::string("f") )
        {
            format = Formatted;
            qfilename += 'F';
        }
        else if ( stuff == std::string("Binary") ||
                 stuff == std::string("binary") ||
                 stuff == std::string("BINARY") ||
                 stuff == std::string("B") ||
                 stuff == std::string("b") )
        {
            format = Binary;
            qfilename += 'B';
        }
        else
        {
            format = Unformatted;
            qfilename += 'U';
        }
        
        stuff = "";
        
        stuff = question("is your Plot3D Q-file single or double precision", std::string("Single"));
        if ( stuff == std::string("double") ||
            stuff == std::string("Double") ||
            stuff == std::string("DOUBLE") ||
            stuff == std::string("d") ||
            stuff == std::string("D") )
        {
            fsize = Double;
        }
        
        qfilename += "S";
        if ( fsize == Double )
        {
            qfilename += "D";
        }
        switch(format)
        {
            case Unformatted: case Binary:
                qfile.open(qfilename.c_str(),std::ios_base::binary);
                break;
            case Formatted:
                qfile.open(qfilename.c_str());
                break;
        }
        if ( !qfile )
        {
            ft -= 5;
            qscalar = false;
            message("could not find scalar Q-file `" + qfilename + "'");
        }
        else
        {
            qscalar = true;
        }
        qfile.close();
        if ( fsize == Double )
        {
            qfilename.resize(qfilename.size()-1);
        }
        qfilename.resize(qfilename.size()-1);
        qfilename += "V";
        if ( fsize == Double )
        {
            qfilename += "D";
        }
        switch(format)
        {
            case Unformatted: case Binary:
                qfile.open(qfilename.c_str(),std::ios_base::binary);
                break;
            case Formatted:
                qfile.open(qfilename.c_str());
                break;
        }
        if ( !qfile )
        {
            ft -= 5;
            qvector = false;
            debug("Plot3D::read_interface","could not find vector Q-file '" + qfilename + "'");
        }
        else
        {
            qvector = true;
        }
        qfile.close();
        vfile.open(vfilename.c_str());
        if ( vfile && ft != 0 )
        {
            numvar = 0;
            std::getline(vfile,vars[numvar]);
            while ( !vfile.eof() )
            {
                numvar++;
                if ( vars[numvar-1][0] == ' ' )
                {
                    for (i=0; i<vars[numvar-1].size()-1; i++)
                    {
                        vars[numvar-1][i] = vars[numvar-1][i+1];
                    }
                    vars[numvar-1].resize(vars[numvar-1].size()-1);
                }
                if ( vars[numvar-1] == std::string("rho") ||
                    vars[numvar-1] == std::string("density") ||
                    vars[numvar-1] == std::string("Rho") ||
                    vars[numvar-1] == std::string("Density") ||
                    vars[numvar-1] == std::string("RHO") ||
                    vars[numvar-1] == std::string("DENSITY") )
                {
                    variable = numvar;
                }
                if ( numvar + 1 > vars.size() )
                {
                    vars.resize(vars.size()+1);
                }
                std::getline(vfile, vars[numvar]);
            }
        }
        else if ( ft > 0 )
        {
            numvar = ft;
        }
        vfile.close();
        if ( numvar != vars.size() )
        {
            vars.resize(numvar);
        }
        if ( ft > 0 )
        {
            std::string number_word;
            number_word = question("Please select a variable to extract", vars, vars[variable - 1]);
            i = variable;
            variable = 0;
            if ( number_word == vars[0] )
            {
                variable = 1;
            }
            else if ( number_word == vars[1] )
            {
                variable = 2;
            }
            else if (  number_word == vars[2] )
            {
                variable = 3;
            }
            else if ( number_word == vars[3] )
            {
                variable = 4;
            }
            else if ( number_word == vars[4] )
            {
                variable = 5;
            }
            else if ( ft > 5 )
            {
                if ( number_word == vars[5] )
                {
                    variable = 6;
                }
                else if ( number_word == vars[6] )
                {
                    variable = 7;
                }
                else if ( number_word == vars[7] )
                {
                    variable = 8;
                }
                else if ( number_word == vars[8] )
                {
                    variable = 9;
                }
                else if ( number_word == vars[9] )
                {
                    variable = 10;
                }
            }
            if ( variable == 0 )
            {
                variable = i;
            }
        }
        
        stuff = "";
        
        multidomain = yesno("is your Plot3D Q-file multidomain", true);
        
        if ( variable <= 5 && qscalar )
        {
            if ( fsize == Double )
            {
                qfilename.resize(qfilename.size()-1);
            }
            qfilename.resize(qfilename.size()-1);
            qfilename += "S";
            if ( fsize == Double )
            {
                qfilename += "D";
            }
        }
        
        message("will read variable " + std::to_string(variable) + " (" + vars[variable - 1] + ") from the Q-file `" + qfilename + "'");
        
        if ( variable > 5 )
        {
            variable -= 5;
        }
        if ( Plot3D::extract_data(data, qfilename, format, fsize, multidomain, variable) )
        {
            open = true;
        }
        else
        {
            open = false;
            message("Could not open the Plot3D Q-file named `" + qfilename + stuff + "', please re-enter your information...");
        }
    }
    
    dataname = vars[variable-1];
    vars.clear();
}





/* ---------- read_data ---------- */
/// reads Plot3D data (or Q-file)
template<class T> bool Plot3D::read_data(std::valarray<T>& data1,
                                         std::valarray<T>& data2,
                                         std::valarray<T>& data3,
                                         std::valarray<T>& data4,
                                         std::valarray<T>& data5,
                                         const std::string& qfilename,
                                         const dataformat& format,
                                         const dataprecision& fsize,
                                         const bool& multidomain)
{
    std::ifstream qfile;
    
    switch(format)
    {
        case Unformatted: case Binary:
            qfile.open(qfilename.c_str(),std::ios_base::binary);
            break;
        case Formatted:
            qfile.open(qfilename.c_str());
            break;
    }
    if ( !qfile )
    {
        debug("Plot3D::read_data", "can't open Plot3D Q-file '" + qfilename + "'");
        throw;
        return false;
    }
    
    std::vector< std::valarray<T>* > data(5);
    data[0] = &data1;
    data[1] = &data2;
    data[2] = &data3;
    data[3] = &data4;
    data[4] = &data5;
    
    // byte swapping facilities
    bool byte_swapping = (std::endian::native == std::endian::little);
    
    // get information about number of domains
    std::size_t domains = 1;
    std::valarray<std::size_t> Nx;
    std::valarray<std::size_t> Ny;
    std::valarray<std::size_t> Nz;
    if ( multidomain )
    {
        Fortran::iread_line(qfile, domains, byte_swapping, format);
    }
    if ( Nx.size() == 0 && Ny.size() == 0 && Nz.size() == 0 )
    {
        Nx.resize(domains);
        Ny.resize(domains);
        Nz.resize(domains);
    }
    
    // get information about number of nodes
    std::size_t isize, itemp;
    Fortran::read_byte_info(qfile, isize, byte_swapping, format);
    if ( format == Unformatted && isize != domains * 3 * 4 )
    {
        debug("Plot3D::read_data", "error reading Plot3D Q-file `" + qfilename + "'");
        throw;
        return false;
    }
    std::size_t numnodes = 0;
    std::size_t loop;
    for (loop=0; loop<domains; loop++)
    {
        Fortran::iread_item(qfile, Nx[loop], byte_swapping, format);
        Fortran::iread_item(qfile, Ny[loop], byte_swapping, format);
        Fortran::iread_item(qfile, Nz[loop], byte_swapping, format);
        if ( Nx[loop] < 1 )
        {
            Nx[loop] = 1;
        }
        if ( Ny[loop] < 1 )
        {
            Ny[loop] = 1;
        }
        if ( Nz[loop] < 1 )
        {
            Nx[loop] = 1;
        }
        numnodes += Nx[loop] * Ny[loop] * Nz[loop];
    }
    Fortran::read_byte_info(qfile, itemp, byte_swapping, format);
    if ( itemp != isize )
    {
        debug("Plot3D::read_data", "inconsistent Plot3D Q-file `" + qfilename + "'");
        throw;
        return false;
    }
    
    // allocate memory
    for (std::size_t i=0; i<5; i++)
    {
        (data[i])->resize(numnodes);
    }
    
    // read actual data values
    std::size_t i, j, k;
    std::size_t node = 0, zero = 0;
    for (loop=0; loop<domains; loop++)
    {
        // read unused real variables in header
        Fortran::read_byte_info(qfile, isize, byte_swapping, format);
        if ( format == Unformatted && isize != 16 )
        {
            for (std::size_t i=0; i<5; i++)
            {
                (data[i])->resize(0);
            }
            debug("Plot3D::read_data", "inconsistent Plot3D Q-file `" + qfilename + "'");
            throw;
            return false;
        }
        Fortran::fread_blank(qfile, format, fsize);
        Fortran::fread_blank(qfile, format, fsize);
        Fortran::fread_blank(qfile, format, fsize);
        Fortran::fread_blank(qfile, format, fsize);
        Fortran::read_byte_info(qfile, itemp, byte_swapping, format);
        if ( itemp != isize )
        {
            for (std::size_t i=0; i<5; i++)
            {
                (data[i])->resize(0);
            }
            debug("Plot3D::read_data", "inconsistent Plot3D Q-file `" + qfilename + "'");
            throw;
            return false;
        }
        // get unformated header info
        Fortran::read_byte_info(qfile, isize, byte_swapping, format);
        // get data
        for (std::size_t qs=0; qs<5; qs++)
        {
            for (k=0; k<Nz[loop]; k++)
            {
                for (j=0; j<Ny[loop]; j++)
                {
                    for (i=0; i<Nx[loop]; i++)
                    {
                        node = zero + k * Ny[loop] * Nx[loop] + j * Nx[loop] + i;
                        Fortran::fread_item(qfile, (*(data[qs]))[node], byte_swapping, format, fsize);
                    }
                }
            }
        }
        // get unformated tail info
        Fortran::read_byte_info(qfile, itemp, byte_swapping, format);
        if ( itemp != isize )
        {
            for (std::size_t i=0; i<5; i++)
            {
                (data[i])->resize(0);
            }
            debug("Plot3D::read_data", "inconsistent Plot3D Q-file `" + qfilename + "'");
            throw;
            return false;
        }
        zero += Nx[loop] * Ny[loop] * Nz[loop];
    }
    
    qfile.close();
    
    return true;
}





/* ---------- extract_data ---------- */
/// reads Plot3D data (or Q-file)
template<class T> bool Plot3D::extract_data(std::valarray<T>& data,
                                            const std::string& qfilename,
                                            const dataformat& format,
                                            const dataprecision& fsize,
                                            const bool& multidomain,
                                            const unsigned short& variable)
{
    std::vector< std::valarray<T> > tdata(5);
    bool edtf = Plot3D::read_data(tdata[0], tdata[1], tdata[2], tdata[3], tdata[4], qfilename, format, fsize, multidomain);
    
    if ( edtf )
    {
        data.resize(tdata[variable-1].size());
        data = tdata[variable-1];
    }
    return edtf;
}





/* ---------- write_interface ---------- */
template<class T> bool Plot3D::write_interface(const std::vector< std::valarray<T> >& nodes,
                                               const std::valarray<std::size_t>& Nx,
                                               const std::valarray<std::size_t>& Ny,
                                               const std::valarray<std::size_t>& Nz,
                                               const std::string& gfilename)
{
    if ( Nx.size() <= 0 || Ny.size() <= 0 || Nz.size() <= 0 )
    {
        debug("Plot3D::write_interface", "can't write this grid");
        throw;
        return false;
    }
    
    // sort out dimensions
    std::size_t numnodes = 0;
    for (std::size_t loop=0; loop<Nx.size(); loop++)
    {
        numnodes += Nx[loop] * Ny[loop] * Nz[loop];
    }
    if ( nodes.size() != numnodes ||
        ( ( nodes[0].size() != 1 ) &&
         ( nodes[0].size() != 2 ) &&
         ( nodes[0].size() != 3 ) ) )
    {
        debug("Plot3D::write_interface", "can't write this grid");
        throw;
        return false;
    }
    
    std::string filename;
    bool open = false;
    
    // get file details
    dataformat format = Unformatted;
    
    while ( !open )
    {
        std::string stuff = "";
        while ( filename == std::string("") ||
               filename == std::string(" ") ||
               filename == std::string("\n") )
        {
            filename = question("what do you want to call your Plot3D file (without the extension)", gfilename);
        }
        filename += ".P";
        stuff = question("do you wish to save it as Unformatted, Formatted or Binary", std::string("Unformatted"));
        if ( stuff == std::string("Formatted") ||
            stuff == std::string("formatted") ||
            stuff == std::string("FORMATTED") ||
            stuff == std::string("F") ||
            stuff == std::string("f") )
        {
            format = Formatted;
            filename += "F";
        }
        else if ( stuff == std::string("Binary") ||
                 stuff == std::string("binary") ||
                 stuff == std::string("BINARY") ||
                 stuff == std::string("B") ||
                 stuff == std::string("b") )
        {
            format = Binary;
            filename += "B";
        }
        else if ( stuff == std::string("EXIT") ||
                 stuff == std::string("Exit") ||
                 stuff == std::string("exit") ||
                 stuff == std::string("QUIT") ||
                 stuff == std::string("Quit") ||
                 stuff == std::string("quit") )
        {
            return false;
        }
        else
        {
            format = Unformatted;
            filename += "U";
        }
        filename += "G";
        dataprecision fsize = Single;
        if ( sizeof(T) > 4 )
        {
            filename += 'D';
            fsize = Double;
        }
        
        std::valarray<T> X(nodes.size()), Y(nodes.size()), Z(nodes.size());
        for (std::size_t i=0; i<nodes.size(); i++)
        {
            X[i] = nodes[i][0];
            Y[i] = nodes[i][1];
            Z[i] = nodes[i][2];
        }
        std::valarray<bool> B(true,nodes.size());
        bool multidomain_ = false;
        bool blanking_ = false;
        if ( Plot3D::write(X, Y, Z, B, Nx, Ny, Nz, filename, format, fsize, multidomain_, blanking_) )
        {
            open = true;
        }
        else
        {
            debug("Plot3D::write_interface", "Could not open the file named `" + filename + "', please re-enter your information...");
            open = false;
        }
    }
    return true;
}





/* ---------- write ---------- */
template<class T> bool Plot3D::write(const std::valarray<T>& X,
                                     const std::valarray<T>& Y,
                                     const std::valarray<T>& Z,
                                     const std::valarray<bool>& B,
                                     const std::size_t& Nx,
                                     const std::size_t& Ny,
                                     const std::size_t& Nz,
                                     const std::string& filename,
                                     const dataformat& format)
{
    std::valarray<std::size_t> vNx(Nx, 1);
    std::valarray<std::size_t> vNy(Ny, 1);
    std::valarray<std::size_t> vNz(Nz, 1);
    return Plot3D::write(X, Y, Z, B, vNx, vNy, vNz, filename, format, Single, false, false);
}





/* ---------- write ---------- */
template<class T> bool Plot3D::write(const std::valarray<T>& X,
                                     const std::valarray<T>& Y,
                                     const std::valarray<T>& Z,
                                     const std::valarray<bool>& B,
                                     const std::valarray<std::size_t>& Nx,
                                     const std::valarray<std::size_t>& Ny,
                                     const std::valarray<std::size_t>& Nz,
                                     const std::string& filename,
                                     const dataformat& format,
                                     const dataprecision& fsize,
                                     const bool& multidomain,
                                     const bool& blanking)
{
    // sort out dimensions
    std::size_t numnodes = 0;
    std::size_t loop;
    std::size_t domains = Nx.size();
    if ( !multidomain && domains > 2 )
    {
        debug("Plot3D::write", "bad number of domains, try turning multidomain on");
        throw;
        return false;
    }
    for (loop=0; loop<domains; loop++)
    {
        numnodes += Nx[loop] * Ny[loop] * Nz[loop];
    }
    if ( X.size() != numnodes || Y.size() != numnodes || Z.size() != numnodes )
    {
        debug("Plot3D::write", "inconsistent Plot3D grid");
        throw;
        return false;
    }
    
    std::ofstream gridfile;
    switch(format)
    {
        case Unformatted: case Binary:
            gridfile.open(filename.c_str(),std::ios_base::binary);
            break;
        default:
            gridfile.open(filename.c_str());
            break;
    }
    if ( !gridfile )
    {
        debug("Plot3D::write", "can't open Plot3D grid file `" + filename + "'");
        throw;
        return false;
    }
    
    // byte swapping facilities
    bool byte_swapping = (std::endian::native == std::endian::little);
    
    // write information about number of domains
    if ( multidomain )
    {
        Fortran::iwrite_line(gridfile, domains, byte_swapping, format);
    }
    
    // write information about number of nodes
    std::size_t isize = 4 * 3 * domains;
    Fortran::write_byte_info(gridfile, isize, byte_swapping, format);
    for (loop=0; loop<domains; loop++)
    {
        Fortran::iwrite_item(gridfile, Nx[loop], byte_swapping, format);
        Fortran::iwrite_item(gridfile, Ny[loop], byte_swapping, format);
        Fortran::iwrite_item(gridfile, Nz[loop], byte_swapping, format);
    }
    if ( format == Formatted )
    {
        gridfile << std::endl;
    }
    
    Fortran::write_byte_info(gridfile, isize, byte_swapping, format);
    
    // write actual grid nodes
    std::size_t node = 0, izero = 0;;
    std::size_t counter = 0;
    std::size_t i, j, k;
    T zero = 0;
    switch(fsize)
    {
        case Double:
            isize = 8 * 3 * numnodes;
            break;
        default:
            isize = 4 * 3 * numnodes;
            break;
    }
    if ( blanking && B.size() == numnodes )
    {
        isize += numnodes * 4;
    }
    for (loop=0; loop<domains; loop++)
    {
        Fortran::write_byte_info(gridfile, isize, byte_swapping, format);
        // I values
        counter = 0;
        for (k=0; k<Nz[loop]; k++)
        {
            for (j=0; j<Ny[loop]; j++)
            {
                for (i=0; i<Nx[loop]; i++)
                {
                    Fortran::fwrite_item(gridfile, X[izero + k * Ny[loop] * Nx[loop] + j * Nx[loop] + i], byte_swapping, format, fsize);
                    if ( format == Formatted )
                    {
                        counter++;
                        if ( counter == 5 && node != numnodes - 1 )
                        {
                            gridfile << std::endl;
                            counter = 0;
                        }
                    }
                }
            }
        }
        if ( format == Formatted )
        {
            gridfile << std::endl;
        }
        // J values
        counter = 0;
        for (k=0; k<Nz[loop]; k++)
        {
            for (j=0; j<Ny[loop]; j++)
            {
                for (i=0; i<Nx[loop]; i++)
                {
                    Fortran::fwrite_item(gridfile, Y[izero + k * Ny[loop] * Nx[loop] + j * Nx[loop] + i], byte_swapping, format, fsize);
                    if ( format == Formatted )
                    {
                        counter++;
                        if ( counter == 5 && node != numnodes - 1 )
                        {
                            gridfile << std::endl;
                            counter = 0;
                        }
                    }
                }
            }
        }
        if ( format == Formatted )
        {
            gridfile << std::endl;
        }
        // K values
        counter = 0;
        for (k=0; k<Nz[loop]; k++)
        {
            for (j=0; j<Ny[loop]; j++)
            {
                for (i=0; i<Nx[loop]; i++)
                {
                    Fortran::fwrite_item(gridfile, Z[izero + k * Ny[loop] * Nx[loop] + j * Nx[loop] + i], byte_swapping, format, fsize);
                    if ( format == Formatted )
                    {
                        counter++;
                        if ( counter == 5 && node != numnodes - 1 )
                        {
                            gridfile << std::endl;
                            counter = 0;
                        }
                    }
                }
            }
        }
        // blanking
        if ( blanking && B.size() == numnodes )
        {
            if ( format == Formatted )
            {
                gridfile << std::endl;
            }
            counter = 0;
            for (k=0; k<Nz[loop]; k++)
            {
                for (j=0; j<Ny[loop]; j++)
                {
                    for (i=0; i<Nx[loop]; i++)
                    {
                        Fortran::bwrite_item(gridfile, B[izero + k * Ny[loop] * Nx[loop] + j * Nx[loop] + i], byte_swapping, format);
                        if ( format == Formatted )
                        {
                            counter++;
                            if ( counter == 38 )
                            {
                                gridfile << std::endl;
                                counter = 0;
                            }
                        }
                    }
                }
            }
        }
        Fortran::write_byte_info(gridfile, isize, byte_swapping, format);
        izero += Nx[loop] * Ny[loop] * Nz[loop];
    }
    
    if ( format == Formatted )
    {
        gridfile << std::endl;
    }
    
    gridfile.close();
    
    return true;
}





/* ---------- write_data_interface ---------- */
/// writes a Plot3D data, or Q, file
template<class T> bool Plot3D::write_data_interface(const std::valarray<T>& data,
                                                    const std::string& dataname,
                                                    const std::valarray<std::size_t>& Nx,
                                                    const std::valarray<std::size_t>& Ny,
                                                    const std::valarray<std::size_t>& Nz,
                                                    const std::string& gfilename)
{
    if ( Nx.size() <= 0 || Ny.size() <= 0 || Nz.size() <= 0 )
    {
        return false;
    }
    std::size_t numnodes = 0;
    for (std::size_t i=0; i<Nx.size(); i++)
    {
        numnodes += Nx[i] * Ny[i] * Nz[i];
    }
    if ( data.size() != numnodes )
    {
        return false;
    }
    
    std::string qfilename;
    std::ofstream qfile;
    bool open = false;
    dataformat format = Unformatted;
    
    // get file details
    while ( !open )
    {
        std::string stuff = "";
        qfilename = "";
        
        while ( qfilename == std::string("") ||
               qfilename == std::string(" ") ||
               qfilename == std::string("\n") )
        {
            qfilename = question("what do you wish to call your Plot3D Q-file (without the extension)", gfilename);
        }
        std::string ext = ".P";
        stuff = question("do you wish to save as Unformatted, Formatted or Binary", std::string("Unformatted"));
        if ( stuff == std::string("Formatted") ||
            stuff == std::string("formatted") ||
            stuff == std::string("FORMATTED") ||
            stuff == std::string("F") ||
            stuff == std::string("f") )
        {
            ext += "FS";
            format = Formatted;
        }
        else if ( stuff == std::string("Binary") ||
                 stuff == std::string("binary") ||
                 stuff == std::string("BINARY") ||
                 stuff == std::string("B") ||
                 stuff == std::string("b") )
        {
            ext += "BS";
            format = Binary;
        }
        else
        {
            ext += "US";
            format = Unformatted;
        }
        
        // which float to use
        dataprecision fsize = Single;
        if ( sizeof(T) > 4 )
        {
            fsize = Double;
            ext += "D";
        }
        std::valarray<T> nothing;
        bool multidomain_ = false;
        if ( Plot3D::write_data(data, nothing, nothing, nothing, nothing, Nx, Ny, Nz, qfilename+ext, format, fsize, multidomain_) )
        {
            Plot3D::write_var(qfilename+".VAR", dataname);
            open = true;
        }
        else
        {
            debug("Plot3D::write_data_interface", "Could not write the Plot3D Q-file named `" + qfilename + "', please re-enter your information...");
            open = false;
        }
    }
    
    return true;
}





/* ---------- write_data ---------- */
template<class T> bool Plot3D::write_data(const std::valarray<T>& data1,
                                          const std::valarray<T>& data2,
                                          const std::valarray<T>& data3,
                                          const std::valarray<T>& data4,
                                          const std::valarray<T>& data5,
                                          const std::size_t& Nx,
                                          const std::size_t& Ny,
                                          const std::size_t& Nz,
                                          const std::string& qfilename,
                                          const dataformat& format)
{
    std::valarray<std::size_t> vNx(Nx, 1);
    std::valarray<std::size_t> vNy(Ny, 1);
    std::valarray<std::size_t> vNz(Nz, 1);
    return Plot3D::write_data(data1, data2, data3, data4, data5, vNx, vNy, vNz, qfilename, format, Single, false);
}





/* ---------- write_data ---------- */
/// writes a Plot3D data (or Q-file)
template<class T> bool Plot3D::write_data(const std::valarray<T>& data1,
                                          const std::valarray<T>& data2,
                                          const std::valarray<T>& data3,
                                          const std::valarray<T>& data4,
                                          const std::valarray<T>& data5,
                                          const std::valarray<std::size_t>& Nx,
                                          const std::valarray<std::size_t>& Ny,
                                          const std::valarray<std::size_t>& Nz,
                                          const std::string& qfilename,
                                          const dataformat& format,
                                          const dataprecision& fsize,
                                          const bool& multidomain)
{
    if ( Nx.size() == 0 || Ny.size() == 0 || Nz.size() == 0 )
    {
        debug("Plot3D::write_data", "can't write this Plot3D grid data");
        throw;
        return false;
    }
    std::size_t domains = Nx.size();
    std::size_t numnodes = 0;
    for (std::size_t i=0; i<domains; i++)
    {
        numnodes += Nx[i] * Ny[i] * Nz[i];
    }
    
    std::ofstream qfile;
    switch(format)
    {
        case Unformatted: case Binary:
            qfile.open(qfilename.c_str(),std::ios_base::binary);
            break;
        default:
            qfile.open(qfilename.c_str());
            break;
    }
    if ( !qfile )
    {
        debug("Plot3D::write_data", "can't open Plot3D Q-file `" + qfilename + "'");
        throw;
        return false;
    }
    
    // byte swapping facilities
    bool byte_swapping = (std::endian::native == std::endian::little);
    
    // write information about number of domains
    if ( multidomain )
    {
        Fortran::iwrite_line(qfile, domains, byte_swapping, format);
    }
    
    // write information about number of nodes
    std::size_t isize = domains * 3 * 4;
    Fortran::write_byte_info(qfile, isize, byte_swapping, format);
    std::size_t loop;
    for (loop=0; loop<domains; loop++)
    {
        Fortran::iwrite_item(qfile, Nx[loop], byte_swapping, format);
        Fortran::iwrite_item(qfile, Ny[loop], byte_swapping, format);
        Fortran::iwrite_item(qfile, Nz[loop], byte_swapping, format);
        if ( format == Formatted )
        {
            if ( loop == domains - 1 )
            {
                qfile << std::endl;
            }
        }
    }
    Fortran::write_byte_info(qfile, isize, byte_swapping, format);
    
    // write actual data values
    std::size_t node = 0, zero = 0;
    std::size_t counter = 0;
    T rtemp;
    for (loop=0; loop<domains; loop++)
    {
        // write unused real variables in header
        switch(fsize)
        {
            case Double:
                isize = 8 * 4;
                break;
            default:
                isize = 4 * 4;
                break;
        }
        Fortran::write_byte_info(qfile, isize, byte_swapping, format);
        rtemp = T(1);
        Fortran::fwrite_item(qfile, rtemp, byte_swapping, format, fsize);
        rtemp = T(0);
        Fortran::fwrite_item(qfile, rtemp, byte_swapping, format, fsize);
        rtemp = T(1);
        Fortran::fwrite_item(qfile, rtemp, byte_swapping, format, fsize);
        rtemp = T(0);
        Fortran::fwrite_item(qfile, rtemp, byte_swapping, format, fsize);
        if ( format == Formatted )
        {
            qfile << std::endl;
        }
        Fortran::write_byte_info(qfile, isize, byte_swapping, format);
        // write unformated header info
        switch(fsize)
        {
            case Double:
                isize = numnodes * 8 * 5;
                break;
            default:
                isize = numnodes * 4 * 5;
                break;
        }
        Fortran::write_byte_info(qfile, isize, byte_swapping, format);
        
        // write data
        rtemp = T(0);
        std::valarray<T> nothing(T(0), numnodes);
        for (unsigned short qs=1; qs<=5; qs++)
        {
            counter = 0;
            for (std::size_t k=0; k<Nz[loop]; k++)
            {
                for (std::size_t j=0; j<Ny[loop]; j++)
                {
                    for (std::size_t i=0; i<Nx[loop]; i++)
                    {
                        node = zero + k * Ny[loop] * Nx[loop] + j * Nx[loop] + i;
                        const std::valarray<T>* pnt = &nothing;
                        switch(qs)
                        {
                            case 1:
                                if ( data1.size() == numnodes )
                                {
                                    pnt = &data1;
                                }
                                break;
                            case 2:
                                if ( data2.size() == numnodes )
                                {
                                    pnt = &data2;
                                }
                                break;
                            case 3:
                                if ( data3.size() == numnodes )
                                {
                                    pnt = &data3;
                                }
                                break;
                            case 4:
                                if ( data4.size() == numnodes )
                                {
                                    pnt = &data4;
                                }
                                break;
                            case 5:
                                if ( data5.size() == numnodes )
                                {
                                    pnt = &data5;
                                }
                                break;
                        }
                        Fortran::fwrite_item(qfile, (*pnt)[node], byte_swapping, format, fsize);
                        if ( format == Formatted )
                        {
                            counter++;
                            if ( counter == 5 && node != numnodes - 1 )
                            {
                                qfile << std::endl;
                                counter = 0;
                            }
                        }
                    }
                }
            }
            if ( format == Formatted )
            {
                qfile << std::endl;
            }
        }
        // write unformatted tail info
        Fortran::write_byte_info(qfile, isize, byte_swapping, format);
        zero += Nx[loop] * Ny[loop] * Nz[loop];
    }
    qfile.close();
    
    return true;
}





/* ---------- write_var ---------- */
bool Plot3D::write_var(const std::string& vfilename,
                       const std::string& var1,
                       const std::string& var2,
                       const std::string& var3,
                       const std::string& var4,
                       const std::string& var5)
{
    std::ofstream vfile(vfilename.c_str());
    vfile << " "
    << var1 << "\n "
    << var2 << "\n "
    << var3 << "\n "
    << var4 << "\n "
    << var5 << std::endl;
    vfile.close();
    return true;
}





/* ---------- get_vertices ---------- */
/// finds the vertices for the given Plot3D cell
void Plot3D::get_vertices(const unsigned short& dim,
                          const std::valarray<std::size_t>& Nx,
                          const std::valarray<std::size_t>& Ny,
                          const std::valarray<std::size_t>& Nz,
                          std::valarray<std::size_t>& vertices,
                          const std::size_t& cell)
{
    std::size_t point_zero = 0;
    // find which domain the cell lies in
    std::size_t domain = Nx.size();
    std::size_t zero = 0;
    std::size_t nnx, nny, nnz;
    for (std::size_t i=0; i<Nx.size(); i++)
    {
        nnx = ( Nx[i] < 3 ) ? 1 : Nx[i] - 1;
        nny = ( Ny[i] < 3 ) ? 1 : Ny[i] - 1;
        nnz = ( Nz[i] < 3 ) ? 1 : Nz[i] - 1;
        if ( zero <= cell && cell < zero + nnx * nny * nnz )
        {
            domain = i;
            break;
        }
        zero += nnx * nny * nnz;
        point_zero += Nx[i] * Ny[i] * Nz[i];
    }
    if ( domain >= Nx.size() )
    {
        debug("Plot3D::get_vertices", "counting error");
        vertices.resize(0);
        throw;
        return;
    }
    
    // find which K layer our cell is in
    std::size_t which = Nz[domain];
    for (std::size_t i=0; i<nnz; i++)
    {
        if ( zero <= cell && cell < zero + nnx * nny )
        {
            which = i;
            break;
        }
        zero += nnx * nny;
        point_zero += Nx[domain] * Ny[domain];
    }
    if ( which >= Nz[domain] )
    {
        debug("Plot3D::get_vertices", "counting error");
        vertices.resize(0);
        throw;
        return;
    }
    
    // find which J row our cell is in
    which = Ny[domain];
    for (std::size_t i=0; i<nny; i++)
    {
        if ( zero <= cell && cell < zero + nnx )
        {
            which = i;
            break;
        }
        zero += nnx;
        point_zero += Nx[domain];
    }
    if ( which >= Ny[domain] )
    {
        debug("Plot3D::get_vertices", "counting error");
        vertices.resize(0);
        throw;
        return;
    }
    
    // find which I column our cell is in
    which = cell - zero;
    //if ( which < 0 )
    //{
    //  debug("Plot3D::get_vertices", "counting error");
    //  vertices.resize(0);
    //  throw;
    //  return;
    //}
    point_zero += which;
    
    // allocate memory for cell vertices
    switch(dim)
    {
        case 1:
            vertices.resize(2);
            break;
        case 2:
            vertices.resize(4);
            break;
        case 3:
            vertices.resize(8);
            break;
        default:
            break;
    }
    
    // get the vertices for the cell
    vertices[0] = point_zero;
    vertices[1] = point_zero + 1;
    if ( dim >= 2 )
    {
        vertices[2] = point_zero + Nx[domain];
        vertices[3] = point_zero + Nx[domain] + 1;
    }
    if ( dim >= 3 )
    {
        point_zero += Nx[domain] * Ny[domain];
        vertices[4] = point_zero;
        vertices[5] = point_zero + 1;
        vertices[6] = point_zero + Nx[domain];
        vertices[7] = point_zero + Nx[domain] + 1;
    }
}





/* ---------- get_surrounds ---------- */
/// finds the surrounds cells for a given node
void Plot3D::get_surrounds(const unsigned short& dim,
                           const std::valarray<std::size_t>& Nx,
                           const std::valarray<std::size_t>& Ny,
                           const std::valarray<std::size_t>& Nz,
                           std::valarray<std::size_t>& surrounds,
                           const std::size_t& node)
{
    surrounds.resize(0);
    std::vector<std::size_t> surr;
    
    std::size_t zero = 0;
    std::size_t zero_cell = 0;
    // find which domain the node lies in
    std::size_t domain = Nx.size();
    std::size_t nnx = 0, nny = 0, nnz = 0;
    for (std::size_t i=0; i<Nx.size(); i++)
    {
        nnx = ( Nx[i] < 3 ) ? 1 : Nx[i] - 1;
        nny = ( Ny[i] < 3 ) ? 1 : Ny[i] - 1;
        nnz = ( Nz[i] < 3 ) ? 1 : Nz[i] - 1;
        if ( zero <= node && node < zero + Nx[i] * Ny[i] * Nz[i] )
        {
            domain = i;
            break;
        }
        zero += Nx[i] * Ny[i] * Nz[i];
        zero_cell += nnx * nny * nnz;
    }
    if ( domain >= Nx.size() )
    {
        debug("Plot3D::get_surrounds", "counting error");
        throw;
        return;
    }
    
    // find which K layer our node is in
    std::size_t Kl = Nz[domain];
    for (std::size_t i=0; i<Nz[domain]; i++)
    {
        if ( zero <= node && node < zero + Nx[domain] * Ny[domain] )
        {
            Kl = i;
            break;
        }
        zero += Nx[domain] * Ny[domain];
        if ( i < Nz[domain] - 2 )
        {
            zero_cell += nnx * nny;
        }
    }
    if ( Kl >= Nz[domain] )
    {
        debug("Plot3D::get_surrounds", "counting error");
        throw;
        return;
    }
    
    // find which J row our node is in
    std::size_t Jr = Ny[domain];
    for (std::size_t i=0; i<Ny[domain]; i++)
    {
        if ( zero <= node && node < zero + Nx[domain] )
        {
            Jr = i;
            break;
        }
        zero += Nx[domain];
        if ( i < Ny[domain] - 2 )
        {
            zero_cell += nnx;
        }
    }
    if ( Jr >= Ny[domain] )
    {
        debug("Plot3D::get_surrounds", "counting error");
        throw;
        return;
    }
    
    // find which I column our node is in
    std::size_t Ic = node - zero;
    //if ( Ic < 0 )
    //{
    //  debug("Plot3D::get_surrounds", "counting error");
    //  throw;
    //  return;
    //}
    zero += Ic;
    zero_cell += ( Ic == Nx[domain] - 1 ) ? Ic - 1 : Ic;
    
    switch(dim)
    {
        case 1:
            if ( Ic == 0 || Ic == Nx[domain] - 1 )
            {
                surr.push_back(zero_cell);
            }
            else
            {
                surr.push_back(zero_cell);
                surr.push_back(zero_cell-1);
            }
            break;
        case 2:
            if ( ( Ic == 0 || Ic == Nx[domain] - 1 ) &&
                ( Jr == 0 || Jr == Ny[domain] - 1 ) )
            {
                surr.push_back(zero_cell);
            }
            else if ( Ic == 0 || Ic == Nx[domain] - 1 )
            {
                surr.push_back(zero_cell);
                surr.push_back(zero_cell-nnx);
            }
            else if ( Jr == 0 || Jr == Ny[domain] - 1 )
            {
                surr.push_back(zero_cell);
                surr.push_back(zero_cell-1);
            }
            else
            {
                surr.push_back(zero_cell);
                surr.push_back(zero_cell-1);
                surr.push_back(zero_cell-nnx);
                surr.push_back(zero_cell-nnx-1);
            }
            break;
        case 3:
            if ( ( Ic == 0 || Ic == Nx[domain] - 1 ) &&
                ( Jr == 0 || Jr == Ny[domain] - 1 ) &&
                ( Kl == 0 || Kl == Nz[domain] - 1 ) )
            {
                surr.push_back(zero_cell);
            }
            else if ( ( Ic == 0 || Ic == Nx[domain] - 1 ) &&
                     ( Jr == 0 || Jr == Ny[domain] - 1 ) )
            {
                surr.push_back(zero_cell);
                surr.push_back(zero_cell-nnx*nny);
            }
            else if ( ( Jr == 0 || Jr == Ny[domain] - 1 ) &&
                     ( Kl == 0 || Kl == Nz[domain] - 1 ) )
            {
                surr.push_back(zero_cell);
                surr.push_back(zero_cell-1);
            }
            else if ( ( Kl == 0 || Kl == Nz[domain] - 1 ) &&
                     ( Ic == 0 || Ic == Nx[domain] - 1 ) )
            {
                surr.push_back(zero_cell);
                surr.push_back(zero_cell-nnx);
            }
            else if ( Ic == 0 || Ic == Nx[domain] - 1 )
            {
                surr.push_back(zero_cell);
                surr.push_back(zero_cell-nnx);
                surr.push_back(zero_cell-nnx*nny);
                surr.push_back(zero_cell-nnx*nny-nnx);
            }
            else if ( Jr == 0 || Jr == Ny[domain] - 1 )
            {
                surr.push_back(zero_cell);
                surr.push_back(zero_cell-1);
                surr.push_back(zero_cell-nnx*nny);
                surr.push_back(zero_cell-nnx*nny-1);
            }
            else if ( Kl == 0 || Kl == Nz[domain] - 1 )
            {
                surr.push_back(zero_cell);
                surr.push_back(zero_cell-1);
                surr.push_back(zero_cell-nnx);
                surr.push_back(zero_cell-nnx-1);
            }
            else
            {
                surr.push_back(zero_cell);
                surr.push_back(zero_cell-1);
                surr.push_back(zero_cell-nnx);
                surr.push_back(zero_cell-nnx-1);
                surr.push_back(zero_cell-nnx*nny);
                surr.push_back(zero_cell-nnx*nny-1);
                surr.push_back(zero_cell-nnx*nny-nnx);
                surr.push_back(zero_cell-nnx*nny-nnx-1);
            }
            break;
    }
    surrounds.resize(surr.size());
    std::copy(&surr[0], &surr[surr.size()], &surrounds[0]);
}





/* ---------- n2c ---------- */
template<class T> void Plot3D::n2c(const std::vector< std::valarray<T> >& nodes,
                                   const std::valarray<std::size_t>& Nx,
                                   const std::valarray<std::size_t>& Ny,
                                   const std::valarray<std::size_t>& Nz,
                                   const std::valarray<T>& ndata,
                                   std::valarray<T>& cdata)
{
    unsigned short dim = nodes[0].size();
    
    if ( ndata.size() != nodes.size() )
    {
        debug("Plot3D::n2c", "data doesn't match grid");
        cdata.resize(0);
        throw;
        return;
    }
    std::size_t cellsz = 0;
    for (std::size_t i=0; i<Nx.size(); i++)
    {
        cellsz += ( ( Nx[i] < 3 ) ? 1 : Nx[i] - 1 )
        * ( ( Ny[i] < 3 ) ? 1 : Ny[i] - 1 )
        * ( ( Nz[i] < 3 ) ? 1 : Nz[i] - 1 );
    }
    cdata.resize(cellsz, T(0));
    
    for (std::size_t i=0; i<cdata.size(); i++)
    {
        std::valarray<std::size_t> vertices;
        Plot3D::get_vertices(dim, Nx, Ny, Nz, vertices, i);
        T total_length = T(0);
        std::valarray<T> cntr(T(0),nodes[0].size());
        for (std::size_t vertex=0; vertex<vertices.size(); vertex++)
        {
            cntr += nodes[vertices[vertex]];
        }
        cntr /= T(vertices.size());
        std::valarray<T> diff(T(0),nodes[0].size());
        for (std::size_t vertex=0; vertex<vertices.size(); vertex++)
        {
            diff = nodes[vertices[vertex]] - cntr;
            total_length += norm(diff);
        }
        for (std::size_t vertex=0; vertex<vertices.size(); vertex++)
        {
            diff = nodes[vertices[vertex]] - cntr;
            T length = norm(diff);
            cdata[i] += ( length / total_length ) * ndata[vertices[vertex]];
        }
    }
}





/* ---------- c2n ---------- */
template<class T> void Plot3D::c2n(const std::vector< std::valarray<T> >& nodes,
                                   const std::valarray<std::size_t>& Nx,
                                   const std::valarray<std::size_t>& Ny,
                                   const std::valarray<std::size_t>& Nz,
                                   const std::valarray<T>& cdata,
                                   std::valarray<T>& ndata)
{
    unsigned short dim = nodes[0].size();
    std::size_t cellsz = 0;
    for (std::size_t i=0; i<Nx.size(); i++)
    {
        cellsz += ( ( Nx[i] < 3 ) ? 1 : Nx[i] - 1 )
        * ( ( Ny[i] < 3 ) ? 1 : Ny[i] - 1 )
        * ( ( Nz[i] < 3 ) ? 1 : Nz[i] - 1 );
    }
    ndata.resize(nodes.size());
    if ( cdata.size() != cellsz )
    {
        debug("Plot3D::c2n", "data doesn't match grid");
        ndata.resize(0);
        throw;
        return;
    }
    
    for (std::size_t i=0; i<ndata.size(); i++)
    {
        std::valarray<std::size_t> surr;
        Plot3D::get_surrounds(dim, Nx, Ny, Nz, surr, i);
        T total_length = T(0);
        std::valarray<T> lengths(T(0),surr.size());
        for (std::size_t cell=0; cell<surr.size(); cell++)
        {
            std::valarray<T> cntr(T(0),nodes[0].size());
            std::valarray<std::size_t> vertices;
            Plot3D::get_vertices(dim, Nx, Ny, Nz, vertices, surr[cell]);
            for (std::size_t vert=0; vert<vertices.size(); vert++)
            {
                cntr += nodes[vertices[vert]];
            }
            cntr /= T(vertices.size());
            std::valarray<T> tmp = cntr - nodes[i];
            lengths[cell] = norm(tmp);
            total_length += lengths[cell];
        }
        for (std::size_t cell=0; cell<surr.size(); cell++)
        {
            ndata[i] += ( lengths[cell] / total_length ) * cdata[surr[cell]];
        }
    }
}





/* ---------- write_neighbours ---------- */
template<class T> bool Plot3D::write_neighbours(const std::vector< std::valarray<T> >& nodes,
                                                const std::valarray<std::size_t>& Nx,
                                                const std::valarray<std::size_t>& Ny,
                                                const std::valarray<std::size_t>& Nz,
                                                const std::vector< std::valarray<std::size_t> >& neighbours,
                                                const std::string& filename,
                                                const dataformat& format)
{
    if ( neighbours.size() == 0 )
    {
        debug("Plot3D::write_neighbours", "no neighbours to write");
        throw;
        return false;
    }
    
    // byte swapping facilities
    bool byte_swapping = (std::endian::native == std::endian::little);
    
    std::ofstream nfile;
    switch(format)
    {
        case Unformatted: case Binary:
            nfile.open(filename.c_str(),std::ios_base::binary);
            break;
        default:
            nfile.open(filename.c_str());
            break;
    }
    
    // write number of domains
    Fortran::iwrite_line(nfile, Nx.size(), byte_swapping, format);
    
    // write domain sizes
    std::size_t isize = Nx.size() * 4 * 3;
    Fortran::write_byte_info(nfile, isize, byte_swapping, format);
    std::size_t loop;
    for (loop=0; loop<Nx.size(); loop++)
    {
        Fortran::iwrite_item(nfile, Nx[loop], byte_swapping, format);
        Fortran::iwrite_item(nfile, Ny[loop], byte_swapping, format);
        Fortran::iwrite_item(nfile, Nz[loop], byte_swapping, format);
        if ( format == Formatted )
        {
            if ( loop == Nx.size() - 1 )
            {
                nfile << std::endl;
            }
        }
    }
    Fortran::write_byte_info(nfile, isize, byte_swapping, format);
    
    // write neighbours
    for (loop=0; loop<neighbours.size(); loop++)
    {
        Fortran::iwrite_vector(nfile, neighbours[loop], byte_swapping, format);
    }
    
    nfile.close();
    
    return true;
}





/* ---------- to_PNG ---------- */
template<class T> bool Plot3D::to_PNG(const std::valarray<T>& X,
                                      const std::valarray<T>& Y,
                                      const std::valarray<T>& Z,
                                      const std::valarray<std::size_t>& Nx,
                                      const std::valarray<std::size_t>& Ny,
                                      const std::valarray<std::size_t>& Nz,
                                      const std::valarray<T>& data,
                                      const std::string& pngfilename,
                                      const Plot3D::values& get_switch)
{
    if ( Nx.size() > 1 || Ny.size() > 1 || Nz.size() > 1 )
    {
        debug("Plot3D::to_PNG", "can't convert multidomain Plot3D grids");
        throw;
        return false;
    }
    
    // determine structure of the grid and set up PNG structure
    png_uint_32 Nrows = 1;
    png_uint_32 Ncols = 1;
    std::valarray<T> ddata;
    switch(get_switch)
    {
        case Plot3D::NODES:
            ddata.resize(data.size());
            ddata = data;
            if ( Nx[0] == 1 && Ny[0] == 1 && Nz[0] == 1 )
            { // ZERO
            }
            else if ( Nx[0] > 1 && Ny[0] == 1 && Nz[0] == 1 )
            { // ONE_x
                Ncols = Nx[0];
            }
            else if ( Nx[0] == 1 && Ny[0] > 1 && Nz[0] == 1 )
            { // ONE_y
                Ncols = Ny[0];
            }
            else if ( Nx[0] == 1 && Ny[0] == 1 && Nz[0] > 1 )
            { // ONE_z
                Ncols = Nz[0];
            }
            else if ( Nx[0] > 1 && Ny[0] > 1 && Nz[0] == 1 )
            { // TWO_xy
                Ncols = Nx[0];
                Nrows = Ny[0];
            }
            else if ( Nx[0] == 1 && Ny[0] > 1 && Nz[0] > 1 )
            { // TWO_yz
                Ncols = Ny[0];
                Nrows = Nz[0];
            }
            else if ( Nx[0] > 1 && Ny[0] == 1 && Nz[0] > 1 )
            { // TWO_zx
                Ncols = Nz[0];
                Nrows = Nx[0];
                for (std::size_t j=0; j<Nx[0]; j++)
                {
                    for (std::size_t i=0; i<Nz[0]; i++)
                    {
                        ddata[j*Nz[0]+i] = data[i*Nx[0]+j];
                    }
                }
            }
            else
            {
                debug("Plot3D::to_PNG", "can't determine dimensions");
                throw;
                return false;
            }
            break;
        case Plot3D::CENTROIDS:
            if ( Nx[0] == 1 && Ny[0] == 1 && Nz[0] == 1 )
            { // ZERO
                ddata.resize(1);
                ddata[0] = data[0];
            }
            else if ( ( Nx[0] > 1 && Ny[0] == 1 && Nz[0] == 1 ) ||
                     ( Nx[0] == 1 && Ny[0] > 1 && Nz[0] == 1 ) ||
                     ( Nx[0] == 1 && Ny[0] == 1 && Nz[0] > 1 ) )
            { // ONE
                if ( Nx[0] > 1 )
                { // ONE_x
                    Ncols = Nx[0] - 1;
                }
                if ( Ny[0] > 1 )
                { // ONE_y
                    Ncols = Ny[0] - 1;
                }
                if ( Nz[0] > 1 )
                { // ONE_z
                    Ncols = Nz[0] - 1;
                }
                ddata.resize(Ncols);
                for (std::size_t i=0; i<Ncols; i++)
                {
                    ddata[i] = 0.5 * (data[i] + data[i+1]);
                }
            }
            else if ( ( Nx[0] > 1 && Ny[0] > 1 && Nz[0] == 1 ) ||
                     ( Nx[0] == 1 && Ny[0] > 1 && Nz[0] > 1 ) ||
                     ( Nx[0] > 1 && Ny[0] == 1 && Nz[0] > 1 ) )
            { // TWO
                if ( Nz[0] == 1 )
                { // TWO_xy
                    Ncols = Nx[0] - 1;
                    Nrows = Ny[0] - 1;
                }
                else if ( Nx[0] == 1 )
                { // TWO_yz
                    Ncols = Ny[0] - 1;
                    Nrows = Nz[0] - 1;
                }
                else if ( Ny[0] == 1 )
                { // TWO_zx
                    Ncols = Nz[0] - 1;
                    Nrows = Nx[0] - 1;
                }
                ddata.resize(Ncols*Nrows);
                for (std::size_t j=0; j<Nrows; j++)
                {
                    for (std::size_t i=0; i<Ncols; i++)
                    {
                        std::size_t i1 = j*(Ncols+1)+i;
                        std::size_t i2 = j*(Ncols+1)+i+1;
                        std::size_t i3 = (j+1)*(Ncols+1)+i;
                        std::size_t i4 = (j+1)*(Ncols+1)+i+1;
                        std::valarray<T> p1(3);
                        std::valarray<T> p2(3);
                        std::valarray<T> p3(3);
                        std::valarray<T> p4(3);
                        p1[0] = X[i1];
                        p1[1] = Y[i1];
                        p1[2] = Z[i1];
                        p2[0] = X[i2];
                        p2[1] = Y[i2];
                        p2[2] = Z[i2];
                        p3[0] = X[i3];
                        p3[1] = Y[i3];
                        p3[2] = Z[i3];
                        p4[0] = X[i4];
                        p4[1] = Y[i4];
                        p4[2] = Z[i4];
                        std::valarray<T> pc = 0.25 * ( p1 + p2 + p3 + p4 );
                        T l1 = norm(&pc, &p1);
                        T l2 = norm(&pc, &p2);
                        T l3 = norm(&pc, &p3);
                        T l4 = norm(&pc, &p4);
                        T lc = l1 + l2 + l3 + l4;
                        ddata[j*Ncols+i] = (l1/lc)*data[i1] + (l2/lc)*data[i2] + (l3/lc)*data[i3] + (l4/lc)*data[i4];
                    }
                }
                if ( Ny[0] == 1 )
                { // TWO_zx (flip axis from xz to zx)
                    std::valarray<T> tdata(ddata.size());
                    for (std::size_t j=0; j<Nrows; j++)
                    {
                        for (std::size_t i=0; i<Ncols; i++)
                        {
                            tdata[j*Ncols+i] = ddata[i*Nrows+j];
                        }
                    }
                    ddata = tdata;
                }
            }
            else
            {
                debug("Plot3D::to_PNG", "unsupported dimension");
                throw;
                return false;
            }
            break;
    }
    
    // swap data from Cartesian format to row-col image format
    for (std::size_t i=0; i<std::size_t(std::floor(double(Nrows)/2)); i++)
    {
        for (std::size_t j=0; j<Ncols; j++)
        {
            std::swap(ddata[i*Ncols+j], ddata[(Nrows-1-i)*Ncols+j]);
        }
    }
    
    // write the PNG file
    png_FILE_p fp = std::fopen(pngfilename.c_str(), "wb");
    if ( !fp )
    {
        debug("Plot3D::to_PNG", "can't open PNG file `" + pngfilename + "'");
        std::fclose(fp);
        throw;
        return false;
    }
    png_structp png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING, nullptr, (png_error_ptr)nullptr, (png_error_ptr)nullptr);
    if ( !png_ptr )
    {
        debug("Plot3D::to_PNG", "can't create PNG pointer");
        std::fclose(fp);
        throw;
        return false;
    }
    png_infop info_ptr = png_create_info_struct(png_ptr);
    if ( !info_ptr )
    {
        debug("Plot3D::to_PNG", "can't create PNG info_pointer");
        std::fclose(fp);
        png_destroy_write_struct(&png_ptr, (png_infopp)nullptr);
        throw;
        return false;
    }
    if ( setjmp(png_jmpbuf(png_ptr)) )
    {
        debug("Plot3D::to_PNG", "can't set PNG buffer");
        std::fclose(fp);
        png_destroy_write_struct(&png_ptr, &info_ptr);
        throw;
        return false;
    }
    png_init_io(png_ptr, fp);
    
    // IHDR
    int bit_depth = 16; // 16bit grayscale
    png_set_IHDR(png_ptr, info_ptr, Ncols, Nrows, bit_depth, PNG_COLOR_TYPE_GRAY, PNG_INTERLACE_NONE, PNG_COMPRESSION_TYPE_DEFAULT, PNG_FILTER_TYPE_DEFAULT);
    
    // pCAL
    T dmin = ddata.min();
    T dmax = ddata.max();
    png_uint_32 no_colors = png_uint_32(ipow(2,bit_depth));
    png_int_32 X0 = 0;
    png_int_32 X1 = no_colors - 1;
    int nparams = 2; // linear_data = p0 + p1 * sample / ( X1 - X0 );
    png_charp units = (char*)"";
    T p0 = dmin;
    T p1 = (dmax - dmin)*((T(X1-X0))/(T(X1)));
    std::string p0s = ntoIEEE(p0);
    std::string p1s = ntoIEEE(p1);
    png_charp params[2] = { (char*)p0s.c_str(), (char*)p1s.c_str() };
    png_set_pCAL(png_ptr, info_ptr, (char*)"Plot3D", X0, X1, PNG_EQUATION_LINEAR, nparams, units, params);
    
    // tIME
    png_time mod_time;
    png_convert_from_time_t(&mod_time, std::time(0));
    png_set_tIME(png_ptr, info_ptr, &mod_time);
    
    // necessities for IEEEdata
    std::vector<std::string> IEEEs(0);
    std::vector<std::string> IEEEsd(0);
    std::string sdn = "";
    std::size_t nodier = 2000; // number of data in each row
    for (std::size_t i=0; i<ddata.size(); i+=nodier)
    {
        std::size_t dvs = i + nodier;
        dvs = ( ddata.size() < dvs ) ? ddata.size() - i : nodier;
        std::valarray<T> data_v = ddata[std::slice(i, dvs, 1)];
        IEEEs.push_back(vtoIEEE(data_v));
    }
    IEEEsd.resize(IEEEs.size());
    for (std::size_t i=0; i<IEEEsd.size(); i++)
    {
        IEEEsd[i] = "IEEEdata" + std::to_string(i);
    }
    sdn = std::to_string(IEEEs.size());
    
    // tEXT
    int num_text = 9 + IEEEs.size();
    png_text text_ptr[num_text];
    text_ptr[0].key = (char*)"Title";
    text_ptr[0].text = (char*)"Plot3D_to_PNG";
    text_ptr[0].text_length = 13;
    text_ptr[0].compression = PNG_TEXT_COMPRESSION_NONE;
    text_ptr[1].key = (char*)"Author";
    text_ptr[1].text = (char*)"Rado Faletic";
    text_ptr[1].text_length = 12;
    text_ptr[1].compression = PNG_TEXT_COMPRESSION_NONE;
    text_ptr[2].key = (char*)"Description";
    text_ptr[2].text = (char*)"slice of a Plot3D grid";
    text_ptr[2].text_length = 22;
    text_ptr[2].compression = PNG_TEXT_COMPRESSION_NONE;
    std::string png_cyear = "(c) 2003, ";
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
    text_ptr[5].text = (char*)"plot3dpng";
    text_ptr[5].text_length = 9;
    text_ptr[5].compression = PNG_TEXT_COMPRESSION_NONE;
    text_ptr[6].key = (char*)"Source";
    text_ptr[6].text = (char*)"Plot3D";
    text_ptr[6].text_length = 6;
    text_ptr[6].compression = PNG_TEXT_COMPRESSION_NONE;
    std::string swhere = ( get_switch == NODES ) ? "NODES" : "CENTROIDS";
    text_ptr[7].key = (char*)"where";
    text_ptr[7].text = (char*)swhere.c_str();
    text_ptr[7].text_length = swhere.size();
    text_ptr[7].compression = PNG_TEXT_COMPRESSION_NONE;
    text_ptr[8].key = (char*)"ndatarows"; // stores how many data arrays there are
    text_ptr[8].text = (char*)sdn.c_str();
    text_ptr[8].text_length = sdn.size();
    text_ptr[8].compression = PNG_TEXT_COMPRESSION_NONE;
    for (std::size_t i=0; i<IEEEs.size(); i++)
    {
        text_ptr[9+i].key = (char*)IEEEsd[i].c_str();
        text_ptr[9+i].text = (char*)IEEEs[i].c_str();
        text_ptr[9+i].text_length = IEEEs[i].size();
        text_ptr[9+i].compression = PNG_TEXT_COMPRESSION_zTXt;
        //text_ptr[9+i].compression = PNG_TEXT_COMPRESSION_NONE;
    }
    png_set_text(png_ptr, info_ptr, text_ptr, num_text);
    
    // convert from T to int
    std::valarray<png_uint_16> idata(ddata.size());
    for (std::size_t i=0; i<idata.size(); i++)
    {
        idata[i] = png_uint_16((ddata[i] - p0)*(T(X1-X0))/p1);
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





#endif /* _PLOT3D_ */
