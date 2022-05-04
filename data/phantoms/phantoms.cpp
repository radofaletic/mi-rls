/**
 phantom.cpp
 
 generate a Plot3D file of various phantom objects
 
 Rado Faletic
 Department of Physics
 Faculty of Science
 Australian National University
 ACT  0200  Australia
 email: Rado.Faletic@anu.edu.au
 
 20th July 2004
 4th May 2022
 */

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <numbers>
#include <sstream>
#include <string>
#include <valarray>
#include "angles.h"
#include "plot3d.h"

template<class T> bool is_ellipse(const T& ix, const T& iy,
                                  const T& t,
                                  const T& rx, const T& ry,
                                  const T& cx, const T& cy)
{
    T tox = ix - cx;
    T toy = iy - cy;
    T ox = tox;
    T oy = toy;
    if ( t == T(0) )
    {
        ox /= rx;
        oy /= ry;
    }
    else
    {
        ox = (tox/rx)*std::cos(Angle::pi*t/180) + (toy/rx)*std::sin(Angle::pi*t/180);
        oy = -(tox/ry)*std::sin(Angle::pi*t/180) + (toy/ry)*std::cos(Angle::pi*t/180);
    }
    return ( ox * ox + oy * oy <= T(1) );
}

template<class T> bool is_ellipsoid(const T& ix, const T& iy, const T& iz,
                                    const T& t,
                                    const T& rx, const T& ry, const T& rz,
                                    const T& cx, const T& cy, const T& cz)
{
    T tox = ix - cx;
    T toy = iy - cy;
    T toz = iz - cz;
    T ox = tox;
    T oy = toy;
    T oz = toz;
    if ( t == T(0) )
    {
        ox /= rx;
        oy /= ry;
        oz /= rz;
    }
    else
    {
        ox = (tox/rx)*std::cos(Angle::pi*t/180.0) + -(toz/rx)*std::sin(Angle::pi*t/180.0);
        oy = (toy/ry);
        oz = (tox/rz)*std::sin(Angle::pi*t/180.0) + (toz/rz)*std::cos(Angle::pi*t/180.0);
    }
    return ( ox * ox + oy * oy + oz * oz <= T(1) );
}

int main(int argc, char* argv[])
{
    unsigned int sizeofgrid = 100;
    unsigned int s3 = sizeofgrid;
    std::string mode = "scalar";
    std::string nameoffile = mode;
    unsigned short dim = 3;
    std::string sext = "_3d";
    float sscale = 1;
    
    if ( argc == 1 )
    {
        std::cout << "\nuse the --help option to learn more\n" << std::endl;
        throw;
        return 0;
    }
    for (unsigned int i=1; i<argc; i++)
    {
        std::ostringstream input(argv[i]);
        if ( input.str().substr(0,6) == std::string("--help") )
        {
            std::cout << "\nby Rado Faletic, 2003, 2004, 2022\n\n"
            << "below is a list of flags:\n\n"
            << "--size=x\n\tsize of the Plot3D grid (" << sizeofgrid << ")\n"
            <<"--mode=<phantom>\n\tthe shape to generate (" << mode << ")\n\tavailable phantoms are:\n\t\tsingle\n\t\tscalar\n\t\tdiamond\n\t\tcircle | sphere\n\t\tsquare | cube\n\t\tpyramid\n\t\tbell\n\t\tshepp-logan\n\t\tshepp-logan_bold\n\t\tstark\n\t\tdorn\n"
            << "--output=<outputfilename>\n\tthe Plot3D output file name (" << nameoffile << ")\n"
            <<"--2d\n\tgenerate a 2d, as opposed to 3d, data set"
            <<"--scale=<s>\n\tlength between nodes (" << sscale << ")\n"
            << std::endl;
            return 1;
        }
        else if ( input.str().substr(0,7) == std::string("--size=") && mode != "single" )
        {
            std::istringstream choice(input.str().substr(7,input.str().size()-7));
            choice >> sizeofgrid;
        }
        else if ( input.str().substr(0,9) == std::string("--output=") )
        {
            nameoffile = input.str().substr(9,input.str().size()-9);
        }
        else if ( input.str().substr(0,7) == std::string("--mode=") )
        {
            mode = input.str().substr(7,input.str().size()-7);
            if ( mode == "single" )
            {
                sizeofgrid = 2;
            }
            else if ( mode == "square" )
            {
                mode = "cube";
                dim = 2;
                sext = "_2d";
                s3 = 1;
            }
            else if ( mode == "circle" )
            {
                mode = "sphere";
                dim = 2;
                sext = "_2d";
                s3 = 1;
            }
            else if ( mode == "smooth" || mode == "blob" )
            {
                mode = "bell";
            }
            if ( mode != "single" &&
                mode != "scalar" &&
                mode != "diamond" &&
                mode != "sphere" &&
                mode != "cube" &&
                mode != "pyramid" &&
                mode != "bell" &&
                mode != "shepp-logan" &&
                mode != "shepp-logan_bold" &&
                mode != "stark" &&
                mode != "dorn" )
            {
                mode = "scalar";
            }
            if ( nameoffile == "scalar" )
            {
                nameoffile = mode;
            }
        }
        else if ( input.str().substr(0,4) == std::string("--2d") )
        {
            dim = 2;
            sext = "_2d";
            s3 = 1;
        }
        else if ( input.str().substr(0,4) == std::string("--3d") )
        {
            dim = 3;
            sext = "_3d";
        }
        else if ( input.str().substr(0,8) == std::string("--scale=") )
        {
            std::istringstream choice(input.str().substr(8,input.str().size()-8));
            choice >> sscale;
        }
        else
        {
            std::cout << "unrecognised option '" << input.str() << "'\nuse the --help option to learn more" << std::endl;
            return 0;
        }
        
    }
    if ( mode == "single" )
    {
        switch(dim)
        {
            case 2:
                sizeofgrid = 1;
                break;
            case 3:
                sizeofgrid = 2;
                break;
        }
    }
    if ( dim == 3 )
    {
        sizeofgrid++;
        s3 = sizeofgrid;
    }
    size_t counter = 0;
    std::valarray<float> X(float(0), sizeofgrid*sizeofgrid*s3);
    std::valarray<float> Y(float(0), sizeofgrid*sizeofgrid*s3);
    std::valarray<float> Z(float(0), sizeofgrid*sizeofgrid*s3);
    // write co-ordinates
    if ( dim == 3 )
    {
        std::cout << sizeofgrid << "x" << sizeofgrid << "x" << s3
        << " output to file '" << nameoffile << sext << "'" << std::endl;
        for (unsigned int k=0; k<s3; k++)
        {
            for (unsigned int j=0; j<sizeofgrid; j++)
            {
                for (unsigned int i=0; i<sizeofgrid; i++)
                {
                    X[counter] = i * sscale;
                    Y[counter] = j * sscale;
                    Z[counter] = k * sscale;
                    counter++;
                }
            }
        }
        std::valarray<bool> B(true, sizeofgrid*sizeofgrid*s3);
        Plot3D::write(X, Y, Z, B, sizeofgrid, sizeofgrid, s3, nameoffile+sext+".PBG", Binary);
    }
    
    // write values
    counter = 0;
    std::valarray<float> data(float(0), sizeofgrid*sizeofgrid*s3);
    
    // phantom data
    double scale = 0.01 * ( sizeofgrid - 1.0 );
    double centre = 0.5 * scale;
    double radius = centre * 2 / 3;
    double m = 4.0 / 3.0;
    double c1 = 0.9*scale - m * 0.5*scale;
    double c2 = 0.9*scale + m * 0.5*scale;
    double scalar = 0.0;
    bool sl_bold = ( mode == "shepp-logan_bold" ) ? true : false;
    bool skeep = false;
    for (unsigned int k=0; k<s3; k++)
    {
        double z = k * 0.01;
        for (unsigned int j=0; j<sizeofgrid; j++)
        {
            double y = j * 0.01;
            for (unsigned int i=0; i<sizeofgrid; i++)
            {
                double x = i * 0.01;
                if ( mode == "single" ) // a single grid cell
                {
                    data[counter] = 1;
                }
                else if ( mode == "diamond" ) // an opaque diamond
                {
                    skeep = false;
                    scalar = std::abs( x - centre ) + std::abs( y - centre );
                    if ( dim > 2 )
                    {
                        scalar += std::abs( z - centre );
                    }
                    if ( scalar <= radius )
                    {
                        data[counter] = 1;
                    }
                    else
                    {
                        data[counter] = 0;
                    }
                }
                else if ( mode == "sphere" ) // an opaque sphere
                {
                    skeep = false;
                    scalar = std::pow(x - centre, 2) + std::pow(y - centre, 2);
                    if ( dim > 2 )
                    {
                        scalar += std::pow(z - centre, 2);
                    }
                    scalar = std::sqrt(scalar);
                    if ( scalar <= radius )
                    {
                        data[counter] = 1;
                    }
                    else
                    {
                        data[counter] = 0;
                    }
                }
                else if ( mode == "cube" ) // an opaque cube
                {
                    skeep = false;
                    scalar = std::max(std::abs(x - centre), std::abs(y - centre));
                    if ( dim > 2 )
                    {
                        scalar = std::max(std::abs(scalar), std::abs(z - centre));
                    }
                    if ( scalar <= radius )
                    {
                        data[counter] = 1;
                    }
                    else
                    {
                        data[counter] = 0;
                    }
                }
                else if ( mode == "pyramid" ) // a pyramid
                {
                    skeep = false;
                    c1 = std::pow(centre, 2) + std::pow(centre, 2);
                    if ( dim > 2 )
                    {
                        c1 += std::pow(centre, 2);
                    }
                    c1 = std::sqrt(c1);
                    scalar = ( 1 - std::abs(x - centre) / c1 ) * ( 1 - std::abs(y - centre) / c1 );
                    if ( dim > 2 )
                    {
                        scalar *= ( 1 - std::abs(z - centre) / c1 );
                    }
                    data[counter] = scalar;
                }
                else if ( mode == "bell" )
                {
                    skeep = true;
                    c1 = std::pow(centre, 2) + std::pow(centre, 2);
                    scalar = std::pow(x - centre, 2) + std::pow(y - centre, 2);
                    if ( dim > 2 )
                    {
                        c1 += std::pow(centre, 2);
                        scalar += std::pow(z - centre, 2);
                    }
                    c1 = std::sqrt(c1);
                    scalar = std::sqrt(scalar);
                    scalar = std::cos( Angle::pi * scalar / ( 2.0 * c1 ) );
                    scalar *= scalar;
                }
                else if ( mode == "shepp-logan" || mode == "shepp-logan_bold" ) // L. A. Shepp and B. F. Logan, "Reconstructing Interior Head Tissue from X-Ray Transmissions"
                {
                    skeep = true;
                    scalar = 0.0;
                    if ( dim == 2 )
                    {
                        // a
                        if ( is_ellipse(x-centre, y-centre, 0.0, 0.69*scale/2, 0.92*scale/2, 0.0, 0.0) )
                        {
                            scalar += ( sl_bold ) ? 1.0 : 2.0;
                        }
                        // b
                        if ( is_ellipse(x-centre, y-centre, 0.0, 0.6624*scale/2, 0.874*scale/2, 0.0, -0.0184*scale/2) )
                        {
                            scalar += ( sl_bold ) ? -0.8 : -0.98;
                        }
                        // c
                        if ( is_ellipse(x-centre, y-centre, -18.0, 0.11*scale/2, 0.31*scale/2, 0.22*scale/2, 0.0) )
                        {
                            scalar += ( sl_bold ) ? -0.2 : -0.02;
                        }
                        // d
                        if ( is_ellipse(x-centre, y-centre, 18.0, 0.16*scale/2, 0.41*scale/2, -0.22*scale/2, 0.0) )
                        {
                            scalar += ( sl_bold ) ? -0.2 : -0.02;
                        }
                        // e
                        if ( is_ellipse(x-centre, y-centre, 0.0, 0.21*scale/2, 0.25*scale/2, 0.0, 0.35*scale/2) )
                        {
                            scalar += ( sl_bold ) ? 0.1 : 0.01;
                        }
                        // f
                        if ( is_ellipse(x-centre, y-centre, 0.0, 0.046*scale/2, 0.046*scale/2, 0.0, 0.1*scale/2) )
                        {
                            scalar += ( sl_bold ) ? 0.1 : 0.01;
                        }
                        // g
                        if ( is_ellipse(x-centre, y-centre, 0.0, 0.046*scale/2, 0.046*scale/2, 0.0, -0.1*scale/2) )
                        {
                            scalar += ( sl_bold ) ? 0.2 : 0.02;
                        }
                        // h
                        if ( is_ellipse(x-centre, y-centre, 0.0, 0.046*scale/2, 0.023*scale/2, -0.08*scale/2, -0.605*scale/2) )
                        {
                            scalar += ( sl_bold ) ? 0.1 : 0.01;
                        }
                        // i
                        if ( is_ellipse(x-centre, y-centre, 0.0, 0.023*scale/2, 0.023*scale/2, 0.0, -0.605*scale/2) )
                        {
                            scalar += ( sl_bold ) ? 0.1 : 0.01;
                        }
                        // j
                        if ( is_ellipse(x-centre, y-centre, 0.0, 0.023*scale/2, 0.046*scale/2, 0.06*scale/2, -0.605*scale/2) )
                        {
                            scalar += ( sl_bold ) ? 0.1 : 0.01;
                        }
                    }
                    else
                    {
                        // a
                        if ( is_ellipsoid(x-centre, y-centre, z-centre, 0.0, 0.69*scale/2, 0.9*scale/2, 0.92*scale/2, 0.0, 0.0, 0.0) )
                        {
                            scalar += ( sl_bold ) ? 1.0 : 2.0;
                        }
                        // b
                        if ( is_ellipsoid(x-centre, y-centre, z-centre, 0.0, 0.6624*scale/2, 0.88*scale/2, 0.874*scale/2, 0.0, 0.0, -0.0184*scale/2) )
                        {
                            scalar += ( sl_bold ) ? -0.8 : -0.98;
                        }
                        // c
                        if ( is_ellipsoid(x-centre, y-centre, z-centre, -18.0, 0.11*scale/2, 0.22*scale/2, 0.31*scale/2, 0.22*scale/2, -0.25*scale/2, 0.0) )
                        {
                            scalar += ( sl_bold ) ? -0.2 : -0.02;
                        }
                        // d
                        if ( is_ellipsoid(x-centre, y-centre, z-centre, 18.0, 0.16*scale/2, 0.21*scale/2, 0.41*scale/2, -0.22*scale/2, -0.25*scale/2, 0.0) )
                        {
                            scalar += ( sl_bold ) ? -0.2 : -0.02;
                        }
                        // e
                        if ( is_ellipsoid(x-centre, y-centre, z-centre, 0.0, 0.21*scale/2, 0.35*scale/2, 0.25*scale/2, 0.0, -0.25*scale/2, 0.35*scale/2) )
                        {
                            scalar += ( sl_bold ) ? 0.1 : 0.01;
                        }
                        // f
                        if ( is_ellipsoid(x-centre, y-centre, z-centre, 0.0, 0.046*scale/2, 0.046*scale/2, 0.046*scale/2, 0.0, -0.25*scale/2, 0.1*scale/2) )
                        {
                            scalar += ( sl_bold ) ? 0.1 : 0.01;
                        }
                        // g
                        if ( is_ellipsoid(x-centre, y-centre, z-centre, 0.0, 0.046*scale/2, 0.046*scale/2, 0.046*scale/2, 0.0, -0.25*scale/2, -0.1*scale/2) )
                        {
                            scalar += ( sl_bold ) ? 0.2 : 0.02;
                        }
                        // h
                        if ( is_ellipsoid(x-centre, y-centre, z-centre, 0.0, 0.046*scale/2, 0.02*scale/2, 0.023*scale/2, -0.08*scale/2, -0.25*scale/2, -0.605*scale/2) )
                        {
                            scalar += ( sl_bold ) ? 0.1 : 0.01;
                        }
                        // i
                        if ( is_ellipsoid(x-centre, y-centre, z-centre, 0.0, 0.023*scale/2, 0.023*scale/2, 0.023*scale/2, 0.0, -0.25*scale/2, -0.605*scale/2) )
                        {
                            scalar += ( sl_bold ) ? 0.1 : 0.01;
                        }
                        // j
                        if ( is_ellipsoid(x-centre, y-centre, z-centre, 0.0, 0.023*scale/2, 0.02*scale/2, 0.046*scale/2, 0.06*scale/2, -0.25*scale/2, -0.605*scale/2) )
                        {
                            scalar += ( sl_bold ) ? 0.1 : 0.01;
                        }
                        // k
                        if ( is_ellipsoid(x-centre, y-centre, z-centre, 0.0, 0.04*scale/2, 0.1*scale/2, 0.056*scale/2, 0.06*scale/2, 0.625*scale/2, -0.105*scale/2) )
                        {
                            scalar += ( sl_bold ) ? 0.2 : 0.02;
                        }
                        // l
                        if ( is_ellipsoid(x-centre, y-centre, z-centre, 0.0, 0.056*scale/2, 0.1*scale/2, 0.056*scale/2, 0.0, 0.625*scale/2, 0.1*scale/2) )
                        {
                            scalar += ( sl_bold ) ? -0.2 : -0.02;
                        }
                    }
                }
                else if ( mode == "stark" )
                {
                    skeep = false;
                    if ( dim == 2 )
                    {
                        if ( // test triangle
                            ( x <= 0.2*scale && y <= 0.2*scale ) &&
                            ( y >= 0.35*scale - x )
                            )
                        {
                            data[counter] = 1;
                        }
                        else if ( // circle-1
                                 ( std::pow(x-0.1484375*scale, 2) + std::pow(y-0.5*scale, 2) ) <= std::pow(0.0234375*scale, 2) )
                        {
                            data[counter] = 0.75;
                        }
                        else if ( // circle-2
                                 ( std::pow(x-0.171875*scale, 2) + std::pow(y-0.5*scale, 2) ) <= std::pow(0.046875*scale, 2) )
                        {
                            data[counter] = 0.6;
                        }
                        else if ( // circle-3
                                 ( std::pow(x-0.21875*scale, 2) + std::pow(y-0.5*scale, 2) ) <= std::pow(0.09375*scale, 2) )
                        {
                            data[counter] = 0.45;
                        }
                        else if ( // circle-4
                                 ( std::pow(x-0.3125*scale, 2) + std::pow(y-0.5*scale, 2) ) <= std::pow(0.1875*scale, 2) )
                        {
                            data[counter] = 0.3;
                        }
                        else if ( // circle-5
                                 ( std::pow(x-0.5*scale, 2) + std::pow(y-0.5*scale, 2) ) <= std::pow(0.375*scale, 2) )
                        {
                            data[counter] = 0.15;
                        }
                        else //
                        {
                            data[counter] = 0;
                        }
                    }
                    else // 3d
                    {
                        if ( // test tetrahedron
                            ( x <= 0.2*scale && y <= 0.2*scale && z <= 0.2*scale ) &&
                            ( ( x - 0.175*scale + y - 0.175*scale + z - 0.175*scale ) >= 0.0 )
                            )
                        {
                            data[counter] = 1;
                        }
                        else if ( // sphere-1
                                 ( std::pow(x-0.1484375*scale, 2) + std::pow(y-0.5*scale, 2) + std::pow(z-0.5*scale, 2) ) <= std::pow(0.0234375*scale, 2) )
                        {
                            data[counter] = 0.75;
                        }
                        else if ( // sphere-2
                                 ( std::pow(x-0.171875*scale, 2) + std::pow(y-0.5*scale, 2) + std::pow(z-0.5*scale, 2) ) <= std::pow(0.046875*scale, 2) )
                        {
                            data[counter] = 0.6;
                        }
                        else if ( // sphere-3
                                 ( std::pow(x-0.21875*scale, 2) + std::pow(y-0.5*scale, 2) + std::pow(z-0.5*scale, 2) ) <= std::pow(0.09375*scale, 2) )
                        {
                            data[counter] = 0.45;
                        }
                        else if ( // sphere-4
                                 ( std::pow(x-0.3125*scale, 2) + std::pow(y-0.5*scale, 2) + std::pow(z-0.5*scale, 2) ) <= std::pow(0.1875*scale, 2) )
                        {
                            data[counter] = 0.3;
                        }
                        else if ( // sphere-5
                                 ( std::pow(x-0.5*scale, 2) + std::pow(y-0.5*scale, 2) + std::pow(z-0.5*scale, 2) ) <= std::pow(0.375*scale, 2) )
                        {
                            data[counter] = 0.15;
                        }
                        else //
                        {
                            data[counter] = 0;
                        }
                    }
                }
                else if ( mode == "dorn" ) // Oliver Dorn, et al., in: "Shape reconstruction in 2D from limited-view multifrequency electromagnetic data", 2001
                {
                    skeep = false;
                    if ( dim == 2 )
                    {
                        if ( // test house
                            ( 0.3*scale <= x && x <= 0.41*scale && // left side
                             0.1*scale <= y && y <= 0.5*scale ) ||
                            ( 0.59*scale <= x && x <= 0.7*scale && // right side
                             0.1*scale <= y && y <= 0.5*scale ) ||
                            ( 0.3*scale <= x && x <= 0.7*scale && // bottom side
                             0.1*scale <= y && y <= 0.21*scale ) ||
                            ( 0.3*scale <= x && x <= 0.7*scale && // top side
                             0.39*scale <= y && y <= 0.5*scale ) ||
                            ( 0.2*scale <= x && x <= 0.8*scale && // roof
                             0.5*scale <= y && y <= (m*x+c1) && y <= ((-m)*x+c2) ) ||
                            ( 0.3*scale <= x && x <= 0.41*scale && // chimney
                             0.5*scale <= y && y <= 0.9*scale )
                            )
                        {
                            data[counter] = 1;
                        }
                        else //
                        {
                            data[counter] = 0;
                        }
                    }
                    else // 3d
                    {
                        if ( // test house
                            ( 0.3*scale <= x && x <= 0.41*scale && // left side
                             0.1*scale <= y && y <= 0.5*scale &&
                             0.3*scale <= z && z < 0.7*scale ) ||
                            ( 0.59*scale <= x && x <= 0.7*scale && // right side
                             0.1*scale <= y && y <= 0.5*scale &&
                             0.3*scale <= z && z < 0.7*scale ) ||
                            ( 0.3*scale <= x && x <= 0.7*scale && // bottom side
                             0.1*scale <= y && y <= 0.21*scale &&
                             0.3*scale <= z && z < 0.7*scale ) ||
                            ( 0.3*scale <= x && x <= 0.7*scale && // top side
                             0.39*scale <= y && y <= 0.5*scale &&
                             0.3*scale <= z && z < 0.7*scale ) ||
                            ( 0.2*scale <= x && x <= 0.8*scale && // roof
                             0.5*scale <= y && y <= (m*x+c1) && y <= ((-m)*x+c2) && y <= (m*z+c1) && y <= ((-m)*z+c2) &&
                             0.2*scale <= z && z <= 0.8*scale ) ||
                            ( 0.3*scale <= x && x <= 0.41*scale && // chimney
                             0.5*scale <= y && y <= 0.9*scale &&
                             0.59*scale <= z && z <= 0.7*scale )
                            )
                        {
                            data[counter] = 1;
                        }
                        else //
                        {
                            data[counter] = 0;
                        }
                    }
                }
                else
                    // scalar
                {
                    data[counter] = 1;
                }
                if ( skeep )
                {
                    if ( -std::numeric_limits<double>::epsilon() < scalar && scalar < std::numeric_limits<double>::epsilon() )
                    {
                        scalar = double(0);
                        data[counter] = 0;
                    }
                    else
                    {
                        data[counter] = scalar;
                    }
                }
                counter++;
            }
        }
    }
    
    if ( dim == 2 )
    {
        std::valarray<size_t> Nx(sizeofgrid, 1);
        std::valarray<size_t> Ny(sizeofgrid, 1);
        std::valarray<size_t> Nz(1, 1);
        Plot3D::to_PNG(X, Y, Z, Nx, Ny, Nz, data, nameoffile + sext + ".png", Plot3D::NODES);
    }
    else if ( dim == 3 )
    {
        std::valarray<float> d0(float(0), sizeofgrid * sizeofgrid * s3);
        Plot3D::write_data(data, d0, d0, d0, d0, sizeofgrid, sizeofgrid, s3, nameoffile + sext + ".PBS", Binary);
        Plot3D::write_var(nameoffile + sext + ".VAR", mode, "zero_1", "zero_2", "zero_3", "zero_4");
    }
    
    return 0;
}
