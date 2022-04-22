/**
 interpolation
 
 interpolation functions
 
 Rado Faletic
 Department of Physics
 Faculty of Science
 Australian National University  ACT  0200
 Australia
 
 Rado.Faletic@anu.edu.au
 12th September 2004
 19th April 2022
 */





#ifndef _INTERPOLATION_
#define _INTERPOLATION_





/* ---------- standard header files ---------- */
#include <valarray>
#include <vector>





/* ---------- user header files ---------- */
#include "angles.h"
#include "extra_math.h"
#include "nngridr.h"
#include "tomography.h"





enum interpolation_method { BILINEAR, NN };





/* ---------- function declarations ---------- */

template<class T> void bilinear_interpolate(std::valarray<T>&, const std::valarray<T>&, const std::size_t& = 0, const std::size_t& = 0, const bool& = false);

template<class T> void nn_interpolate(std::valarray<T>&, const std::valarray<T>&, const std::size_t& = 0, const std::size_t& = 0);





/* ---------- function definitions ---------- */





/* ---------- bilinear_interpolate ---------- */
template<class T> void bilinear_interpolate(std::valarray<T>& udata, const std::valarray<T>& iangles,
                                            const std::size_t& xres, const std::size_t& yres, const bool& extend)
{
    // create 360 angles from 180 angles
    std::valarray<T> angles(2*iangles.size()+1);
    angles[std::slice(0,iangles.size(),1)] = iangles;
    angles[std::slice(iangles.size(),iangles.size(),1)] = iangles + T(180);
    angles[angles.size()-1] = T(360) + angles[0];
    
    // create 360 data from 180 data
    std::size_t oresolution = udata.size() / iangles.size();
    std::size_t DC = ( !(oresolution%2) ) ? (oresolution-2)/2 : (oresolution-1)/2;
    std::size_t resolution = oresolution - DC;
    std::valarray<T> data(angles.size()*resolution);
    for (std::size_t i=0; i<iangles.size(); i++)
    {
        // zero & positive components
        data[std::slice(i * resolution, resolution, 1)] = udata[std::slice(i * oresolution + DC, resolution, 1)];
        // zero components for negative side
        data[(i + iangles.size()) * resolution] = udata[i * oresolution + DC];
        // negative components
        for (std::size_t j=1; j<=DC; j++)
        {
            data[(i + iangles.size()) * resolution + j] = udata[i * oresolution + DC - j];
        }
        if ( !(oresolution%2) ) // even
        {
            data[(i + iangles.size()) * resolution + DC + 1] = udata[i * oresolution];
        }
    }
    // ... come full circle
    data[std::slice(data.size() - resolution, resolution, 1)] = data[std::slice(0, resolution, 1)];
    
    // create Cartesian points (x,y) for interpolated values
    std::size_t iresolution = ( xres != oresolution ) ? xres : oresolution;
    std::size_t jresolution = ( yres != oresolution ) ? yres : oresolution;
    std::vector< std::valarray<T> > points(iresolution * jresolution, std::valarray<T>(2));
    for (std::size_t j=0; j<jresolution; j++)
    {
        for (std::size_t i=0; i<iresolution; i++)
        {
            points[j*iresolution+i][0] = ( iresolution != oresolution ) ? (-T(DC) + T(i*oresolution)/T(iresolution))/std::sqrt(T(2)) : -T(DC) + T(i); // x
            points[j*iresolution+i][1] = ( jresolution != oresolution ) ? (-T(DC) + T(j*oresolution)/T(jresolution))/std::sqrt(T(2)) : -T(DC) + T(j); // y
        }
    }
    // now turn these points into polar format (r,theta)
    std::for_each(points.begin(), points.end(), Angle::C2p<T>);
    
    // create new interpolated data array
    udata.resize(points.size(), T(0));
    
    // interpolate
    for (std::size_t i=0; i<points.size(); i++)
    {
        for (std::size_t theta=0; theta<angles.size()-1; theta++)
        {
            if ( eq(points[i][1], angles[theta]) ) // point lies along one of the radial axes
            {
                for (std::size_t r=0; r<resolution-1; r++)
                {
                    if ( eq(points[i][0], T(r)) ) // point lies on old radius
                    {
                        udata[i] = data[theta*resolution+r];
                        break;
                    }
                    else if ( T(r) < points[i][0] && points[i][0] < T(r+1) ) // points lies between old radii
                    {
                        T p1 = points[i][0] - T(r);
                        T p2 = T(r+1) - points[i][0];
                        udata[i] = p1 * data[theta*resolution+r] + p2 * data[theta*resolution+r+1];
                        break;
                    }
                }
                if ( eq(points[i][0], T(resolution-1)) ||
                    ( T(resolution-1) < points[i][0] && extend ) ) // point lies beyond outer radius
                {
                    udata[i] = data[theta*resolution+resolution-1];
                }
            }
            else if ( angles[theta] < points[i][1] && points[i][1] < angles[theta+1] ) // between radial axes
            {
                T a1 = ( points[i][1] - angles[theta] ) / ( angles[theta+1] - angles[theta]);
                T a2 = ( angles[theta+1] - points[i][1] ) / ( angles[theta+1] - angles[theta]);
                for (std::size_t r=0; r<resolution-1; r++)
                {
                    if ( eq(points[i][0], T(r)) ) // point lies on old radius
                    {
                        udata[i] = a1 * data[theta*resolution+r] + a2 * data[(theta+1)*resolution+r];
                        break;
                    }
                    else if ( T(r) < points[i][0] && points[i][0] < T(r+1) ) // points lies between old radii
                    {
                        T interp1 = a1 * data[theta*resolution+r] + a2 * data[(theta+1)*resolution+r];
                        T interp2 = a1 * data[theta*resolution+r+1] + a2 * data[(theta+1)*resolution+r+1];
                        T p1 = points[i][0] - T(r);
                        T p2 = T(r+1) - points[i][0];
                        udata[i] = p1 * interp1 + p2 * interp2;
                        break;
                    }
                }
                if ( eq(points[i][0], T(resolution-1)) ||
                    ( T(resolution-1) < points[i][0] && extend ) ) // point lies beyond outer radius
                {
                    udata[i] = a1 * data[theta*resolution+resolution-1] + a2 * data[(theta+1)*resolution+resolution-1];
                }
            }
        }
    }
}





/* ---------- nn_interpolate ---------- */
template<class T> void nn_interpolate(std::valarray<T>& udata, const std::valarray<T>& angles,
                                      const std::size_t& xres, const std::size_t& yres)
{
    std::size_t resolution = udata.size() / angles.size();
    std::size_t DC = ( !(resolution%2) ) ? (resolution-2)/2 : (resolution-1)/2;
    
    std::vector< std::valarray<T> > ipoints(resolution, std::valarray<T>(2));
    
    for (std::size_t i=0; i<resolution; i++)
    {
        ipoints[i][0] = T(i) - T(DC);
        ipoints[i][1] = 0;
    }
    std::vector< std::valarray<T> > tpoints(resolution*angles.size(), std::valarray<T>(2));
    
    Rotation<T> rot(2);
    for (std::size_t i=0; i<angles.size(); i++)
    {
        rot.reset(angles[i], Angle::XY);
        for (std::size_t j=0; j<resolution; j++)
        {
            tpoints[i*resolution+j] = rot.O(ipoints[j]);
        }
    }
    
    std::size_t Nrows = ( yres ) ? yres : resolution;
    std::size_t Ncols = ( xres ) ? xres : resolution;
    nngridr(tpoints, udata, Nrows, Ncols);
}





#endif /* _INTERPOLATION_ */
