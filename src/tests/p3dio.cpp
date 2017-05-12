/*
  Rado's own little routine for testing Plot3D I/O
*/

#include <iostream>
#include <valarray>
#include <vector>
#include <string>

#if DEBUG >= 5
#ifndef USE_MESSAGES
#define USE_MESSAGES
#endif /* USE_MESSAGES */
#endif /* DEBUG */

#include "front-end.h"
#include "plot3d.h"

#ifdef USE_MESSAGES
#include <fstream>
std::ofstream messages;
#endif /* USE_MESSAGES */

int main(int argc, char* argv[])
{
#ifdef USE_MESSAGES
  messages.open((std::string(argv[0])+std::string(".messages")).c_str());
#endif /* USE_MESSAGES */

  std::vector< std::valarray<float> > nodes;
  std::valarray<size_t> Nx;
  std::valarray<size_t> Ny;
  std::valarray<size_t> Nz;
  std::valarray<float> data;
  std::string gfilename;
  std::string dataname;
  dataformat format;
  dataprecision fsize;
  bool multidomain;

  Plot3D::read_interface(nodes, Nx, Ny, Nz, gfilename, format, fsize, multidomain);
  Plot3D::read_data_interface(data, dataname, Nx, Ny, Nz, gfilename, format, fsize, multidomain);

  std::string oname = "p3dio_output";

  Plot3D::write_interface(nodes, Nx, Ny, Nz, oname);
  Plot3D::write_data_interface(data, dataname, Nx, Ny, Nz, oname);

  message("FINISHED running "+ntos(argv[0]));
}
