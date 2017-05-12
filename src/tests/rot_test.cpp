/*
  for testing rotation
*/

#include <iostream>
#include <string>
#include <valarray>
#include <vector>

#if DEBUG >= 5
#ifndef USE_MESSAGES
#define USE_MESSAGES
#endif /* USE_MESSAGES */
#endif /* DEBUG */

#include "conversions.h"
#include "front-end.h"
#include "matrix.h"
#include "rotation.h"

#ifdef USE_MESSAGES
#include <fstream>
std::ofstream messages;
#endif /* USE_MESSAGES */

int main(int argc, char* argv[])
{
#ifdef USE_MESSAGES
  messages.open((std::string(argv[0])+std::string(".messages")).c_str());
#endif /* USE_MESSAGES */

  float angle = float(135);
  float length = float(1);
  std::valarray<float> origin(float(1),3);

  std::string tab = "\t\t\t\t";

  //
  // the initial point
  //
  std::valarray<float> point(float(0),3);
  if ( argc > 1 )
    {
      if ( std::string(argv[1]) == std::string("x") )
	{
	  message("creating initial point ("+ntos(length)+",0,0)...");
	  point[0] = length;
	}
      else if ( std::string(argv[1]) == std::string("y") )
	{
	  message("creating initial point (0,"+ntos(length)+",0)...");
	  point[1] = length;
	}
      else if ( std::string(argv[1]) == std::string("z") )
	{
	  message("creating initial point (0,0,"+ntos(length)+")...");
	  point[2] = length;
	}
      else if ( std::string(argv[1]) == std::string("xy") )
	{
	  message("creating initial point ("+ntos(length)+","+ntos(length)+",0)...");
	  point[0] = length;
	  point[1] = length;
	}
      else if ( std::string(argv[1]) == std::string("yz") )
	{
	  message("creating initial point (0,"+ntos(length)+","+ntos(length)+")...");
	  point[1] = length;
	  point[2] = length;
	}
      else if ( std::string(argv[1]) == std::string("zx") )
	{
	  message("creating initial point ("+ntos(length)+",0,"+ntos(length)+")...");
	  point[2] = length;
	  point[0] = length;
	}
      else if ( std::string(argv[1]) == std::string("xyz") )
	{
	  message("creating initial point ("+ntos(length)+","+ntos(length)+","+ntos(length)+")...");
	  point[0] = length;
	  point[1] = length;
	  point[2] = length;
	}
      else
	{
	  message("\n* U N K N O W N   S W I T C H *     "+std::string(argv[1]));
	}
    }
  else
    {
      message("creating initial point ("+ntos(length)+",0,0)...");
      point[0] = length;
    }
  message(tab+vtos(point));

  //
  // zero origin
  //
  message("\n*ZERO ORIGIN*");

  //
  // define the rotation element
  //
  message("creating rotation element...");
  Matrix<float> trot(3,3);
  trot(0,0) = cos(angle);
  trot(0,1) = -sin(angle);
  trot(1,0) = sin(angle);
  trot(1,1) = cos(angle);
  message(trot.print());
  Rotation<float> rot(3);
  rot.set(angle);
  message(ctos(typeid(rot).name())+"\n"+rot.print());

  //
  // rotate the initial point +90degrees on the xy-plane
  //
  std::valarray<float> tnewvec(3);
  tnewvec[0] = cos(angle);
  tnewvec[1] = sin(angle);
  message("rotating initial point +z "+ntos(angle)+"...\t"+vtos(tnewvec));
  std::valarray<float> newvec = rot(point);
  message(tab+vtos(newvec));

  //
  // rotate the initial point +90degrees on the zx-plane
  //
  tnewvec[0] = cos(angle);
  tnewvec[1] = float(0);
  tnewvec[2] = -sin(angle);
  message("rotating initial point +y "+ntos(angle)+"...\t"+vtos(tnewvec));
  rot.reset(angle,Angle::Y);
  newvec = rot(point);
  message(tab+vtos(newvec));

  //
  // rotate the initial point +90degrees on the yz-plane
  //
  tnewvec = point;
  message("rotating initial point +x "+ntos(angle)+"...\t"+vtos(tnewvec));
  rot.reset(angle,Angle::X);
  newvec = rot(point);
  message(tab+vtos(newvec));

  //
  // shifted origin
  //
  message("\n*SHIFTED ORIGIN*");

  //
  // set the new origin
  //
  message("setting rotation origin "+vtos(origin)+"...");
  rot.set_origin(origin);
  message(tab+rot.origin_print());

  //
  // rotate the initial point +90degrees on the xy-plane
  //
  tnewvec[0] = float(1)+sin(angle);
  tnewvec[1] = float(1)+sin(angle);
  tnewvec[2] = float(1);
  message("rotating initial point +z "+ntos(angle)+"...\t"+vtos(tnewvec));
  rot.reset(angle);
  newvec = rot(point);
  message(tab+vtos(newvec));


  long double atr = 90.0;
  std::vector< std::valarray<long double> > at(6);
  Rotation<long double> atrot(2);
  atrot.set(atr);
  std::cout << std::endl << std::endl
	    << "ROTATE A SERIES OF POINTS" << std::endl
	    << "initial series of " << at.size() << ":" << std::endl;
  for (size_t i=0; i<at.size(); i++)
    {
      at[i].resize(2);
      at[i][0] = (long double)(i) - ((long double)(at.size()-1)/2.0);
      at[i][1] = 0.0;
      std::cout << at[i][0] << "," << at[i][1] << std::endl;
    }
  std::cout << "series rotated by " << atr << " degrees:" << std::endl;
  for (size_t i=0; i<at.size(); i++)
    {
      at[i] = atrot.O(at[i]);
      std::cout << at[i][0] << "," << at[i][1] << std::endl;
    }
  message("\nFINISHED running "+std::string(argv[0]));
  return 0;
}
