/*
 Rado's own little routine for testing
 Plot3d I/O through the grid interface
 */

#include <iostream>

#if DEBUG >= 5
#ifndef USE_MESSAGES
#define USE_MESSAGES
#endif /* USE_MESSAGES */
#endif /* DEBUG */

#include "front-end.h"
#include "grid.h"

#ifdef USE_MESSAGES
#include <fstream>
std::ofstream messages;
#endif /* USE_MESSAGES */

int main(int argc, char* argv[])
{
#ifdef USE_MESSAGES
	messages.open((std::string(argv[0])+std::string(".messages")).c_str());
#endif /* USE_MESSAGES */
	
	grid_input gridinputs;
	gridinputs.type() = structured;
	gridinputs.format() = Unformatted;
	gridinputs.precision() = Single;
	gridinputs.multidomain() = true;
	gridinputs.blanking() = true;
	gridinputs.gridfile() = "data/BF1/BF1.PUG";
	gridinputs.datafile() = "data/BF1/BF1.PUS";
	gridinputs.qdata() = 1;
	
	grid<float> test_grid(gridinputs);
	test_grid.read_data(gridinputs);
	
	std::string outputname = "gio_output";
	dataformat format = Formatted;
	test_grid.write(outputname,SaveGrid,format);
	test_grid.write(outputname,SaveData,format);
	test_grid.write(outputname,SaveNeighbours,format);
	
	message("FINISHED running "+ntos(argv[0]));
	
}
