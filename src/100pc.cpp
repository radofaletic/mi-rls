
#include <fstream>
#include <iostream>
#include <string>

int main(int argc, char* argv[])
{
	std::ifstream input(argv[1]);
	std::ofstream output(argv[2]);
	
	size_t counter = 0;
	
	std::string buf;
	while ( input >> buf )
	{
		if ( !((counter++)%100) ) output << buf << std::endl;
	}
}
