
#if DEBUG >= 5
#ifndef USE_MESSAGES
#define USE_MESSAGES
#endif /* USE_MESSAGES */
#endif /* DEBUG */

#include <iostream>
#include "argv.h"

#ifdef USE_MESSAGES
#include <fstream>
std::ofstream messages;
#endif /* USE_MESSAGES */

int main(int argc, char* argv[])
{
  std::vector<args> my_args = get_args(argc, argv);
  std::cout << "my_args.size() = " << my_args.size() << std::endl;
  for (size_t i=0; i<my_args.size(); i++)
    {
      std::cout << "[" << i << "].var_ = " << my_args[i].var() << "\t"
		<< "[" << i << "].val_ = " << my_args[i].val() << std::endl;
    }
}
