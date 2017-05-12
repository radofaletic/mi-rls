
/* check whether this machine is little or big endian */

#include <iostream>
#include <string>

int main(int argc, char* argv[])
{
  
  int mynumber = 67305985;
  std::cout << "mynumber = " << mynumber << std::endl;

  char* c = (char*) &mynumber;

  int c1 = int(*(c+0));
  int c2 = int(*(c+1));
  int c3 = int(*(c+2));
  int c4 = int(*(c+3));
  std::cout << "c.1 = " << c1 << " * 1\n";
  std::cout << "c.2 = " << c2 << " * 256\n";
  std::cout << "c.3 = " << c3 << " * 65536\n";
  std::cout << "c.4 = " << c4 << " * 16777216\n";
  int mynew = (c1*1) + (c2*256) + (c3*65536) + (c4*16777216);
  int mynewbs = (c4*1) + (c3*256) + (c2*65536) + (c1*16777216);
  std::cout << " c  = " << mynew;
  if ( mynew == mynumber )
    {
      std::cout << " (LITTLE ENDIAN)";
    }
  std::cout << std::endl;
  std::cout << "cbs = " << mynewbs;
  if ( mynewbs == mynumber )
    {
      std::cout << " (BIG ENDIAN)";
    }
  std::cout << std::endl;
}
