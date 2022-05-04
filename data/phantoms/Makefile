
#CXX = g++ -mcpu=athlon-tbird -march=athlon-tbird -m3dnow
#CXX = g++ -mtune=core2
CXX = c++ -std=c++20

OPT = -I../../includes

default: phantoms

clean:
	rm -rf *~ phantoms

remove: clean
	rm -rf *.png *.PBG *.PBS *.VAR

phantoms: phantoms.cpp ../../includes/angles.h ../../includes/plot3d.h
	$(CXX) $(OPT) phantoms.cpp -lpng -o phantoms
