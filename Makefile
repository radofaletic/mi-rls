# Makefiles for Rado Faletic's tomographic routines
# May 2001
# 			Department of Physics
# 			Faculty of Science
# 			Australian National University
# 			Canberra ACT 0200
# 			Australia
#

.PHONY: all clean tidy
.SUFFIXES:
.SUFFIXES: .cpp .o .h

srcdir = ./src
incdir = ./includes
testsrc = ./src/tests
prefix = .
bindir = $(prefix)
libdir = $(prefix)
mandir = $(prefix)

#
# ------------------------------------------------ #
# -------------------- defines ------------------- #
# ------------------------------------------------ #
#
# DEBUG is a flag for debugging, below are features
#       of each debugging level
#       0 - no debugging options at all (default)
#       1 - crucial debugging only (such as memory
#           allocation checking)
#       2 - (1) + crucial markers and descriptors
#       3 - (2) + I/O messages
#       4 - (3) + calculation messages
#	5 - all debugging information
DEBUG = #-DDEBUG=5 #-DUSE_MESSAGES
# GUI is either defined or not. If it is defined
#     then a GUI will be compiled into the program
GUI = #-DGUI
#
#
DEFINES = $(REAL) $(DEBUG) $(GUI)
#
# ------------------------------------------------ #
# --------------- Compiler options --------------- #
# ------------------------------------------------ #
#
CXXFLAGS = -O $(DEFINES) #-Wall
CPPFLAGS = -I. -I$(incdir)
LDFLAGS = -L.
LDLIBS = -lfftw3 -lm -lpng

# -------------------------------------------------------------- #
# ----------------- DO NOT EDIT BELOW THIS LINE ---------------- #
# -------------------------------------------------------------- #

#
# ------------------------------------------------ #
# ----------------- Make routines ---------------- #
# ------------------------------------------------ #
#

all: tidy

clean: clean_tests clean_progs tidy
	@echo "removing object files"
	@rm -f *.o

tidy: tidy_includes
	@echo "removing old and temporary files"
	@rm -rf ii_files rii_files *~ core

%.o:
	@echo -n "compiling '$@': "
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -c $< -o $@

%: %.o
	@echo -n "creating '$@': "
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) $(LDFLAGS) $< $(OBJS) $(LDLIBS) -o $@

#
# ------------------------------------------------ #
# -------------------- targets ------------------- #
# ------------------------------------------------ #
#

include includes/Makefile

include src/tests/Makefile

include src/Makefile
