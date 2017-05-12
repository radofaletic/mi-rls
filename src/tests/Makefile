
#
# ------------------------------------------------------------ #
# my own little test programs to see if each module works okay #
# ------------------------------------------------------------ #
#

.PHONY: clean_tests tests

TEST_EXECS = p3dio gio mtest rot_test invtest ffttest pngtest nngridr \
	testarg makeshape endian some_maths

OBJS = 

$(TEST_EXECS): $(OBJS)

tests: $(TEST_EXECS)


tidy_tests:
	@echo "removing old and temporary test files"
	@rm -rf $(testsrc)/*~

clean_tests:
	@echo "removing test executables"
	@rm -f $(TEST_EXECS) projected* saved_data.* *.messages pngtest.png


p3dio.o: $(testsrc)/p3dio.cpp \
	$(IFRONT-END) \
	$(IPLOT3D)

gio.o: $(testsrc)/gio.cpp \
	$(IFRONT-END) \
	$(IGRID)

mtest.o: $(testsrc)/mtest.cpp \
	$(ICONVERSIONS) \
	$(IFILE) \
	$(IFRONT-END) \
	$(ISPARSE_MATRIX)

rot_test.o: $(testsrc)/rot_test.cpp \
	$(IANGLES) \
	$(ICONVERSIONS) \
	$(IERROR) \
	$(IFRONT-END) \
	$(IROTATION)

invtest.o: $(testsrc)/invtest.cpp \
	$(ICONVERSIONS) \
	$(IFORTRAN_IO) \
	$(IFRONT-END) \
	$(IGRID) \
	$(ILINE) \
	$(IMATRIX_UTILITIES) \
	$(ISPARSE_MATRIX) \
	$(ITOMOGRAPHY)

ffttest.o: $(testsrc)/ffttest.cpp \
	$(IANGLES) \
	$(IARGV) \
	$(ICONVERSIONS) \
	$(IEXTRA_MATH) \
	$(IFORTRAN_IO) \
	$(IFRONT-END) \
	$(IGRID) \
	$(IIMGLIB) \
	$(ILINE) \
	$(INNGRIDR) \
	$(IROTATION)

pngtest.o: $(testsrc)/pngtest.cpp \
	$(IANGLES) \
	$(IIMGLIB)

nngridr.o: $(testsrc)/nngridr.cpp \
	$(IIMGLIB) \
	$(INNGRIDR)

testarg.o: $(testsrc)/testarg.cpp \
	$(IFRONT-END)

makeshape.o: $(testsrc)/makeshape.cpp \
	$(ITOMOGRAPHY)

endian.o: $(testsrc)/endian.cpp

some_maths.o: $(testsrc)/some_maths.cpp

