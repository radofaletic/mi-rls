
#
# ------------------------------------------------------------ #
# raw source files, and compilation relationships between them #
# ------------------------------------------------------------ #
#

backup_includes: backup_tests backup_progs
	@echo "backing up includes to 'Backups/'"
	@cp $(incdir)/*.h Backups

tidy_includes: tidy_tests tidy_progs
	@echo "removing old and temporary source files"
	@rm -rf $(incdir)/*~

IANGLES =	$(incdir)/angles.h

IARGS =		$(incdir)/args.h

ICONVERSIONS =	$(incdir)/conversions.h

IEXTRA_MATH =	$(incdir)/extra_math.h

IFILE =		$(incdir)/file.h

ITEXT-FRONT =	$(incdir)/text-front.h \
		$(ICONVERSIONS)

IGAUSSIAN =     $(incdir)/gaussian.h

IGUI-FRONT =	

IFRONT-END =	$(incdir)/front-end.h \
		$(ITEXT-FRONT) \
		$(IGUI-FRONT)

IFORTRAN_IO = 	$(incdir)/fortran_io.h \
		$(IFILE) \
		$(IFRONT-END)

IIMGLIB =	$(incdir)/imglib.h \
		$(IANGLES) \
		$(ICONVERSIONS) \
		$(IEXTRA_MATH) \
		$(IFILE) \
		$(IFRONT-END) \
		$(IFORTRAN_IO)

ILSQR =         $(incdir)/lsqr.h

IMATRIX =	$(incdir)/matrix.h \
		$(ICONVERSIONS) \
		$(IEXTRA_MATH) \
		$(IFILE) \
		$(IFORTRAN_IO) \
		$(IFRONT-END)

INNGRIDR = 	$(incdir)/nngridr.h \
		$(IANGLES) \
		$(ICONVERSIONS) \
		$(IFRONT-END)

IINTERPOLATION = $(incdir)/interpolation.h \
		$(IANGLES) \
		$(INNGRIDR)

IROTATION =	$(incdir)/rotation.h \
		$(IANGLES) \
		$(IEXTRA_MATH) \
		$(IMATRIX)

ISPARSE_MATRIX =	$(incdir)/sparse_matrix.h \
		$(ICONVERSIONS) \
		$(IEXTRA_MATH) \
		$(IFILE) \
		$(IFORTRAN_IO) \
		$(IFRONT-END) \
		$(ILSQR) \
		$(IMATRIX)


ISHAPE =	$(incdir)/shape.h

ILINE =		$(incdir)/line.h \
		$(IEXTRA_MATH) \
		$(ISHAPE)

IPLANE =	$(incdir)/plane.h \
		$(IEXTRA_MATH) \
		$(ILINE) \
		$(ISHAPE)

ISEGMENT =	$(incdir)/segment.h \
		$(IEXTRA_MATH) \
		$(ILINE) \
		$(ISHAPE)

ITRIANGLE =	$(incdir)/triangle.h \
		$(IEXTRA_MATH) \
		$(ILINE) \
		$(IPLANE) \
		$(ISEGMENT) \
		$(ISHAPE)

ITETRAHEDRON =	$(incdir)/tetrahedron.h \
		$(IEXTRA_MATH) \
		$(ILINE) \
		$(IPLANE) \
		$(ISHAPE) \
		$(ITRIANGLE)

IPOLYGON =	$(incdir)/polygon.h \
		$(IEXTRA_MATH) \
		$(ILINE) \
		$(ISEGMENT) \
		$(ISHAPE) \
		$(ITETRAHEDRON) \
		$(ITRIANGLE)

IPOLYHEDRON =	$(incdir)/polyhedron.h \
		$(IEXTRA_MATH) \
		$(ILINE) \
		$(IPOLYGON) \
		$(ISHAPE) \
		$(ITETRAHEDRON)

IPLOT3D =	$(incdir)/plot3d.h \
		$(ICONVERSIONS) \
		$(IEXTRA_MATH) \
		$(IFILE) \
		$(IFORTRAN_IO) \
		$(IFRONT-END) \
		$(ILINE) \
		$(IPOLYGON) \
		$(IPOLYHEDRON) \
		$(ISHAPE) \
		$(ISEGMENT)

IGRID =		$(incdir)/grid.h \
		$(ICONVERSIONS) \
		$(IEXTRA_MATH) \
		$(IFILE) \
		$(IFRONT-END) \
		$(ILINE) \
		$(IPLOT3D) \
		$(ISHAPE)

IMATRIX_UTILITIES =	$(incdir)/matrix_utilities.h \
		$(IFRONT-END) \
		$(IGRID)\
		$(ISPARSE_MATRIX)

ITOMOGRAPHY =	$(incdir)/tomography.h \
		$(IANGLES) \
		$(ICONVERSIONS) \
		$(IFRONT-END) \
		$(IGRID) \
		$(ILINE) \
		$(IROTATION) \
		$(ISPARSE_MATRIX)
