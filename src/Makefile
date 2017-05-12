
#
# ------------------------- #
# fully featured programmes #
# ------------------------- #
#

.PHONY: progs tidy_progs clean_progs

PROGS_EXECS = cfdrc_convert plot3d_convert plot3dpng sinogram projection \
	fft2d inv2d projinv2d fft3d inv3d projinv3d png2proj ieee2png \
	inversion inversion2d halvepng error error3d png2plot3d errorise \
	xy2rt flatten transparent cutter pullout expand_csv stat scaledata \
	getreal 100pc

OBJS = 

$(PROGS_EXECS): $(OBJS)

progs: $(PROGS_EXECS)


tidy_progs:
	@echo "removing old and temporary programme files"
	@rm -rf $(srcdir)/*~

clean_progs:
	@echo "removing programme executables"
	@rm -f $(PROGS_EXECS) projected* saved_data.* *.messages

cfdrc_convert.o: $(srcdir)/cfdrc_convert.cpp \
	$(IARGV) \
	$(IFILE) \
	$(IFRONT-END) \
	$(IPLOT3D)

plot3d_convert.o: $(srcdir)/plot3d_convert.cpp \
	$(IARGV) \
	$(IFILE) \
	$(IFRONT-END) \
	$(IPLOT3D)

plot3dpng.o: $(srcdir)/plot3dpng.cpp \
	$(IARGV) \
	$(IFILE) \
	$(IFRONT-END) \
	$(IPLOT3D)

sinogram.o: $(srcdir)/sinogram.cpp \
	$(IANGLES) \
	$(IARGV) \
	$(IFRONT-END) \
	$(IGRID) \
	$(ILINE) \
	$(ISPARSE-MATRIX) \
	$(ITOMOGRAPHY)

projection.o: $(srcdir)/projection.cpp \
	$(IANGLES) \
	$(IARGV) \
	$(IFRONT-END) \
	$(IGRID) \
	$(ILINE) \
	$(ISPARSE-MATRIX) \
	$(ITOMOGRAPHY)

fft2d.o: $(srcdir)/fft2d.cpp \
	$(IANGLES) \
	$(IARGV) \
	$(IFRONT-END) \
	$(IGRID) \
	$(IINTERPOLATION) \
	$(ILINE) \
	$(ISPARSE_MATRIX) \
	$(ITOMOGRAPHY)

inv2d.o: $(srcdir)/inv2d.cpp \
	$(IANGLES) \
	$(IARGV) \
	$(IFRONT-END) \
	$(IGRID) \
	$(ILINE) \
	$(ILSQR) \
	$(IMATRIX_UTILITIES) \
	$(ISPARSE_MATRIX) \
	$(ITOMOGRAPHY)

projinv2d.o: $(srcdir)/projinv2d.cpp \
	$(IANGLES) \
	$(IARGV) \
	$(IFRONT-END) \
	$(IGRID) \
	$(IINTERPOLATION) \
	$(ILINE) \
	$(ILSQR) \
	$(IMATRIX_UTILITIES) \
	$(ISPARSE_MATRIX) \
	$(ITOMOGRAPHY)

fft3d.o: $(srcdir)/fft3d.cpp \
	$(IANGLES) \
	$(IARGV) \
	$(IFRONT-END) \
	$(IGRID) \
	$(IINTERPOLATION) \
	$(ILINE) \
	$(ISPARSE_MATRIX) \
	$(ITOMOGRAPHY)

inv3d.o: $(srcdir)/inv3d.cpp \
	$(IANGLES) \
	$(IARGV) \
	$(IFRONT-END) \
	$(IGRID) \
	$(ILINE) \
	$(ILSQR) \
	$(IMATRIX_UTILITIES) \
	$(ISPARSE_MATRIX) \
	$(ITOMOGRAPHY)

projinv3d.o: $(srcdir)/projinv3d.cpp \
	$(IANGLES) \
	$(IARGV) \
	$(IFRONT-END) \
	$(IGRID) \
	$(IINTERPOLATION) \
	$(ILINE) \
	$(ILSQR) \
	$(IMATRIX_UTILITIES) \
	$(ISPARSE_MATRIX) \
	$(ITOMOGRAPHY)

png2proj.o: $(srcdir)/png2proj.cpp \
	$(IANGLES) \
	$(ITOMOGRAPHY)

ieee2png.o: $(srcdir)/ieee2png.cpp \
	$(IANGLES) \
	$(IFILE) \
	$(IFORTRAN_IO) \
	$(ITOMOGRAPHY)

inversion.o: $(srcdir)/inversion.cpp \
	$(IANGLES) \
	$(IARGV) \
	$(IFRONT-END) \
	$(IGRID) \
	$(IINTERPOLATION) \
	$(ILINE) \
	$(ILSQR) \
	$(IMATRIX_UTILITIES) \
	$(ISPARSE_MATRIX) \
	$(ITOMOGRAPHY)

inversion2d.o: $(srcdir)/inversion2d.cpp \
	$(IANGLES) \
	$(IARGV) \
	$(IFRONT-END) \
	$(IGRID) \
	$(IINTERPOLATION) \
	$(ILINE) \
	$(ILSQR) \
	$(IMATRIX_UTILITIES) \
	$(ISPARSE_MATRIX) \
	$(ITOMOGRAPHY)

halvepng.o: $(srcdir)/halvepng.cpp \
	$(IANGLES) \
	$(ITOMOGRAPHY)

error.o: $(srcdir)/error.cpp \
	$(IANGLES) \
	$(IARGV) \
	$(IEXTRA_MATH) \
	$(IFILE) \
	$(IFRONT-END) \
	$(ITOMOGRAPHY)

error3d.o: $(srcdir)/error3d.cpp \
	$(IARGV) \
	$(IEXTRA_MATH) \
	$(IFILE) \
	$(IFRONT-END) \
	$(IPLOT3D)

png2plot3d.o: $(srcdir)/png2plot3d.cpp \
	$(IANGLES) \
	$(IGRID) \
	$(IPLOT3D) \
	$(ITOMOGRAPHY)

errorise.o: $(srcdir)/errorise.cpp \
	$(IANGLES) \
	$(IGAUSSIAN) \
	$(ITOMOGRAPHY)

xy2rt.o: $(srcdir)/xy2rt.cpp \
	$(IANGLES) \
	$(IARGV) \
	$(IFRONT-END) \
	$(ITOMOGRAPHY)

flatten.o: $(srcdir)/flatten.cpp \
	$(IANGLES) \
	$(ITOMOGRAPHY)

transparent.o: $(srcdir)/transparent.cpp \
	$(IANGLES) \
	$(ITOMOGRAPHY)

cutter.o: $(srcdir)/cutter.cpp \
	$(IANGLES) \
	$(ITOMOGRAPHY)

pullout.o: $(srcdir)/pullout.cpp \
	$(IANGLES) \
	$(ICONVERSIONS) \
	$(ITOMOGRAPHY)

expand_csv.o: $(srcdir)/expand_csv.cpp \
	$(IARGV) \
	$(IFILE) \
	$(IFRONT-END)

stat.o: $(srcdir)/stat.cpp \
	$(IANGLES) \
	$(ITOMOGRAPHY)

scaledata.o: $(srcdir)/scaledata.cpp \
	$(IANGLES) \
	$(ITOMOGRAPHY)

getreal.o: $(srcdir)/getreal.cpp \
	$(IANGLES) \
	$(IARGV) \
	$(IFRONT-END) \
	$(ITOMOGRAPHY)

100pc.o: $(srcdir)/100pc.cpp

