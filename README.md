# README #

Rado Faletič, May 2005

Contains the code for the C++ routines used for my PhD thesis “Tomographic Reconstruction of Shock Layer Flows”, which can be read at http://hdl.handle.net/1885/46916

The results are generally stored in plain PNG files, which should be viewable in any modern web browser or imaging programme. Note, however, that PNG images are integer. Also stored in these files, in textual IEEE number format are the decimal values of the results. These values can be retrieved with the Tomography::pngread function found in the include/tomography.h header file.

This function, along with Tomography::pngwrite, allow for reading and writing of decimal data in PNG files. They are also able to store multiple images (of the same dimensions) in one PNG file, which makes displaying a set of images simply a matter of viewing the PNG image.

Any queries about this work should be directed, initially, to Rado Faletic. Failing to contact Rado, you should try Frank Houwing. Failing that, contact the Department of Physics in the Faculty of Science at the Australian National University.

### byte_swap ###
> Contains routines for byte-swapping all of the standard C/C++ variable types. It has been verified to swap between SGI/IRIX and Intel/Linux platforms (the main motivation for this piece of code since I used an SGI at ANU and had an Intel machine at home

### delaunay (NOT WRITTEN) ###
> Routines for dealing with unstructured grids. In particular, it supports delaunay grids (much nicer than other structures)

### grid ###
> Here are the classes and associated functions for grids stored in system memory. I have defined a grid to consist primarily of nodes. Overlayed is then information about the cell structure, and this is dependant on the type of grid. For _fast_ cell navigation I have also included support for storing the cell neighbours, which will enable the use of the “cell walking” algorithm. The user may also include a-priori data over the grid cells, which can be used as initial iteration values, or residual calculations. Support is included to project a straight line through the grid. (NOT COMPLETED)

### lines_and_planes ###
> Exactly what the title suggests, code for dealing with lines and planes. These are different to the other geometric structures used here, since they have no defined boundary. A line is determined by a slope, and an initial point, and a (hyper)-plane is given by a normal vector and an initial point. These work in arbitrary dimension (at least two). Some of the functions include find the intersection of objects of these type, the determining if a point is contained in one of these objects, and finding the perpendicular distance to one of these objects, given any point in space.

### plot3d ###
> Special routines for handling structured (in particular Plot3D) data. I/O functions are given, which include support for passing data values pre and post calculation.

### shapes ###
> Here we have three classes: triangles, quadrilaterals, and tetrahedra. The later are purely 3D objects, but the others may exist in two or more dimensions. Some of the functions include intersecting with a line, determining if a point lies within an object, making special exception if the point lies on the boundary, finding the area/volume etc. These objects are used extensively throughout the code since all cell structures supported can be described by a combination of these objects.