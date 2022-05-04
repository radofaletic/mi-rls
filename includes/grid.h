/**
 grid
 
 A collection of C++ routine to handle various 2d and 3d grid structure.
 
 Rado Faletic
 Department of Physics
 Faculty of Science
 Australian National University  ACT  0200
 Australia
 
 Rado.Faletic@anu.edu.au
 
 7th November 2004
 19th April 2022, updated to C++20
 */

/**
 NOTES
 - a plot3d grid can be any IJK grid structure
 - specific reference to Plot3D refers to the specific file format used by Plot3D, to find out more about this format please read `Plot3D.txt' or `Plot3D.html'
 - cell 'neighbours', in this code, are defined as cells which share any node with the cell in question
 
 grid arrangements referenced from:
 
 CFD-ACE
 Theory Manual
 Chapter 1: Domain Modeling
 Version 4.0, February 1998
 CFD Research Corporation
 */





#ifndef _GRID_
#define _GRID_





/* ---------- standard header files ---------- */
#include <algorithm>
#include <cmath>
#include <cstring>
#include <ctime>
#include <ostream>
#include <sstream>
#include <string>
#include <typeinfo>
#include <valarray>
#include <vector>





/* ---------- user header files ---------- */
#include "conversions.h"
#include "extra_math.h"
#include "file.h"
#include "front-end.h"
#include "line.h"
#include "plot3d.h"
#include "shape.h"





/* ---------- grid specific enumerators ---------- */

/// the supported types of grid
enum grid_type
{
	empty,
	structured
};

/// items available to save
enum grid_saveitem
{
	SaveGrid,
	SaveData,
	SaveNeighbours
};

/// the supported sets of basis functions
enum grid_basis
{
	unity,
	l2norm
};

/// the supported ray-tracing methods
enum grid_proj_method
{
	brute,
	walk,
	walkfast
};





/* ---------- class & function declaration ---------- */

class grid_input
{
private:
	grid_type type_;
	bool load_grid_;
	dataformat format_;
	dataprecision precision_;
	bool multidomain_;
	bool blanking_;
	std::string gridfile_;
	std::string datafile_;
	unsigned short qdata_;
	std::size_t g_nX_;
    std::size_t g_nY_;
    std::size_t g_nZ_;
	double g_scale_;
public:
	grid_input() : type_(empty), load_grid_(false), format_(Formatted), precision_(Single), multidomain_(false), blanking_(false), gridfile_(""), datafile_(""), qdata_(1), g_nX_(0), g_nY_(0), g_nZ_(0), g_scale_(0) { };
	grid_type type() const { return this->type_; };
	grid_type& type() { return this->type_; };
	bool load_grid() const { return this->load_grid_; };
	bool& load_grid() { return this->load_grid_; };
	dataformat format() const { return this->format_; };
	dataformat& format() { return this->format_; };
	dataprecision precision() const { return this->precision_; };
	dataprecision& precision() { return this->precision_; };
	bool multidomain() const { return this->multidomain_; };
	bool& multidomain() { return this->multidomain_; };
	bool blanking() const { return this->blanking_; };
	bool& blanking() { return this->blanking_; };
	std::string gridfile() const { return this->gridfile_; };
	std::string& gridfile() { return this->gridfile_; };
	std::string datafile() const { return this->datafile_; };
	std::string& datafile() { return this->datafile_; };
	unsigned short qdata() const { return this->qdata_; };
	unsigned short& qdata()  { return this->qdata_; };
    std::size_t g_nX() const { return this->g_nX_; };
    std::size_t& g_nX() { return this->g_nX_; };
    std::size_t g_nY() const { return this->g_nY_; };
    std::size_t& g_nY() { return this->g_nY_; };
    std::size_t g_nZ() const { return this->g_nZ_; };
    std::size_t& g_nZ() { return this->g_nZ_; };
	double g_scale() const { return this->g_scale_; };
	double& g_scale() { return this->g_scale_; };
};

std::ostream& operator<<(std::ostream& s, const grid_input& g)
{
	std::ostringstream o;
	o << "type = ";
	switch (g.type())
	{
		case structured:
			o << "structured\n";
			break;
		default:
			o << "empty\n";
			break;
	}
	o << "load_grid = ";
	o << ( g.load_grid() ? "true" : "false\n" );
	o << "format = ";
	switch (g.format())
	{
		case Formatted:
			o << "Formatted\n";
			break;
		case Unformatted:
			o << "Unformatted\n";
			break;
		case Binary:
			o << "Binary\n";
			break;
	}
	o << "precision = ";
	switch (g.precision())
	{
		case Single:
			o << "Single\n";
			break;
		case Double:
			o << "Double\n";
			break;
	}
	o << "multidomain = ";
	o << ( g.multidomain() ? "on" : "off\n" );
	o << "blanking = ";
	o << ( g.blanking() ? "on" : "off\n" );
	o << "gridfile = " << g.gridfile() + '\n';
	o << "datafile = " << g.datafile() + '\n';
	o << "qdata = " << g.qdata() + '\n';
	o << "g_nX = " << g.g_nX() + '\n';
	o << "g_nY = " << g.g_nY() + '\n';
	o << "g_nZ = " << g.g_nZ() + '\n';
	o << "g_scale = " << g.g_scale() + '\n';
	
	return s << o.str();
};

template<class T> class grid_surface
{
private:
    std::size_t cellnumber_;            // proper grid cell that this surface structure belongs to
	const compact_shape<T>* cell_; // the surface structure
public:
	grid_surface() { };
	~grid_surface() { };
	grid_surface(const std::size_t& cn, const compact_shape<T>* c) : cellnumber_(cn), cell_(c) { };
    std::size_t cellnumber() const { return cellnumber_; };
	const compact_shape<T>* cell() const { return cell_; };
};

template<class T> class grid
{
private:
	// general grid quantities
	grid_type type_;                                  // the type of grid
	std::vector< std::valarray<T> > node_;            // grid nodes
	std::vector< compact_shape<T>* > cell_;           // grid cells
	std::valarray<T> measure_;                        // area/volume of each cell
	std::vector< std::valarray<T> > bnode_;           // nodes forming a bounding cell around the entire grid
	compact_shape<T>* bcell_;                         // bounding cell
	unsigned short dim_;                              // embeded dimension
	std::valarray<std::size_t> Nx_;                        // nodes in the x direction (structured)
	std::valarray<std::size_t> Ny_;                        // nodes in the y direction (structured)
	std::valarray<std::size_t> Nz_;                        // nodes in the z direction (structured)
	std::vector< std::valarray<std::size_t> > neighbour_;  // adjacent neighbours
	std::vector< grid_surface<T> > surface_;          // surface cells
	std::valarray<std::size_t> Ns_a_;                      // surface cells on the 0 and 3 faces (structured)
	std::valarray<std::size_t> Ns_b_;                      // surface cells on the 1 and 4 faces (structured)
	std::valarray<std::size_t> Ns_c_;                      // surface cells on the 2 and 5 faces (structured)
	std::vector< std::valarray<std::size_t> > surface_neighbour_;// surface neighbours
	std::valarray<T> data_;                           // data
	std::valarray<T> adata_;                          // apriori data
	bool aux_;                                        //
    std::size_t cell_walk_;                                // last cut cell
    std::size_t surface_walk_;                             // last cut on the surface
	T cscale_;                                        // characteristic scale length
	std::valarray<T> min_;                            // minimum bounds
	std::valarray<T> max_;                            // maximum bounds
	std::valarray<T> centroid_;                       // weighted center of the grid
	std::valarray<T> center_;                         // geometric center of the grid (max_+min_)/2
	grid_basis basis_;                                // type of measure on the cells
	// grid file quantities
	std::string dataname_;                            // name of the data
	std::string gfilename_;                           // grid filename
	dataformat format_;                               // file format
	dataprecision precision_;                         // size of floats used
	bool byte_swapping_;                              // boolean flag (byte swap)
	// private functions
	void generate_grid(const grid_input&);
	void generate_cells();
	void generate_measures();
	void generate_neighbours();
	void generate_surface();
	void generate_scale();
	void generate_bcell();
public:
    // user functions
	grid(const grid_input&, const bool& = true);
	void init(const grid_input&, const bool& = true);
	void auxs();
	~grid();
	void clear();
	void clear_cells();
	void read(const grid_input&);
	void write(const std::string& = "grid_output", const grid_saveitem& = SaveGrid, const dataformat& = Formatted) const;
	void read_data(const grid_input&);
	void put_adata(const std::valarray<T>&);
	std::valarray<T> get_adata() const;
	void write_data(const std::string& = "grid_output", const dataformat& = Formatted, const std::valarray<T>& = std::valarray<T>(0)) const;
	bool check() const;
	void aux();
	void clear_aux();
	unsigned short dim() const;
    std::size_t ncells() const;
	std::valarray<T> center() const;
	std::valarray<T> wcenter() const;
	std::valarray<T> min() const;
	std::valarray<T> max() const;
	T raytrace(const line<T>&, std::vector< two_numbers<T> >&, const grid_proj_method& = walkfast);
	void shift(const std::valarray<T>&);
	void give_dataname(const std::string&);
	std::string get_dataname() const;
	std::string get_filename() const;
	T& operator[] (const std::size_t&);
	T operator[] (const std::size_t&) const;
	void clear_data();
	void set_basis(const grid_basis&);
	T applybasis(const std::size_t&, const T& = 1) const;
	T scale() const;
	std::valarray<std::size_t> get_neighbours(const std::size_t&) const;
    std::size_t nex() const { return this->Nx_[0]; };
    std::size_t ney() const { return this->Ny_[0]; };
    std::size_t nez() const { return this->Nz_[0]; };
};





/* ---------- function definitions ---------- */





/* ---------- generate_grid ---------- */
template<class T> void grid<T>::generate_grid(const grid_input& thisgrid)
{
	this->dim_ = 0;
	if ( thisgrid.g_nX() )
	{
		this->dim_++;
	}
	if ( thisgrid.g_nY() )
	{
		this->dim_++;
	}
	if ( thisgrid.g_nZ() )
	{
		this->dim_++;
	}
	if ( this->dim_ == 0 )
	{
		return;
	}
	this->type_ = thisgrid.type();
	this->format_ = thisgrid.format();
	this->precision_ = thisgrid.precision();
	this->gfilename_ = thisgrid.gridfile();
	this->Nx_.resize(1);
	this->Ny_.resize(1);
	this->Nz_.resize(1);
	this->Nx_[0] = 1;
	this->Ny_[0] = 1;
	this->Nz_[0] = 1;
	this->cscale_ = std::abs(thisgrid.g_scale()) ? std::abs(T(thisgrid.g_scale())) : T(1);
	if ( this->dim_ == 2 )
	{
		if ( thisgrid.g_nX() && thisgrid.g_nY() )
		{
			this->Nx_[0] += thisgrid.g_nX();
			this->Ny_[0] += thisgrid.g_nY();
		}
		else if ( thisgrid.g_nY() && thisgrid.g_nZ() )
		{
			this->Nx_[0] += thisgrid.g_nY();
			this->Ny_[0] += thisgrid.g_nZ();
		}
		else if ( thisgrid.g_nZ() && thisgrid.g_nX() )
		{
			this->Nx_[0] += thisgrid.g_nZ();
			this->Ny_[0] += thisgrid.g_nX();
		}
		this->node_.resize(this->Nx_[0]*this->Ny_[0]);
		for (std::size_t j=0; j<this->Ny_[0]; j++)
		{
			for (std::size_t i=0; i<this->Nx_[0]; i++)
			{
				this->node_[j*this->Nx_[0]+i].resize(2);
				this->node_[j*this->Nx_[0]+i][0] = i * this->cscale_;
				this->node_[j*this->Nx_[0]+i][1] = j * this->cscale_;
			}
		}
	}
	else if ( this->dim_ == 3 )
	{
		this->Nx_[0] += thisgrid.g_nX();
		this->Ny_[0] += thisgrid.g_nY();
		this->Nz_[0] += thisgrid.g_nZ();
		this->node_.resize(this->Nx_[0]*this->Ny_[0]*this->Nz_[0]);
		for (std::size_t k=0; k<this->Nz_[0]; k++)
		{
			for (std::size_t j=0; j<this->Ny_[0]; j++)
			{
				for (std::size_t i=0; i<this->Nx_[0]; i++)
				{
					node_[k*this->Nx_[0]*this->Ny_[0]+j*this->Nx_[0]+i].resize(3);
					node_[k*this->Nx_[0]*this->Ny_[0]+j*this->Nx_[0]+i][0] = i * this->cscale_;
					node_[k*this->Nx_[0]*this->Ny_[0]+j*this->Nx_[0]+i][1] = j * this->cscale_;
					node_[k*this->Nx_[0]*this->Ny_[0]+j*this->Nx_[0]+i][2] = k * this->cscale_;
				}
			}
		}
	}
}





/* ---------- generate_cells ---------- */
template<class T> void grid<T>::generate_cells()
{
	this->cell_.clear();
    std::size_t total_cells = 0;
	switch(this->type_)
	{
        case structured: {
            std::size_t z = 0;
			for (std::size_t i=0; i<Nx_.size(); i++)
			{
				z += ( ( this->Nx_[i] < 3 ) ? 1 : this->Nx_[i] - 1 )
				* ( ( this->Ny_[i] < 3 ) ? 1 : this->Ny_[i] - 1 )
				* ( ( this->Nz_[i] < 3 ) ? 1 : this->Nz_[i] - 1 );
			}
			std::string vtxt = " " + std::string(typeid(compact_shape<T>).name());
			switch(this->dim_)
			{
				case 2:
					vtxt = " " + std::string(typeid(polygon<T>).name());
					break;
				case 3:
					vtxt = " " + std::string(typeid(polyhedron<T>).name());
					break;
			}
			vtxt += " cells";
			total_cells = z;
			debug("grid<T>::generate_cells", "generating " + std::to_string(total_cells) + vtxt);
			
			this->cell_.resize(z);
			z = 0;
            std::size_t pz = 0;
			for (std::size_t d=0; d<this->Nx_.size(); d++)
			{
                std::size_t nnx = ( this->Nx_[d] < 3 ) ? 1 : this->Nx_[d] - 1;
                std::size_t nny = ( this->Ny_[d] < 3 ) ? 1 : this->Ny_[d] - 1;
                std::size_t nnz = ( this->Nz_[d] < 3 ) ? 1 : this->Nz_[d] - 1;
				for (std::size_t k=0; k<nnz; k++)
				{
					for (std::size_t j=0; j<nny; j++)
					{
						for (std::size_t i=0; i<nnx; i++)
						{
                            std::size_t cell = z + k * nnx * nny + j * nnx + i;
                            std::size_t node = pz + k * Nx_[d] * Ny_[d] + j * Nx_[d] + i;
							// assign the vertices for the cell
							switch(this->dim_)
							{
								case 2: // polygon
									this->cell_[cell] = new polygon<T>(&this->node_[node], &this->node_[node+1], &this->node_[node+this->Nx_[d]+1], &this->node_[node+this->Nx_[d]]);
									break;
								case 3: // polyhedron
                                    std::size_t nextnode = node + this->Nx_[d] * this->Ny_[d];
									this->cell_[cell] = new polyhedron<T>(&this->node_[node], &this->node_[node+1], &this->node_[node+this->Nx_[d]], &this->node_[node+this->Nx_[d]+1], &this->node_[nextnode], &this->node_[nextnode+1], &this->node_[nextnode+this->Nx_[d]], &this->node_[nextnode+this->Nx_[d]+1]);
									break;
							}
							counter("cell", total_cells, cell+1);
						}
					}
				}
				z += nnx * nny * nnz;
				pz += Nx_[d] * Ny_[d] * Nz_[d];
			}
			break;
        }
        default:
            break;
	}
}





/* ---------- generate_measures ---------- */
template<class T> void grid<T>::generate_measures()
{
	debug("grid<T>::generate_measures", "generating " + std::to_string(this->cell_.size()) + " measures");
	this->measure_.resize(this->cell_.size());
	for (std::size_t i=0; i<this->measure_.size(); i++)
	{
		this->measure_[i] = (this->cell_[i])->measure();
		counter("measure", this->measure_.size(), i+1);
	}
}





/* ---------- generate_neighbours ---------- */
template<class T> void grid<T>::generate_neighbours()
{
	debug("grid<T>::generate_neighbours", "generating cell neighbours");
	switch(this->type_)
	{
        case structured: {
			this->neighbour_.resize(this->cell_.size());
            std::size_t z = 0;
			for (std::size_t d=0; d<this->Nx_.size(); d++)
			{
                std::size_t nnx = ( this->Nx_[d] < 3 ) ? 1 : this->Nx_[d] - 1;
                std::size_t nny = ( this->Ny_[d] < 3 ) ? 1 : this->Ny_[d] - 1;
                std::size_t nnz = ( this->Nz_[d] < 3 ) ? 1 : this->Nz_[d] - 1;
				for (std::size_t k=0; k<nnz; k++)
				{
					for (std::size_t j=0; j<nny; j++)
					{
						for (std::size_t i=0; i<nnx; i++)
						{
                            std::size_t cell = z + k * nnx * nny + j * nnx + i;
							std::vector<std::size_t> neigh_tmp(0);
							if ( 0 < k )
							{
								if ( 0 < j )
								{
									if ( 0 < i )
									{
										neigh_tmp.push_back(cell-nnx*nny-nnx-1);
									}
									neigh_tmp.push_back(cell-nnx*nny-nnx);
									if ( i < nnx - 1 )
									{
										neigh_tmp.push_back(cell-nnx*nny-nnx+1);
									}
								}
								if ( 0 < i )
								{
									neigh_tmp.push_back(cell-nnx*nny-1);
								}
								neigh_tmp.push_back(cell-nnx*nny);
								if ( i < nnx - 1 )
								{
									neigh_tmp.push_back(cell-nnx*nny+1);
								}
								if ( j < nny - 1 )
								{
									if ( 0 < i )
									{
										neigh_tmp.push_back(cell-nnx*nny+nnx-1);
									}
									neigh_tmp.push_back(cell-nnx*nny+nnx);
									if ( i < nnx - 1 )
									{
										neigh_tmp.push_back(cell-nnx*nny+nnx+1);
									}
								}
							}
							if ( 0 < j )
							{
								if ( 0 < i )
								{
									neigh_tmp.push_back(cell-nnx-1);
								}
								neigh_tmp.push_back(cell-nnx);
								if ( i < nnx - 1 )
								{
									neigh_tmp.push_back(cell-nnx+1);
								}
							}
							if ( 0 < i )
							{
								neigh_tmp.push_back(cell-1);
							}
							if ( i < nnx - 1 )
							{
								neigh_tmp.push_back(cell+1);
							}
							if ( j < nny - 1 )
							{
								if ( 0 < i )
								{
									neigh_tmp.push_back(cell+nnx-1);
								}
								neigh_tmp.push_back(cell+nnx);
								if ( i < nnx - 1 )
								{
									neigh_tmp.push_back(cell+nnx+1);
								}
							}
							if ( k < nnz - 1 )
							{
								if ( 0 < j )
								{
									if ( 0 < i )
									{
										neigh_tmp.push_back(cell+nnx*nny-nnx-1);
									}
									neigh_tmp.push_back(cell+nnx*nny-nnx);
									if ( i < nnx - 1 )
									{
										neigh_tmp.push_back(cell+nnx*nny-nnx+1);
									}
								}
								if ( 0 < i )
								{
									neigh_tmp.push_back(cell+nnx*nny-1);
								}
								neigh_tmp.push_back(cell+nnx*nny);
								if ( i < nnx - 1 )
								{
									neigh_tmp.push_back(cell+nnx*nny+1);
								}
								if ( j < nny - 1 )
								{
									if ( 0 < i )
									{
										neigh_tmp.push_back(cell+nnx*nny+nnx-1);
									}
									neigh_tmp.push_back(cell+nnx*nny+nnx);
									if ( i < nnx - 1 )
									{
										neigh_tmp.push_back(cell+nnx*nny+nnx+1);
									}
								}
							}
							(this->neighbour_[cell]).resize(neigh_tmp.size());
							std::copy(neigh_tmp.begin(), neigh_tmp.end(), &((this->neighbour_[cell])[0]));
						}
					}
				}
				z += nnx * nny * nnz;
			}
			break;
        }
        default:
            break;
	}
}





/* ---------- generate_surface ---------- */
template<class T> void grid<T>::generate_surface()
{
    this->surface_.clear();
    this->surface_neighbour_.clear();
    std::size_t total_surfs = 0;
    switch(this->type_)
    {
        case structured: {
            std::size_t z = 0; // surface domain counter
            for (std::size_t d=0; d<Nx_.size(); d++)
            {
                std::size_t nnx = ( this->Nx_[d] < 3 ) ? 1 : this->Nx_[d] - 1;
                std::size_t nny = ( this->Ny_[d] < 3 ) ? 1 : this->Ny_[d] - 1;
                std::size_t nnz = ( this->Nz_[d] < 3 ) ? 1 : this->Nz_[d] - 1;
                switch(this->dim_)
                {
                    case 2:
                        z += 2 * ( nnx + nny );
                        break;
                    default:
                        z += 2 * ( nnx * nny + nny * nnz + nnz * nnx );
                        break;
                }
            }
            std::string vtxt = " " + std::string(typeid(compact_shape<T>).name());
            switch(this->dim_)
            {
                case 2:
                    vtxt = " " + std::string(typeid(segment<T>).name());
                    break;
                case 3:
                    vtxt = " " + std::string(typeid(polygon<T>).name());
                    break;
                default:
                    break;
            }
            vtxt += " surface cells";
            total_surfs = z;
            debug("grid<T>::generate_surface", "calculating " + std::to_string(total_surfs) + vtxt);
            
            z = 0; // cell domain counter
            std::size_t snc = 0; // number of surface neighbours cells
            for (std::size_t d=0; d<this->Nx_.size(); d++)
            {
                std::size_t sz = this->surface_.size();
                std::size_t nnx = ( this->Nx_[d] < 3 ) ? 1 : this->Nx_[d] - 1;
                std::size_t nny = ( this->Ny_[d] < 3 ) ? 1 : this->Ny_[d] - 1;
                std::size_t nnz = ( this->Nz_[d] < 3 ) ? 1 : this->Nz_[d] - 1;
                std::valarray<std::size_t> sntmp;
                switch(this->dim_)
                {
                    case 2: {
                        // face 0
                        for (std::size_t i=z; (z<=i)&&(i<=z+nnx-1); i+=1)
                        {
                            this->surface_.push_back(grid_surface<T>(i, (this->cell_[i])->face(0)));
                            counter("surface", total_surfs, this->surface_.size());
                        }
                        // face 1
                        for (std::size_t i=z+nnx-1; (z+nnx-1<=i)&&(i<=z+nnx*nny-1); i+=nnx)
                        {
                            this->surface_.push_back(grid_surface<T>(i, (this->cell_[i])->face(1)));
                            counter("surface", total_surfs, this->surface_.size());
                        }
                        // face 2
                        for (std::size_t i=z+nnx*nny-1; (z+nnx*nny-nnx<=i)&&(i<=z+nnx*nny-1); i-=1)
                        {
                            this->surface_.push_back(grid_surface<T>(i, (this->cell_[i])->face(2)));
                            counter("surface", total_surfs, this->surface_.size());
                        }
                        // face 3
                        for (std::size_t i=z+nnx*nny-nnx; (z<=i)&&(i<=z+nnx*nny-nnx); i-=nnx)
                        {
                            this->surface_.push_back(grid_surface<T>(i, (this->cell_[i])->face(3)));
                        }
                        // surface neighbours
                        total_surfs = 4*(nnx+nny);
                        debug("grid<T>::generate_surface","calculating " + std::to_string(total_surfs) + " surface neighbours");
                        sntmp.resize(2);
                        sntmp[0] = sz + 2*(nnx+nny)-1;
                        sntmp[1] = sz + 1;
                        this->surface_neighbour_.push_back(sntmp);
                        counter("surface_neighbour", total_surfs, 2 * this->surface_neighbour_.size());
                        for (std::size_t i=1; i<2*(nnx+nny); i++)
                        {
                            sntmp[0] = sz + i-1;
                            sntmp[1] = sz + ((i+1)%(2*(nnx+nny)));
                            this->surface_neighbour_.push_back(sntmp);
                            counter("surface_neighbour", total_surfs, 2 * this->surface_neighbour_.size());
                        }
                        break;
                    }
                    case 3: {
                        std::size_t tsc = 0;
                        // face 0
                        for (std::size_t i=z; i<=z+nnx*nny-1; i+=1)
                        {
                            this->surface_.push_back(grid_surface<T>(i, (this->cell_[i])->face(0)));
                            counter("surface", total_surfs, this->surface_.size());
                        }
                        // face 1
                        for (std::size_t i=z+nnx-1; i<=z+nnx*nny*nnz-1; i+=nnx)
                        {
                            this->surface_.push_back(grid_surface<T>(i, (this->cell_[i])->face(1)));
                            counter("surface", total_surfs, this->surface_.size());
                        }
                        // face 2
                        tsc = 0;
                        for (std::size_t i=z+nnx*nny-nnx; i<=z+nnx*nny*nnz-1; i+=1)
                        {
                            this->surface_.push_back(grid_surface<T>(i, (this->cell_[i])->face(2)));
                            counter("surface", total_surfs, this->surface_.size());
                            tsc++;
                            if ( tsc == nnx )
                            {
                                i += nnx*nny-nnx;
                                tsc = 0;
                            }
                        }
                        // face 3
                        for (std::size_t i=z+nnx*nny*nnz-nnx*nny; i<=z+nnx*nny*nnz-1; i+=1)
                        {
                            this->surface_.push_back(grid_surface<T>(i, (this->cell_[i])->face(3)));
                            counter("surface", total_surfs, this->surface_.size());
                        }
                        // face 4
                        for (std::size_t i=z; i<=z+nnx*nny*nnz-nnx; i+=nnx)
                        {
                            this->surface_.push_back(grid_surface<T>(i, (this->cell_[i])->face(4)));
                            counter("surface", total_surfs, this->surface_.size());
                        }
                        // face 5
                        tsc = 0;
                        for (std::size_t i=z; i<=z+nnx*nny*nnz-nnx*nny+nnx-1; i+=1)
                        {
                            this->surface_.push_back(grid_surface<T>(i, (this->cell_[i])->face(5)));
                            counter("surface", total_surfs, this->surface_.size());
                            tsc++;
                            if ( tsc == nnx )
                            {
                                i += nnx*nny-nnx;
                                tsc = 0;
                            }
                        }
                        // surface neighbours
                        total_surfs = 16 * (nnx*nny+nny*nnz+nnz*nnx) - 24;
                        debug("grid<T>::generate_surface", "calculating " + std::to_string(total_surfs) + " surface neighbours");
                        for (std::size_t i=sz; i<this->surface_.size(); i++)
                        {
                            std::vector<std::size_t> sntmpv(0);
                            for (std::size_t j=sz; j<this->surface_.size(); j++)
                            {
                                if ( j == i )
                                {
                                    continue;
                                }
                                else if ( (this->surface_[i].cell())->adjacent(*(this->surface_[j].cell())) )
                                {
                                    sntmpv.push_back(j);
                                }
                            }
                            sntmp.resize(sntmpv.size());
                            std::copy(sntmpv.begin(), sntmpv.end(), &(sntmp[0]));
                            this->surface_neighbour_.push_back(sntmp);
                            snc += sntmp.size();
                            counter("surface_neighbour", total_surfs, snc);
                        }
                        break;
                    }
                    default:
                        break;
                }
                z += nnx * nny * nnz;
            }
            break;
        }
        default:
            break;
    }
}





/* ---------- generate_scale ---------- */
template<class T> void grid<T>::generate_scale()
{
	// calculate average scale
	this->cscale_ = T(0);
	for (std::size_t i=0; i<this->cell_.size(); i++)
	{
		this->cscale_ += (this->cell_[i])->scale();
	}
	this->cscale_ /= T((this->cell_).size());
	
	// calculate key geometric points
	this->min_.resize(this->dim_);
	min_ = node_[0];
	this->max_.resize(this->dim_);
	max_ = node_[0];
	this->centroid_.resize(this->dim_, T(0));
	for (std::size_t i=0; i<this->node_.size(); i++)
	{
		for (unsigned short j=0; j<this->dim_; j++)
		{
			if ( this->node_[i][j] < this->min_[j] )
			{
				this->min_[j] = this->node_[i][j];
			}
			if ( this->node_[i][j] > this->max_[j] )
			{
				this->max_[j] = this->node_[i][j];
			}
		}
		this->centroid_ += this->node_[i];
	}
	this->centroid_ /= T(this->node_.size());
	this->center_.resize(this->dim_);
	this->center_ = ( this->max_ + this->min_ ) / T(2);
}





/* ---------- generate_bcell ---------- */
template<class T> void grid<T>::generate_bcell()
{
	this->bnode_.resize(ipow(2, this->dim_));
	for (std::size_t i=0; i<this->bnode_.size(); i++)
	{
		this->bnode_[i].resize(this->dim_);
	}
	switch(this->dim_)
	{
		case 2:
			this->bnode_[0] = this->min_;
			this->bnode_[1] = this->min_;
			this->bnode_[1][0] = this->max_[0];
			this->bnode_[2] = this->max_;
			this->bnode_[3] = this->max_;
			this->bnode_[3][0] = this->min_[0];
			bcell_ = new polygon<T>(&this->bnode_[0], &this->bnode_[1], &this->bnode_[2], &this->bnode_[3]);
			break;
		case 3:
			this->bnode_[0] = this->min_;
			this->bnode_[1] = this->min_;
			this->bnode_[1][0] = this->max_[0];
			this->bnode_[2] = this->min_;
			this->bnode_[2][1] = this->max_[1];
			this->bnode_[3] = this->min_;
			this->bnode_[3][0] = this->max_[0];
			this->bnode_[3][1] = this->max_[1];
			this->bnode_[4] = this->max_;
			this->bnode_[4][0] = this->min_[0];
			this->bnode_[4][1] = this->min_[1];
			this->bnode_[5] = this->max_;
			this->bnode_[5][1] = this->min_[1];
			this->bnode_[6] = this->max_;
			this->bnode_[6][0] = this->min_[0];
			this->bnode_[7] = this->max_;
			bcell_ = new polyhedron<T>(&this->bnode_[0], &this->bnode_[1], &this->bnode_[2], &this->bnode_[3],
									   &this->bnode_[4], &this->bnode_[5], &this->bnode_[6], &this->bnode_[7]);
			break;
        default:
            break;
	}
}





/* ---------- clear_cells ---------- */
template<class T> void grid<T>::clear_cells()
{
	for (std::size_t i=this->cell_.size(); i>0; i--)
	{
		delete this->cell_[i-1];
	}
	this->cell_.clear();
}





/* ---------- grid ---------- */
template<class T> inline grid<T>::grid(const grid_input& thisgridinput, const bool& auxs)
{
	this->init(thisgridinput, auxs);
}





/* ---------- init ---------- */
template<class T> void grid<T>::init(const grid_input& thisgrid, const bool& auxs)
{
	debug("grid<T>::init", "initialising grid");
	if ( thisgrid.load_grid() && thisgrid.gridfile() == std::string("") )
	{
		throw;
        return;
	}
	
	// clear any existing grid
	this->clear();
	
	// load the grid
	if ( thisgrid.load_grid() )
	{
		this->read(thisgrid);
	}
	else
	{
		this->generate_grid(thisgrid);
	}
	
	if ( !this->dim_ )
	{
		this->dim_ = this->node_[0].size();
	}
	
	// check the grid
	if ( !this->check() )
	{
		this->clear();
		throw;
        return;
	}
	
	// create the grid
	this->generate_cells();
	
	if ( auxs )
	{
		this->auxs();
	}
	
	this->basis_ = unity;
}





/* ---------- auxs ---------- */
template<class T> void grid<T>::auxs()
{
	this->generate_measures();
	this->aux();
	
	// create characteristic scale length
	this->generate_scale();
	
	// creating bounding cell around entire grid
	this->generate_bcell();
}





/* ---------- ~grid ---------- */
template<class T> inline grid<T>::~grid()
{
	this->clear();
}





/* ---------- clear ---------- */
template<class T> inline void grid<T>::clear()
{
	// remove the stack
	this->clear_aux();
	this->clear_cells();
	delete this->bcell_;
	this->bnode_.clear();
	
	// zero all elements
	this->type_ = empty;
	this->node_.clear();
	this->Nx_.resize(0);
	this->Ny_.resize(0);
	this->Nz_.resize(0);
	this->data_.resize(0);
	this->adata_.resize(0);
	this->aux_ = false;
	this->cell_walk_ = 0;
	this->surface_walk_ = 0;
	this->cscale_ = T(0);
	this->min_.resize(0);
	this->max_.resize(0);
	this->centroid_.resize(0);
	this->basis_ = unity;
	
	this->dataname_.clear();
	this->gfilename_.clear();
	this->format_ = Formatted;
	this->precision_ = Single;
	this->byte_swapping_ = false;
}





/* ---------- read ---------- */
template<class T> void grid<T>::read(const grid_input& thisgrid)
{
	debug("grid<T>::read", "loading grid");
	std::valarray<T> X, Y, Z;
	std::valarray<bool> B;
	bool bx = false;
	bool by = false;
	bool bz = false;
	std::valarray<T>* Xp = &X;
	std::valarray<T>* Yp = &Y;
	switch(thisgrid.type())
	{
		case structured:
			this->type_ = thisgrid.type();
			this->format_ = thisgrid.format();
			this->precision_ = thisgrid.precision();
			this->gfilename_ = thisgrid.gridfile();
			Plot3D::read(X, Y, Z, B, this->Nx_, this->Ny_, this->Nz_, this->gfilename_, this->format_, this->precision_, thisgrid.multidomain(), thisgrid.blanking());
			for (std::size_t i=0; i<X.size()-1; i++)
			{
				if ( X[i] != X[i+1] )
				{
					bx = true;
					break;
				}
			}
			for (std::size_t j=0; j<Y.size()-1; j++)
			{
				if ( Y[j] != Y[j+1] )
				{
					by = true;
					break;
				}
			}
			for (std::size_t k=0; k<Z.size()-1; k++)
			{
				if ( Z[k] != Z[k+1] )
				{
					bz = true;
					break;
				}
			}
			if ( bx & !by & !bz )
			{
				this->dim_ = 1;
				Xp = &X;
				Y.resize(0);
				Z.resize(0);
			}
			else if ( !bx & by & !bz )
			{
				this->dim_ = 1;
				Xp = &Y;
				Z.resize(0);
				X.resize(0);
			}
			else if ( !bx & !by & bz )
			{
				this->dim_ = 1;
				Xp = &Z;
				X.resize(0);
				Y.resize(0);
			}
			else if ( bx & by & !bz )
			{
				this->dim_ = 2;
				Xp = &X;
				Yp = &Y;
				Z.resize(0);
			}
			else if ( !bx & by & bz )
			{
				this->dim_ = 2;
				Xp = &Y;
				Yp = &Z;
				X.resize(0);
			}
			else if ( bx & !by & bz )
			{
				this->dim_ = 2;
				Xp = &Z;
				Yp = &X;
				Y.resize(0);
			}
			else
			{
				this->dim_ = 3;
			}
			this->node_.resize((*Xp).size(), std::valarray<T>(this->dim_));
			if ( this->dim_ == 1 )
			{
				for (std::size_t loop=0; loop<this->node_.size(); loop++)
				{
					this->node_[loop][0] = (*Xp)[loop];
				}
			}
			else if ( this->dim_ == 2 )
			{
				for (std::size_t loop=0; loop<this->node_.size(); loop++)
				{
					this->node_[loop][0] = (*Xp)[loop];
					this->node_[loop][1] = (*Yp)[loop];
				}
			}
			else
			{
				for (std::size_t loop=0; loop<this->node_.size(); loop++)
				{
					this->node_[loop][0] = X[loop];
					this->node_[loop][1] = Y[loop];
					this->node_[loop][2] = Z[loop];
				}
			}
			break;
        default:
            break;
	}
	debug("grid<T>::read", "loaded " + std::to_string(this->dim_) + " dimensional grid");
}





/* ---------- write ---------- */
template<class T> void grid<T>::write(const std::string& filename,
                                      const grid_saveitem& what,
                                      const dataformat& writeformat) const
{
	std::string function = "grid<T>::write";
	
	std::string ext = ".P";
	switch(writeformat)
	{
		case Formatted:
			ext += "F";
			break;
		case Unformatted:
			ext += "U";
			break;
		case Binary:
			ext += "B";
			break;
        default:
            break;
	}
	
	dataprecision fsize = Single;
	std::string pext = "";
	//if ( sizeof(T) > 4 )
	//{
	//  fsize = Double;
	//  pext = "D";
	//}
	
	std::valarray<T> X(this->node_.size()), Y(this->node_.size()), Z(this->node_.size());
	for (std::size_t i=0; i<this->node_.size(); i++)
	{
		X[i] = this->node_[i][0];
		Y[i] = this->node_[i][1];
		if ( this->dim_ > 2 )
		{
			Z[i] = this->node_[i][2];
		}
		else
		{
			Z[i] = T(0);
		}
	}
	std::valarray<bool> B(true, this->node_.size());
	
	std::valarray<T> output(0);
	
	std::valarray<T> nothing;
	
	switch(what)
	{
		case SaveGrid:
			debug(function,"saving grid `"+filename+ext+"G"+pext+"'");
			switch(this->type_)
			{
				case structured:
					Plot3D::write(X, Y, Z, B, this->Nx_, this->Ny_, this->Nz_, filename+ext+"G"+pext, writeformat, fsize, false, false);
					break;
                default:
                    break;
			}
			break;
		case SaveData:
			write_data(filename, writeformat);
			break;
		case SaveNeighbours:
			debug(function,"saving grid cell neighbours `"+filename+ext+"N"+"'");
			switch(this->type_)
			{
				case structured:
					Plot3D::write_neighbours(this->node_, this->Nx_, this->Ny_, this->Nz_, this->neighbour_, filename+ext+"N", writeformat);
					break;
                default:
                    break;
			}
			break;
        default:
            break;
	}
}





/* ---------- read_data ---------- */
template<class T> void grid<T>::read_data(const grid_input& thisgrid)
{
	debug("grid<T>::read_data", "reading grid data");
	if ( thisgrid.datafile() == std::string("") )
	{
		debug("grid<T>::read_data", "no data file given");
		throw; return;
	}
	switch(this->type_)
	{
        case structured: {
			std::valarray<T> p3d_data;
			Plot3D::extract_data(p3d_data, thisgrid.datafile(), this->format_, this->precision_, thisgrid.multidomain(), thisgrid.qdata());
			this->adata_.resize(this->cell_.size());
			Plot3D::n2c(this->node_, this->Nx_, this->Ny_, this->Nz_, p3d_data, this->adata_);
			break;
        }
        default:
            break;
	}
}





/* ---------- put_adata ---------- */
template<class T> inline void grid<T>::put_adata(const std::valarray<T>& user_adata)
{
	this->adata_.resize(this->cell_.size());
	this->adata_ = user_adata;
}





/* ---------- get_adata ---------- */
template<class T> inline std::valarray<T> grid<T>::get_adata() const
{
	return this->adata_;
}





/* ---------- write_data ---------- */
template<class T> void grid<T>::write_data(const std::string& filename,
                                           const dataformat& writeformat,
                                           const std::valarray<T>& extradata) const
{
	std::string function = "grid<T>::write_data";
	
	std::string ext = ".P";
	switch(writeformat)
	{
		case Formatted:
			ext += "FS";
			break;
		case Unformatted:
			ext += "US";
			break;
		case Binary:
			ext += "BS";
			break;
        default:
            break;
	}
	
	dataprecision fsize = Single;
	//if ( sizeof(T) > 4 )
	//{
	//  fsize = Double;
	//  ext += "D";
	//}
	
	std::valarray<T> output1, output2, output3;
	std::string dataname1 = "ZERO_1", dataname2 = "ZERO_2", dataname3 = "ZERO_3";
	
	std::valarray<T> nothing;
	
	debug(function,"saving grid data `" + filename + ext + "'");
	
	// data array 1
	if ( this->data_.size() )
	{
		switch(this->type_)
		{
			case structured:
				Plot3D::c2n(this->node_, this->Nx_, this->Ny_, this->Nz_, this->data_, output1);
				break;
            default:
                break;
		}
		dataname1 = this->dataname_;
		debug(dataname1);
	}
	
	// data array 2
	if ( this->adata_.size() )
	{
		if ( output1.size() )
		{
			switch(this->type_)
			{
				case structured:
					Plot3D::c2n(this->node_, this->Nx_, this->Ny_, this->Nz_, this->adata_, output2);
					break;
                default:
                    break;
			}
			dataname2 = this->dataname_ + "_apriori";
			debug(dataname2);
		}
		else
		{
			switch(this->type_)
			{
				case structured:
					Plot3D::c2n(this->node_, this->Nx_, this->Ny_, this->Nz_, this->adata_, output1);
					break;
                default:
                    break;
			}
			dataname1 = this->dataname_+"_apriori";
			debug(dataname1);
		}
	}
	
	// data array 3
	if ( extradata.size() )
	{
		if ( output2.size() )
		{
			switch(this->type_)
			{
				case structured:
					Plot3D::c2n(this->node_, this->Nx_, this->Ny_, this->Nz_, extradata, output3);
					break;
                default:
                    break;
			}
			dataname3 = "extra_data";
			debug(dataname3);
		}
		else if ( output1.size() )
		{
			switch(this->type_)
			{
				case structured:
					Plot3D::c2n(this->node_, this->Nx_, this->Ny_, this->Nz_, extradata, output2);
					break;
                default:
                    break;
			}
			dataname2 = "extra_data";
			debug(dataname2);
		}
		else
		{
			switch(this->type_)
			{
				case structured:
					Plot3D::c2n(this->node_, this->Nx_, this->Ny_, this->Nz_, extradata, output1);
					break;
                default:
                    break;
			}
			dataname1 = "extra_data";
			debug(dataname1);
		}
	}
	
	// no data
	if ( !output1.size() )
	{
		debug(function, "no data to save");
	}
	
	switch(this->type_)
	{
		case structured:
			Plot3D::write_data(output1, output2, output3, nothing, nothing, this->Nx_, this->Ny_, this->Nz_, filename+ext, writeformat, fsize, false);
			Plot3D::write_var(filename + ".VAR", dataname1, dataname2, dataname3, "ZERO_4", "ZERO_5");
			break;
        default:
            break;
	}
}





/* ---------- check ---------- */
template<class T> bool grid<T>::check() const
{
	std::string fname = "grid<T>::check";
	switch(this->type_)
	{
        case structured: {
			fname += " structured";
			// check array sizes
			if ( this->Nx_.size() != this->Ny_.size() || this->Ny_.size() != this->Nz_.size() ||
				this->Nz_.size() != this->Nx_.size() || this->node_.size() < 2 || this->node_[0].size() < 1 )
			{
				return false;
			}
			// check grid properties
            std::size_t numnodes = 0;
			for (std::size_t i=0; i<this->Nx_.size(); i++)
			{
				if ( this->Nx_[i] < 1 || this->Ny_[i] < 1 || this->Nz_[i] < 1 )
				{
					return false;
				}
				numnodes += this->Nx_[i] * this->Ny_[i] * this->Nz_[i];
			}
			if ( numnodes != this->node_.size() )
			{
				return false;
			}
			// check node dimensions
			for (std::size_t i=0; i<this->node_.size(); i++)
			{
				if ( this->node_[i].size() != this->dim_ )
				{
					return false;
				}
			}
			return true;
            break;
        }
        default:
            break;
	}
	return false;
}





/* ---------- aux ---------- */
template<class T> void grid<T>::aux()
{
	this->generate_neighbours();
	this->generate_surface();
	this->aux_ = true;
}





/* ---------- clear_aux ---------- */
template<class T> void grid<T>::clear_aux()
{
	this->neighbour_.clear();
	this->surface_.clear();
	this->surface_neighbour_.clear();
	this->aux_ = false;
}





/* ---------- dim ---------- */
template<class T> inline unsigned short grid<T>::dim() const
{
	return this->dim_;
}





/* ---------- ncells ---------- */
template<class T> inline std::size_t grid<T>::ncells() const
{
	return this->cell_.size();
}





/* ---------- center ---------- */
/* finds the center of the grid */
template<class T> inline std::valarray<T> grid<T>::center() const
{
	return this->center_;
}





/* ---------- min ---------- */
//// finds the minimum bound for the grid
template<class T> inline std::valarray<T> grid<T>::min() const
{
	return this->min_;
}





/* ---------- max ---------- */
/// finds the maximum bound for the grid
template<class T> inline std::valarray<T> grid<T>::max() const
{
	return this->max_;
}





/* ---------- raytrace ---------- */
/// finds the projection of a line/ray through a grid
template<class T> inline T grid<T>::raytrace(const line<T>& rayline,
                                             std::vector< two_numbers<T> >& raydata,
                                             const grid_proj_method& method)
{
	raydata.clear();
	T projected_sum = T(0);
	bool parallel = false;
	
	if ( !this->bcell_->intersect(rayline, parallel) )
	{
		return projected_sum;
	}
	
	std::valarray<T> cuts;
	T raylen = T(0);
	
	switch(method)
	{
		case brute:
			// loop over ALL cells in the grid (hence brute)
			for (std::size_t i=0; i<this->cell_.size(); i++)
			{
				raylen = this->cell_[i]->intersect(rayline, cuts, parallel);
				if ( raylen && cuts.size() )
				{
					raydata.push_back(two_numbers<T>(i, this->applybasis(i, raylen)));
					if ( this->adata_.size() )
					{
						projected_sum += ( parallel ) ? raylen * this->adata_[i] / T(2) : raylen * this->adata_[i];
					}
					else
					{
						projected_sum += ( parallel ) ? raylen / T(2) : raylen;
					}
				}
			}
			break;
			
		case walk: case walkfast:
			// find a face which the line passes through, so we have a
			// starting cell that is on the boundary.
			// we also need to keep a copy so that we make sure that
			// the ray tracing does not miss any non-convex regions
			
			// next in line to be checked
			std::vector<std::size_t> nil(0);
			
			// check local surface neighbours, for faster ray-tracing
			if ( method == walkfast )
			{
				if ( (this->surface_[this->surface_walk_]).cell()->intersect(rayline, parallel) )
				{
					nil.push_back((this->surface_[this->surface_walk_]).cellnumber());
				}
				else
				{
					for (std::size_t i=0; i<this->surface_neighbour_[this->surface_walk_].size(); i++)
					{
                        std::size_t tsn = this->surface_neighbour_[this->surface_walk_][i];
						if ( (this->surface_[tsn]).cell()->intersect(rayline, parallel) )
						{
							nil.push_back((this->surface_[tsn]).cellnumber());
							this->surface_walk_ = tsn;
							break;
						}
					}
				}
			}
			
			// check all surface cells (in the case that the quicksearch failed)
			if ( !nil.size() )
			{
				for (std::size_t i=0; i<this->surface_.size(); i++)
				{
					if ( (this->surface_[i]).cell()->intersect(rayline, cuts, parallel) )
					{
						nil.push_back((this->surface_[i]).cellnumber());
						this->surface_walk_ = i;
						if ( method == walkfast )
						{
							break;
						}
					}
				}
			}
			
			if ( !nil.size() )
			{
				return T(0);
			}
			remove_duplicates(nil);
			
			// set up an array of flags to store information about
			// which cells have already been dealt with
			std::valarray<bool> checked(false, this->cell_.size());
			for (std::size_t i=0; i<nil.size(); i++)
			{
				checked[nil[i]] = true;
			}
			
			// Now for the actual ray tracing.
			// We scan only the natural neighbours, and follow EVERY
			// single path. This will be okay since most paths will
			// end on that neighbour, and there should only be one
			// path which continues along the ray.
			
			while ( nil.size() )
			{
				std::vector<std::size_t> hits(0);
				for (std::size_t i=0; i<nil.size(); i++)
				{
					raylen = (this->cell_[nil[i]])->intersect(rayline, cuts, parallel);
					if ( cuts.size() )
					{
						this->cell_walk_ = nil[i];
						hits.push_back(nil[i]);
						if ( raylen )
						{
							raydata.push_back(two_numbers<T>(nil[i],this->applybasis(nil[i],raylen)));
							if ( this->adata_.size() )
							{
								projected_sum += ( parallel ) ? raylen * this->adata_[nil[i]] / T(2) : raylen * this->adata_[nil[i]];
							}
							else
							{
								projected_sum += ( parallel ) ? raylen / T(2): raylen;
							}
						}
					}
				}
				nil.clear();
				for (std::size_t i=0; i<hits.size(); i++)
				{
					for (std::size_t j=0; j<this->neighbour_[hits[i]].size(); j++)
					{
						if ( !checked[this->neighbour_[hits[i]][j]] )
						{
							nil.push_back(this->neighbour_[hits[i]][j]);
							checked[this->neighbour_[hits[i]][j]] = true;
						}
					}
				}
			}
			break;
	}
	
	if ( !raydata.size() )
	{
		return T(0);
	}
	if ( method != brute )
	{
		std::sort(raydata.begin(), raydata.end());
	}
	
	return projected_sum;
}





/* ---------- shift ---------- */
/// shift the entire grid by the ammount defined in the arguement vector
template<class T> void grid<T>::shift(const std::valarray<T>& shift_ammount)
{
	// check dimensions
	if ( this->dim_ != shift_ammount.size() )
	{
		throw; return;
	}
	
	bool thisaux = this->aux_;
	if ( thisaux )
	{
		this->clear_aux();
	}
	this->clear_cells();
	std::for_each(this->node_.begin(), this->node_.end(), AddVector<T>(shift_ammount));
	this->center_ += shift_ammount;
	this->centroid_ += shift_ammount;
	this->min_ += shift_ammount;
	this->max_ += shift_ammount;
	this->generate_cells();
	if ( thisaux )
	{
		this->aux();
	}
}





/* ---------- give_dataname ---------- */
template<class T> inline void grid<T>::give_dataname(const std::string& nm)
{
	this->dataname_ = nm;
}





/* ---------- get_dataname ---------- */
template<class T> inline std::string grid<T>::get_dataname() const
{
	return this->dataname_;
}





/* ---------- get_filename ---------- */
template<class T> inline std::string grid<T>::get_filename() const
{
	std::string::size_type start = this->gfilename_.find_last_of(std::string("/")) + 1;
	return this->gfilename_.substr(start, this->gfilename_.size() - start);
}





/* ---------- operator[] ---------- */
template<class T> inline T& grid<T>::operator[] (const std::size_t& tcell)
{
	if ( this->data_.size() != this->cell_.size() )
	{
		this->data_.resize(this->cell_.size(), T(0));
	}
	return this->data_[tcell];
}





/* ---------- operator[] ---------- */
template<class T> inline T grid<T>::operator[] (const std::size_t& tcell) const
{
	if ( this->data_.size() != this->cell_.size() )
	{
		this->data_.resize(this->cell_.size(), T(0));
	}
	return this->data_[tcell];
}





/* ---------- clear_data ---------- */
template<class T> inline void grid<T>::clear_data()
{
	this->data_.resize(0);
}





/* ---------- set_basis ---------- */
template<class T> inline void grid<T>::set_basis(const grid_basis& this_basis)
{
	this->basis_ = this_basis;
}





/* ---------- applybasis ---------- */
template<class T> inline T grid<T>::applybasis(const std::size_t& thiscell,
                                               const T& value) const
{
	switch(this->basis_)
	{
		case l2norm:
			return value / std::sqrt((this->cell_[thiscell])->measure());
			break;
        default:
            break;
	}
	return value;
}





/* ---------- scale ---------- */
template<class T> inline T grid<T>::scale() const
{
	return this->cscale_;
}





/* ---------- get_neighbours ---------- */
template<class T> inline std::valarray<std::size_t> grid<T>::get_neighbours(const std::size_t& cell) const
{
	return this->neighbour_[cell];
}





#endif /* _GRID_ */
