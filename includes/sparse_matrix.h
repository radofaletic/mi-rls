/**
 sparse_matrix
 
 sparse matrix class
 
 Rado Faletic
 Department of Physics
 Faculty of Science
 Australian National University  ACT  0200
 Australia
 
 Rado.Faletic@anu.edu.au
 17th April 2005
 22nd April 2022, updated to C++20
 */





#ifndef _SPARSE_MATRIX_
#define _SPARSE_MATRIX_





/* ---------- standard header files ---------- */
#include <algorithm>
#include <bit>
#include <fstream>
#include <ios>
#include <cmath>
#include <sstream>
#include <string>
#include <valarray>
#include <vector>





/* ---------- user header files ---------- */
#include "conversions.h"
#include "extra_math.h"
#include "file.h"
#include "fortran_io.h"
#include "front-end.h"
#include "matrix.h"





/* ---------- class & function declaration ---------- */

template<class T> class SparseMatrix : public Matrix<T>
{
private:
	// std::size_t rows_;                    // from Matrix<T> class
	// std::size_t cols_;                    // from Matrix<T> class
	// std::valarray<T> data_;               // from Matrix<T> class
	std::vector<T> vdata_;                // data in vector format
	std::valarray<std::size_t> elements_; // positions of non-zero elements
	std::vector<std::size_t> velements_;  // positions of non-zero elements
	std::valarray<std::size_t> ne_;       // number of non-zero elements in each row
	std::vector<std::size_t> vne_;        // number of non-zero elements in each row
	std::valarray<bool> referenced_;      // columns numbers referenced
	// NOTE: With this SparseMatrix we want to do two incompatible operations. Mainly, dynamically add to the SparseMatrix, in it's construction phase, and maintain a highly efficient numerical algorithm. The solution presented here is to use std::vector in the construction phase, to enable fast and easy dynamic changes to the size and contents of the SparseMatrix. Then, we switch to std::valarray once the SparseMatrix is ready for calculations. The speed efficiencies in doing so are considerable, compared with either case separately.
	bool fixed_;
	void fix();
	void unfix();
public:
	SparseMatrix(const std::vector< std::valarray<std::size_t> >&, const std::vector< std::valarray<T> >&);
	SparseMatrix(const std::size_t&, const std::size_t&, const T& = T(1));
	void clear();
	~SparseMatrix() { };
	std::string print(const bool& = false);
	bool write(const std::string&, const dataformat& = Unformatted);
	bool write(const std::string&, const std::valarray<T>&, const dataformat& = Unformatted);
	bool read(const std::string&, const dataformat& = Unformatted);
	bool read(const std::string&, std::valarray<T>&, const dataformat& = Unformatted);
	void AddRow(const std::valarray<std::size_t>&, const std::valarray<T>&);
	void AddRow(const std::vector< two_numbers<T> >&);
	const std::valarray<bool>* Referenced(const bool& = false);
	void AddNonReferenced(std::valarray<T>&, const T& = T(0));
	const std::valarray<bool>* Compress(std::valarray<T>&);
	const std::valarray<bool>* Uncompress(std::valarray<T>&, const T& = T(0));
	bool AddDiagonal(const std::size_t& = 0, const T& = T(1));
	bool Add1Reverse(const std::valarray<std::size_t>&, const std::size_t&, const T& = T(1));
	bool AddNReverse(const std::vector< std::valarray<std::size_t> >&, const std::vector<std::size_t>&, const T& = T(1));
	T operator()(const std::size_t&, const std::size_t&);
	std::valarray<T> operator*(const std::valarray<T>&);
	void multiply(std::valarray<T>&, const std::valarray<T>&);
	void multiplyT(std::valarray<T>&, const std::valarray<T>&);
	T norm() const;
    std::size_t non_zeros() const;
	// std::valarray<float> x = A.iterate(u); // solves x from A * x = u
	std::valarray<T> iterate(std::valarray<T>&, const std::size_t& = 100, const std::size_t& = 10);
};





/* ---------- function definitions ---------- */





/* ---------- fix() ---------- */
template<class T> inline void SparseMatrix<T>::fix()
{
	if ( !this->fixed_ )
	{
		this->data_.resize(this->vdata_.size());
		std::copy(this->vdata_.begin(), this->vdata_.end(), &(this->data_[0]));
		this->vdata_.clear();
		
		this->elements_.resize(this->velements_.size());
		std::copy(this->velements_.begin(), this->velements_.end(), &(this->elements_[0]));
		this->velements_.clear();
		
		this->ne_.resize(this->vne_.size());
		std::copy(this->vne_.begin(), this->vne_.end(), &(this->ne_[0]));
		this->vne_.clear();
	}
	this->fixed_ = true;
}





/* ---------- unfix() ---------- */
template<class T> inline void SparseMatrix<T>::unfix()
{
	if ( this->fixed_ )
	{
		this->vdata_.resize(this->data_.size());
		std::copy(&(this->data_[0]), &(this->data_[this->data_.size()]), this->vdata_.begin());
		this->data_.resize(0);
		
		this->velements_.resize(this->elements_.size());
		std::copy(&(this->elements_[0]), &(this->elements_[this->elements_.size()]), this->velements_.begin());
		this->elements_.resize(0);
		
		this->vne_.resize(this->ne_.size());
		std::copy(&(this->ne_[0]), &(this->ne_[this->ne_.size()]), this->vne_.begin());
		this->ne_.resize(0);
	}
	this->fixed_ = false;
}





/* ---------- SparseMatrix ---------- */
template<class T> SparseMatrix<T>::SparseMatrix(const std::vector< std::valarray<std::size_t> >& elements,
                                                const std::vector< std::valarray<T> >& data)
{
	this->fixed_ = false;
	if ( !data.size() || elements.size() != data.size() )
	{
		debug("SparseMatrix<T>::SparseMatrix", "inconsistent data");
		this->clear();
		throw;
        return;
	}
	
	this->rows_ = data.size();
	this->vne_.resize(this->rows_);
	
    std::size_t nze = 0; // number of total non-zero elements
	for (std::size_t i=0; i<this->rows_; i++)
	{
		if ( elements[i].size() != data[i].size() )
		{
			debug("SparseMatrix<T>::SparseMatrix", "inconsistent data");
			this->clear();
			throw;
            return;
		}
		this->vne_[i] = data[i].size();
		nze += this->vne_[i];
	}
	this->velements_.resize(nze);
	this->vdata_.resize(nze);
	nze = 0; // re-used dummy variable, acts as a counter
	for (std::size_t i=0; i<this->rows_; i++)
	{
		this->cols_ = std::max(this->cols_, elements[i].max() + 1);
		std::copy(&(elements[i][0]), &(elements[i][elements[i].size()]), &(this->velements_[nze]));
		std::copy(&(data[i][0]), &(data[i][data[i].size()]), &(this->vdata_[nze]));
		nze += this->vne_[i];
	}
	
	this->referenced_.resize(0);
}





/* ---------- SparseMatrix ---------- */
template<class T> SparseMatrix<T>::SparseMatrix(const std::size_t& rows, const std::size_t& cols, const T& init_value)
{
	this->fixed_ = false;
	if ( !rows && !cols )
	{
		this->clear();
		return;
	}
	this->rows_ = rows;
	this->cols_ = cols;
	this->vne_.resize(this->rows_, 1);
	if ( this->cols_ < this->rows_ )
	{
		std::fill(this->vne_.begin()+this->cols_-1, this->vne_.end(), 0);
		this->velements_.resize(this->cols_);
		this->vdata_.resize(this->cols_, init_value);
	}
	else
	{
		this->velements_.resize(this->rows_);
		this->vdata_.resize(this->rows_, init_value);
	}
	for (std::size_t i=0; i<this->velements_.size(); i++)
	{
		this->velements_[i] = i;
	}
	
	this->referenced_.resize(0);
}





/* ---------- clear() ---------- */
template<class T> inline void SparseMatrix<T>::clear()
{
	Matrix<T>::clear();
	this->vne_.clear();
	this->velements_.clear();
	this->vdata_.clear();
	this->ne_.resize(0);
	this->elements_.resize(0);
	this->fixed_ = false;
	this->referenced_.resize(0);
}





/* ---------- print() ---------- */
template<class T> std::string SparseMatrix<T>::print(const bool& CSV)
{
	this->unfix();
	
	std::ostringstream tmp;
	
    std::size_t rposition = 0;
	for (std::size_t i=0; i<this->rows_; i++)
	{
        std::size_t cposition = 0;
		for (std::size_t j=0; j<this->cols_; j++)
		{
			if ( j == this->velements_[rposition+cposition] && cposition < this->vne_[i] )
			{
				tmp << this->vdata_[rposition+cposition];
				cposition++;
			}
			else
			{
				tmp << "0";
			}
			if ( j < this->cols_ - 1 )
			{
				tmp << ( CSV ? "," : "\t" );
			}
			else if ( i < this->rows_ - 1 )
			{
				tmp << "\n";
			}
		}
		rposition += this->vne_[i];
	}
	return tmp.str();
}





/* ---------- write ---------- */
template<class T> inline bool SparseMatrix<T>::write(const std::string& filename, const dataformat& format)
{
	std::valarray<T> b(0);
	return this->write(filename, b, format);
}





/* ---------- write ---------- */
template<class T> bool SparseMatrix<T>::write(const std::string& filename, const std::valarray<T>& b, const dataformat& format)
{
	this->fix();
	
	std::ofstream matrixfile;
	switch(format)
	{
		case Unformatted: case Binary:
			matrixfile.open(filename.c_str(),std::ios_base::binary);
			break;
		default:
			matrixfile.open(filename.c_str());
			break;
	}
	
	bool byte_swapping = (std::endian::native == std::endian::little) ? true : false;
	
	// write the number of rows, columns and total number of non-zero elements
	int hb = ( b.size() ) ? 16 : 12;
	Fortran::write_byte_info(matrixfile, hb, byte_swapping, format);
	Fortran::iwrite_item(matrixfile, this->rows_, byte_swapping, format);
	Fortran::iwrite_item(matrixfile, this->cols_, byte_swapping, format);
	Fortran::iwrite_item(matrixfile, this->elements_.size(), byte_swapping, format);
	if ( b.size() )
	{
		Fortran::iwrite_item(matrixfile, b.size(), byte_swapping, format);
	}
	Fortran::write_byte_info(matrixfile, hb, byte_swapping, format);
	
	// write the elements and data
	Fortran::iwrite_vector(matrixfile, this->ne_, byte_swapping, format);
	Fortran::iwrite_vector(matrixfile, this->elements_, byte_swapping, format);
	if ( Fortran::is_integer(T(0)) )
	{
		Fortran::iwrite_vector(matrixfile, this->data_, byte_swapping, format);
		if ( b.size() )
		{
			Fortran::iwrite_vector(matrixfile, b, byte_swapping, format);
		}
	}
	else if ( Fortran::is_real(T(0)) )
	{
		dataprecision precision = ( sizeof(T) > 4 ) ? Double : Single;
		Fortran::fwrite_vector(matrixfile, this->data_, byte_swapping, format, precision);
		if ( b.size() )
		{
			Fortran::fwrite_vector(matrixfile, b, byte_swapping, format, precision);
		}
	}
	else
	{
		return false;
	}
	
	matrixfile.close();
	
	return true;
}





/* ---------- read ---------- */
template<class T> inline bool SparseMatrix<T>::read(const std::string& filename, const dataformat& format)
{
	std::valarray<T> b(0);
	return this->read(filename, b, format);
}





/* ---------- read ---------- */
template<class T> bool SparseMatrix<T>::read(const std::string& filename, std::valarray<T>& b, const dataformat& format)
{
	// firstly remove any "residual" matrix
	this->clear();
	b.resize(0);
	
	// open the file
	std::ifstream matrixfile;
	switch(format)
	{
		case Unformatted: case Binary:
			matrixfile.open(filename.c_str(),std::ios_base::binary);
			break;
		default:
			matrixfile.open(filename.c_str());
			break;
	}
	
	// read the number of rows, columns and total number of non-zero elements
	bool byte_swapping = (std::endian::native == std::endian::little) ? true : false;
	dataprecision precision = Single;
	int hb = 12;
	if ( !Fortran::fcheckset(matrixfile, hb, byte_swapping, format, precision) )
	{
		hb = 16;
	}
	if ( !Fortran::fcheckset(matrixfile, hb, byte_swapping, format, precision) )
	{
		debug("SparseMatrix<T>::read", "bad checkset");
		matrixfile.close();
		this->clear();
		throw;
        return false;
	}
    std::size_t isize;
	Fortran::read_byte_info(matrixfile, isize, byte_swapping, format);
	Fortran::iread_item(matrixfile, this->rows_, byte_swapping, format);
	Fortran::iread_item(matrixfile, this->cols_, byte_swapping, format);
    std::size_t nonzero;
	Fortran::iread_item(matrixfile, nonzero, byte_swapping, format);
	int bsize = 0;
	if ( hb == 16 )
	{
		Fortran::iread_item(matrixfile, bsize, byte_swapping, format);
	}
    std::size_t itemp;
	Fortran::read_byte_info(matrixfile, itemp, byte_swapping, format);
	if ( itemp != isize )
	{
		debug("SparseMatrix<T>::read", "inconsistent matrix file '"+filename+"'");
		matrixfile.close();
		this->clear();
		throw;
        return false;
	}
	
	// read the elements
	Fortran::iread_vector(matrixfile, this->rows_, this->vne_, byte_swapping, format);
	Fortran::iread_vector(matrixfile, nonzero, this->velements_, byte_swapping, format);
	
	// read the data (automatic double precision)
    std::size_t fsize = ( precision == Double ) ? 8 : 4;
	fsize *= nonzero;
	if ( !Fortran::fcheckset(matrixfile, fsize, byte_swapping, format, precision) )
	{
		debug("SparseMatrix<T>::read", "bad checkset");
		matrixfile.close();
		this->clear();
		throw;
        return false;
	}
	if ( Fortran::is_integer(T(0)) )
	{
		Fortran::iread_vector(matrixfile, nonzero, this->vdata_, byte_swapping, format);
		if ( bsize )
		{
			Fortran::iread_vector(matrixfile, bsize, b, byte_swapping, format);
		}
	}
	else if ( Fortran::is_real(T(0)) )
	{
		Fortran::fread_vector(matrixfile, nonzero, this->vdata_, byte_swapping, format, precision);
		if ( bsize )
		{
			Fortran::fread_vector(matrixfile, bsize, b, byte_swapping, format, precision);
		}
	}
	else
	{
		debug("SparseMatrix<T>::read", "can't determine numerics");
		matrixfile.close();
		this->clear();
		throw;
        return false;
	}
	
	matrixfile.close();
	
	return true;
}





/* ---------- AddRow ---------- */
template<class T> void SparseMatrix<T>::AddRow(const std::valarray<std::size_t>& elements,
                                               const std::valarray<T>& data)
{
	this->unfix();
	
	if ( !elements.size() || elements.size() != data.size() )
	{
		this->vne_.push_back(0);
		this->rows_++;
		return;
	}
	this->velements_.insert(this->velements_.end(), &(elements[0]), &(elements[elements.size()]));
	this->vdata_.insert(this->vdata_.end(), &(data[0]), &(data[data.size()]));
	this->vne_.push_back(elements.size());
	this->cols_ = std::max(this->cols_, elements.max()+1);
	this->rows_++;
}





/* ---------- AddRow ---------- */
template<class T> void SparseMatrix<T>::AddRow(const std::vector< two_numbers<T> >& rowdata)
{
	this->unfix();
	
	if ( !rowdata.size() )
	{
		this->vne_.push_back(0);
		this->rows_++;
		return;
	}
	
    std::size_t os = this->velements_.size();
	
	this->velements_.resize(this->velements_.size()+rowdata.size());
	for (std::size_t i=0; i<rowdata.size(); i++)
	{
		this->velements_[os+i] = rowdata[i].n_integer;
		this->cols_ = std::max(this->cols_, rowdata[i].n_integer+1);
	}
	
	this->vdata_.resize(this->vdata_.size()+rowdata.size());
	for (std::size_t i=0; i<rowdata.size(); i++)
	{
		this->vdata_[os+i] = rowdata[i].n_real;
	}
	
	this->vne_.push_back(rowdata.size());
	
	this->rows_++;
}





/* ---------- Referenced ---------- */
template<class T> const std::valarray<bool>* SparseMatrix<T>::Referenced(const bool& generate)
{
	if ( !this->referenced_.size() || generate )
	{
		this->unfix();
		this->referenced_.resize(this->cols_, false);
		for (std::size_t i=0; i<this->velements_.size(); i++)
		{
			this->referenced_[this->velements_[i]] = true;
		}
	}
	return &(this->referenced_);
}





/* ---------- AddNonReferenced ---------- */
template<class T> void SparseMatrix<T>::AddNonReferenced(std::valarray<T>& b, const T& value)
{
	this->Referenced();
	for (std::size_t i=0; i<this->referenced_.size(); i++)
	{
		if ( !this->referenced_[i] )
		{
			std::valarray<std::size_t> e(i,1);
			std::valarray<T> d(T(1),1);
			this->AddRow(e,d);
		}
	}
	std::valarray<T> b_temp = b;
	b.resize(this->rows_, value);
	b[std::slice(0,b_temp.size(),1)] = b_temp;
}





/* ---------- Compress ---------- */
template<class T> const std::valarray<bool>* SparseMatrix<T>::Compress(std::valarray<T>& orig_x)
{
	this->Referenced();
	this->unfix();
	std::vector<T> x(orig_x.size());
	std::copy(&(orig_x[0]), &(orig_x[orig_x.size()]), x.begin());
	
	for (std::size_t i=this->referenced_.size(); 1<=i; i--)
	{
		counter("compression", this->referenced_.size(), this->referenced_.size()-i);
		if ( !this->referenced_[i -1] )
		{
			std::for_each(this->velements_.begin(), this->velements_.end(), decreaser<std::size_t>(i -1));
			x.erase(x.begin()+(i -1));
			this->cols_--;
		}
	}
	
	orig_x.resize(x.size());
	std::copy(x.begin(), x.end(), &(orig_x[0]));
	
	return &(this->referenced_);
}





/* ---------- Uncompress ---------- */
template<class T> const std::valarray<bool>* SparseMatrix<T>::Uncompress(std::valarray<T>& orig_x, const T& filler)
{
	this->unfix();
	std::vector<T> x(orig_x.size());
	std::copy(&(orig_x[0]), &(orig_x[orig_x.size()]), x.begin());
	
	for (std::size_t i=0; i<this->referenced_.size(); i++)
	{
		counter("decompression", this->referenced_.size(), i);
		if ( !this->referenced_[i] )
		{
			std::for_each(this->velements_.begin(), this->velements_.end(), increaser<std::size_t>(i));
			x.insert(x.begin()+i, filler);
			this->cols_++;
		}
	}
	
	orig_x.resize(x.size());
	std::copy(x.begin(), x.end(), &(orig_x[0]));
	
	return &(this->referenced_);
}





/* ---------- AddDiagonal ---------- */
template<class T> bool SparseMatrix<T>::AddDiagonal(const std::size_t& dsize, const T& dvalue)
{
	this->unfix();
	
	if ( !dvalue || !dsize )
	{
		return true;
	}
	
    std::size_t ddsize = dsize;
	if ( this->referenced_.size() )
	{
		ddsize = 0;
		for (std::size_t i=0; i<this->referenced_.size(); i++)
		{
			if ( this->referenced_[i] )
			{
				ddsize++;
			}
		}
	}
	
    std::size_t os = this->velements_.size();
	this->velements_.resize(this->velements_.size()+ddsize);
	if ( this->referenced_.size() )
	{
		for (std::size_t i=0; i<dsize; i++)
		{
			if ( this->referenced_[i] )
			{
				this->velements_[os++] = i;
			}
		}
	}
	else
	{
		for (std::size_t i=0; i<ddsize; i++)
		{
			this->velements_[os+i] = i;
		}
	}
	
	this->vdata_.resize(this->vdata_.size()+ddsize, dvalue);
	
	this->vne_.resize(this->vne_.size()+ddsize, 1);
	
	this->cols_ = std::max(this->cols_, dsize);
	this->rows_ += ddsize;
	
	return true;
}





/* ---------- Add1Reverse ---------- */
template<class T> bool SparseMatrix<T>::Add1Reverse(const std::valarray<std::size_t>& els,
                                                    const std::size_t& reverse_element,
                                                    const T& value)
{
	this->unfix();
	
	if ( !value || !els.size() )
	{
		return true;
	}
	
    std::size_t os = this->velements_.size();
	
	this->velements_.insert(this->velements_.end(), reverse_element);
	this->velements_.insert(this->velements_.end(), &(els[0]), &(els[els.size()]));
	std::sort(&(this->velements_[os]), &(this->velements_[this->velements_.size()]));
	
    std::size_t pos = 0;
	for (std::size_t i=os; i<this->velements_.size(); i++)
	{
		if ( this->velements_[i] == reverse_element )
		{
			pos = i;
			break;
		}
	}
	
	this->vdata_.resize(this->velements_.size(), value/T(els.size()));
	this->vdata_[pos] = -value;
	
	this->vne_.insert(this->vne_.end(), els.size()+1);
	
	this->cols_ = std::max(this->cols_, els.max()+1);
	this->rows_++;
	
	return true;
}





/* ---------- AddNReverse ---------- */
template<class T> bool SparseMatrix<T>::AddNReverse(const std::vector< std::valarray<std::size_t> >& els,
                                                    const std::vector<std::size_t>& r_e,
                                                    const T& value)
{
	this->unfix();
	
	if ( !value || !r_e.size() || els.size() != r_e.size() )
	{
		return true;
	}
    std::size_t total_new_els = 0;
	for (std::size_t i=0; i<els.size(); i++)
	{
		total_new_els += els[i].size() + 1;
	}
	
    std::size_t os = this->velements_.size();
	this->velements_.resize(this->velements_.size()+total_new_els);
	std::valarray<std::size_t> pos(std::size_t(0), r_e.size());
	for (std::size_t i=0; i<els.size(); i++)
	{
		this->velements_[os] = r_e[i];
		std::copy(&(els[i][0]), &(els[i][els[i].size()]), &(this->velements_[os+1]));
		std::sort(&(this->velements_[os]), &(this->velements_[os+els[i].size()+1]));
		for (std::size_t j=os; j<os+els[i].size()+1; j++)
		{
			if ( this->velements_[j] == r_e[i] )
			{
				pos[i] = j;
				break;
			}
		}
		os += els[i].size() + 1;
		this->cols_ = std::max(this->cols_, els[i].max()+1);
	}
	this->rows_ += r_e.size();
	
	os = this->vdata_.size();
	this->vdata_.resize(this->velements_.size(), T(1));
	for (std::size_t i=0; i<els.size(); i++)
	{
		std::fill(this->vdata_.begin()+os-1, this->vdata_.begin()+os+els[i].size()+1, -value);
		this->vdata_[pos[i]] = value * T(els[i].size());
		os += els[i].size() + 1;
	}
	
	os = this->vne_.size();
	this->vne_.resize(this->vne_.size()+els.size(), 0);
	for (std::size_t i=0; i<els.size(); i++)
	{
		this->vne_[os+i] = els[i].size() + 1;
	}
	
	return true;
}





/* ---------- operator() ---------- */
template<class T> T SparseMatrix<T>::operator()(const std::size_t& irow, const std::size_t& icol)
{
	this->fix();
    std::size_t pos = 0;
	for (std::size_t i=0; i<irow; i++)
	{
		pos += this->ne_[i];
	}
	std::valarray<std::size_t> els = this->elements_[std::slice(pos,this->ne_[irow],1)];
	for (std::size_t i=0; i<els.size(); i++)
	{
		if ( els[i] == icol )
		{
			return this->data_[pos+i];
		}
	}
	return T(0);
}

/* ---------- operator* ---------- */
/// compute vec = A*v
template<class T> std::valarray<T> SparseMatrix<T>::operator*(const std::valarray<T>& v)
{
	this->fix();
	
	std::valarray<T> b(T(0), this->rows_);
	
    std::size_t position = 0;
	for (std::size_t i=0; i<this->rows_; i++)
	{
		if ( !this->ne_[i] )
		{
			continue;
		}
        const std::valarray<T> tmp_ne_ = this->data_[std::slice(position, this->ne_[i], 1)];
		std::valarray<T> tmp = v[this->elements_[std::slice(position, this->ne_[i], 1)]];
		tmp *= tmp_ne_;
		b[i] = tmp.sum();
		position += this->ne_[i];
	}
	
	return b;
}





/* ---------- multiply ---------- */
/// compute u += A*v
template<class T> void SparseMatrix<T>::multiply(std::valarray<T>& u, const std::valarray<T>& v)
{
	this->fix();
	
    std::size_t position = 0;
	for (std::size_t i=0; i<u.size(); i++)
	{
		if ( !this->ne_[i] )
		{
			continue;
		}
		std::valarray<T> tmp = std::valarray<T>(this->data_[std::slice(position, this->ne_[i], 1)]) *
		std::valarray<T>(v[this->elements_[std::slice(position, this->ne_[i], 1)]]);
		u[i] += tmp.sum();
		position += this->ne_[i];
	}
}





/* ---------- multiplyT ---------- */
/// compute v += A(transpose) * u
template<class T> void SparseMatrix<T>::multiplyT(std::valarray<T>& v, const std::valarray<T>& u)
{
	this->fix();
	
    std::size_t position = 0;
	for (std::size_t i=0; i<this->rows_; i++)
	{
		if ( !this->ne_[i] )
		{
			continue;
		}
		for (std::size_t j=0; j<this->ne_[i]; j++)
		{
			v[this->elements_[position]] += this->data_[position] * u[i];
			position++;
		}
	}
}





/* ---------- norm ---------- */
/// compute the norm of the values of the matrix
template<class T> T SparseMatrix<T>::norm() const
{
	T ave = T(0);
	if ( this->fixed_ )
	{
		ave = std::sqrt( ( this->data_ * this->data_ ).sum() );
	}
	else
	{
		for (std::size_t i=0; i<this->vdata_.size(); i++)
		{
			ave += this->vdata_[i] * this->vdata_[i];
		}
		ave = std::sqrt(ave);
	}
	return ave;
}





/* ---------- non_zeros ---------- */
/// number of non-zeros is the sparse matrix
template<class T> std::size_t SparseMatrix<T>::non_zeros() const
{
	return ( this->fixed_ ) ? this->data_.size() : this->vdata_.size();
}





/* ---------- iterate ---------- */
/// iteratively solve A*x = u where itmax is the maximum number of iterations to perform and conv determines the order of magnitude convergence before the system is considered stable
template<class T> std::valarray<T> SparseMatrix<T>::iterate(std::valarray<T>& u, const std::size_t& itmax, const std::size_t& conv)
{
    /**
     subroutine to solve the linear problem Ax=u using the lsqr algorithm.
     
     reference: C.C.Paige and M.A.Saunders, ACM Trans.Math.Softw. 8, 43-71, 1982 and ACM Trans.Math.Softw. 8, 195-209, 1982.
     See also A.v.d.Sluis and H.v.d. Vorst in: G. Nolet (ed.), Seismic Tomography, Reidel, 1987.
     
     Input: u contains the data (is overwritten), itmax is the number of iterations.
     Output: x is the solution;
     Scratch: arrays v(n) and w(n)
     */
    
	this->fix();
	
	std::valarray<T> x(T(0), this->cols_);
	
	if ( u.size() != this->rows_ )
	{
		x.resize(0);
		debug("SparseMatrix::iterate","bad dimensions, cannot invert");
		throw;
        return x;
	}
	
	debug(" iter\t rms    \t rel   ");
	
	std::valarray<T> v(T(0), this->cols_);
	
	T beta = normalise(u);
	T b1 = beta;
	this->multiplyT(v,u);
	T alpha = normalise(v);
	T rhobar = alpha;
	T phibar = beta;
	
	std::valarray<T> w = v;
	
	debugnr(" ");
	std::string dout = "\t" + std::to_string(0) + "  " + std::to_string(beta) + "     " + std::to_string(T(1));
	debugnr(dout);
	
	T initial_rms = T(0);
	
	for (std::size_t iter=0; iter<itmax; iter++)
	{ // repeat for fixed number of iterations itmax
		u *= -alpha; // bidiagonalisation
		this->multiply(u,v);
		beta = normalise(u);
		v *= -beta;
		this->multiplyT(v,u);
		alpha = normalise(v);
		// modified QR factorisation
		T rho = std::sqrt(rhobar * rhobar + beta * beta);
		T c = rhobar / rho;
		T s = beta / rho;
		rhobar = -c * alpha;
		T phi = c * phibar;
		phibar *= s;
		
		// update solution x and storage vector w
		x += w * (phi / rho);
		w = w * (-(s * alpha) / rho) + v;
		
		debugnr(dout.size(), "\b");
		dout = "\t" + std::to_string(iter+1) + "  " + std::to_string(phibar) + "     " + std::to_string(phibar/b1);
		debugnr(dout);
		
		if ( !iter )
		{
			initial_rms = phibar;
		}
		else
		{
			if ( phibar < initial_rms * std::pow(T(10),-int(conv)) )
			{
				debug("");
				message("iterations finished due to order " + std::to_string(conv) + " convergence");
				break;
			}
		}
	}
	debug("");
	return x;
}





#endif /* _SPARSE_MATRIX_ */
