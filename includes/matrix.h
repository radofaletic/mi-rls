/**
 matrix
 
 matrix class
 
 Rado Faletic
 Department of Physics
 Faculty of Science
 Australian National University  ACT  0200
 Australia
 
 Rado.Faletic@anu.edu.au
 5th August 2004
 19th April 2022, updated to C++20
 */





#ifndef _MATRIX_
#define _MATRIX_





/* ---------- standard header files ---------- */
#include <algorithm>
#include <bit>
#include <fstream>
#include <ios>
#include <sstream>
#include <string>
#include <valarray>





/* ---------- user header files ---------- */
#include "conversions.h"
#include "extra_math.h"
#include "file.h"
#include "fortran_io.h"
#include "front-end.h"





/* ---------- class & function declaration ---------- */

template<class T> class Matrix
{
protected:
	std::size_t rows_;
    std::size_t cols_;
	std::valarray<T> data_;
public:
	Matrix() : rows_(0), cols_(0) { };
	Matrix(const std::size_t&, const std::size_t&, const std::valarray<T>&);
	Matrix(const std::size_t&);
	Matrix(const std::size_t&, const std::size_t&, const T& = T(0));
	void init(const std::size_t&, const T& = T(1));
	void init(const std::size_t&, const std::size_t&, const T& = T(0));
	void clear();
	~Matrix() { };
    std::size_t dim() const;
    std::size_t row_dim() const;
    std::size_t num_rows() const;
    std::size_t rows() const;
    std::size_t col_dim() const;
    std::size_t num_cols() const;
    std::size_t cols() const;
	std::string print(const bool& = false) const;
	bool write(const std::string&, const dataformat& = Unformatted) const;
	bool write(const std::string&, const std::valarray<T>&, const dataformat& = Unformatted) const;
	bool read(const std::string&, const dataformat& = Unformatted);
	bool read(const std::string&, std::valarray<T>&, const dataformat& = Unformatted);
	T& operator()(const std::size_t& position);
	T operator()(const std::size_t& position) const;
	T& operator()(const std::size_t& row, const std::size_t& col);
	T operator()(const std::size_t& row, const std::size_t& col) const;
	std::valarray<T> operator[](const std::slice&) const;
	void AddRow(const T& value = T(0));
	void AddRow(const std::valarray<T>&);
	void AddDiagonal(const T& = T(1));
	std::valarray<T> operator*(const std::valarray<T>&) const;
	Matrix<T> operator*(const Matrix<T>&) const;
	T norm() const;
};





/* ---------- function definitions ---------- */





/* ---------- Matrix ---------- */
template<class T> Matrix<T>::Matrix(const std::size_t& rows, const std::size_t& cols, const std::valarray<T>& data)
{
	if ( !rows || !cols || rows * cols < data.size() )
	{
		throw;
        return;
	}
	this->rows_ = rows;
	this->cols_ = cols;
	this->data_.resize(this->rows_ * this->cols_);
	this->data_[std::slice(0, data.size(), 1)] = data;
}





/* ---------- Matrix ---------- */
template<class T> Matrix<T>::Matrix(const std::size_t& sz)
{
	this->init(sz);
}





/* ---------- Matrix ---------- */
template<class T> Matrix<T>::Matrix(const std::size_t& rows, const std::size_t& cols, const T& init_value)
{
	this->init(rows, cols, init_value);
}





/* ---------- init ---------- */
template<class T> void Matrix<T>::init(const std::size_t& sz, const T& init_value)
{
	if ( !sz )
	{
		this->clear();
		throw;
        return;
	}
	this->rows_ = sz;
	this->cols_ = sz;
	this->data_.resize(this->rows_ * this->cols_,T(0));
	this->data_[std::slice(0, this->rows_, this->cols_+1)] = init_value;
}





/* ---------- init ---------- */
template<class T> void Matrix<T>::init(const std::size_t& rows, const std::size_t& cols, const T& init_value)
{
	if ( !rows || !cols )
	{
		this->clear();
		throw;
        return;
	}
	this->rows_ = rows;
	this->cols_ = cols;
	this->data_.resize(this->rows_ * this->cols_,T(0));
	this->data_[std::slice(0,std::min(this->rows_,this->cols_), this->cols_+1)] = init_value;
}





/* ---------- clear() ---------- */
template<class T> inline void Matrix<T>::clear()
{
	this->rows_ = 0;
	this->cols_ = 0;
	this->data_.resize(0);
}





/* ---------- dim() ---------- */
template<class T> inline std::size_t Matrix<T>::dim() const
{
	return std::max(this->rows_, this->cols_);
}





/* ---------- row_dim() ---------- */
template<class T> inline std::size_t Matrix<T>::row_dim() const
{
	return this->rows_;
}





/* ---------- num_rows() ---------- */
template<class T> inline std::size_t Matrix<T>::num_rows() const
{
	return this->rows_;
}





/* ---------- rows() ---------- */
template<class T> inline std::size_t Matrix<T>::rows() const
{
	return this->rows_;
}





/* ---------- col_dim() ---------- */
template<class T> inline std::size_t Matrix<T>::col_dim() const
{
	return this->cols_;
}





/* ---------- num_cols() ---------- */
template<class T> inline std::size_t Matrix<T>::num_cols() const
{
	return this->cols_;
}





/* ---------- cols() ---------- */
template<class T> inline std::size_t Matrix<T>::cols() const
{
	return this->cols_;
}





/* ---------- print ---------- */
template<class T> std::string Matrix<T>::print(const bool& CSV) const
{
	std::ostringstream tmp;
	for (std::size_t i=0; i<this->rows_; i++)
	{
		for (std::size_t j=0; j<this->cols_; j++)
		{
			tmp << this->data_[(i*(this->cols_))+j];
			if ( j < (this->cols_) - 1 )
			{
				tmp << ( CSV ? "," : "\t" );
			}
			else if ( i < (this->rows_) - 1 )
			{
				tmp << "\n";
			}
		}
	}
	return tmp.str();
}





/* ---------- write ---------- */
template<class T> inline bool Matrix<T>::write(const std::string& filename, const dataformat& format) const
{
	std::valarray<T> b(0);
	return this->write(filename, b, format);
}





/* ---------- write ---------- */
template<class T> bool Matrix<T>::write(const std::string& filename, const std::valarray<T>& b, const dataformat& format) const
{
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
	
	int hb = ( b.size() ) ? 12 : 8;
	
	// write the number of rows and columns
	Fortran::write_byte_info(matrixfile, hb, byte_swapping, format);
	Fortran::iwrite_item(matrixfile, this->rows_, byte_swapping, format);
	Fortran::iwrite_item(matrixfile, this->cols_, byte_swapping, format);
	if ( b.size() )
	{
		Fortran::iwrite_item(matrixfile, b.size(), byte_swapping, format);
	}
	Fortran::write_byte_info(matrixfile, hb, byte_swapping, format);
	
	// write the data
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
template<class T> inline bool Matrix<T>::read(const std::string& filename, const dataformat& format)
{
	std::valarray<T> b(0);
	return this->read(filename, b, format);
}





/* ---------- read ---------- */
template<class T> bool Matrix<T>::read(const std::string& filename, std::valarray<T>& b, const dataformat& format)
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
	if ( !matrixfile )
	{
		debug("Matrix<T>::read", "can't open matrix file '" + filename + "'");
		throw;
        return false;
	}
	
	// read the number of rows, columns and total number of non-zero elements
	bool byte_swapping = (std::endian::native == std::endian::little) ? true : false;
	dataprecision precision = Single;
	int hb = 8;
	if ( !Fortran::fcheckset(matrixfile, hb, byte_swapping, format, precision) )
	{
		hb = 12;
	}
	if ( !Fortran::fcheckset(matrixfile, hb, byte_swapping, format, precision) )
	{
		debug("Matrix<T>::read", "bad checkset");
		matrixfile.close();
		throw;
        return false;
	}
    std::size_t isize;
	Fortran::read_byte_info(matrixfile, isize, byte_swapping, format);
	Fortran::iread_item(matrixfile, this->rows_, byte_swapping, format);
	Fortran::iread_item(matrixfile, this->cols_, byte_swapping, format);
	int bsize = 0;
	if ( hb == 12 )
	{
		Fortran::iread_item(matrixfile, bsize, byte_swapping, format);
	}
    std::size_t itemp;
	Fortran::read_byte_info(matrixfile, itemp, byte_swapping, format);
	if ( itemp != isize )
	{
		debug("Matrix<T>::read", "error reading matrix file '"+filename+"'");
		matrixfile.close();
		this->clear();
		throw;
        return false;
	}
	
	// read the data
    std::size_t fsize = ( precision == Double ) ? 8 : 4;
	fsize *= rows_ * cols_;
	if ( !Fortran::fcheckset(matrixfile, fsize, byte_swapping, format, precision) )
	{
		debug("Matrix<T>::read", "bad checkset");
		matrixfile.close();
		this->clear();
		throw;
        return false;
	}
	if ( Fortran::is_integer(T(0)) )
	{
		Fortran::iread_vector(matrixfile, this->rows_ * this->cols_, this->data_, byte_swapping, format);
		if ( bsize )
		{
			Fortran::iread_vector(matrixfile, bsize, b, byte_swapping, format);
		}
	}
	else if ( Fortran::is_real(T(0)) )
	{
		Fortran::fread_vector(matrixfile, this->rows_ * this->cols_, this->data_, byte_swapping, format, precision);
		if ( bsize )
		{
			Fortran::fread_vector(matrixfile, bsize, b, byte_swapping, format, precision);
		}
	}
	else
	{
		debug("Matrix<T>::read", "can't determine numerical type");
		matrixfile.close();
		this->clear();
		throw;
        return false;
	}
	
	matrixfile.close();
	
	return true;
}





/* ---------- operator() ---------- */
template<class T> inline T& Matrix<T>::operator()(const std::size_t& position)
{
	return this->data_[position];
}





/* ---------- operator() ---------- */
template<class T> inline T Matrix<T>::operator()(const std::size_t& position) const
{
	return this->data_[position];
}





/* ---------- operator() ---------- */
template<class T> inline T& Matrix<T>::operator()(const std::size_t& row, const std::size_t& col)
{
	return this->data_[(row * (this->cols_)) + col];
}





/* ---------- operator() ---------- */
template<class T> inline T Matrix<T>::operator()(const std::size_t& row, const std::size_t& col) const
{
	return this->data_[(row * (this->cols_)) + col];
}





/* ---------- operator[] ---------- */
template<class T> inline std::valarray<T> Matrix<T>::operator[](const std::slice& tsl) const
{
	return this->data_[tsl];
}





/* ---------- AddRow ---------- */
template<class T> void Matrix<T>::AddRow(const T& value)
{
	std::valarray<T> tmp = this->data_;
	this->data_.resize(tmp.size() + this->cols_, value);
	this->data_[std::slice(0, tmp.size(), 1)] = tmp;
	this->rows_++;
}





/* ---------- AddRow ---------- */
template<class T> void Matrix<T>::AddRow(const std::valarray<T>& value)
{
	if ( this->cols_ < value.size() )
	{
		throw;
        return;
	}
	std::valarray<T> tmp = this->data_;
	this->data_.resize(tmp.size()+cols_, T(0));
	this->data_[std::slice(0, tmp.size(), 1)] = tmp;
	this->data_[std::slice(tmp.size(), value.size(), 1)] = value;
	this->rows_++;
}





/* ---------- AddDiagonal ---------- */
template<class T> void Matrix<T>::AddDiagonal(const T& value)
{
	std::valarray<T> tmp = this->data_;
	this->data_.resize(tmp.size() + this->cols_ * this->cols_, T(0));
	this->data_[std::slice(0, tmp.size(), 1)] = tmp;
	this->data_[std::slice(tmp.size(), this->cols_, this->cols_ + 1)] = value;
	this->rows_ += this->cols_;
}





/* ---------- operator* ---------- */
/// compute b = A*v
template<class T> std::valarray<T> Matrix<T>::operator*(const std::valarray<T>& v) const
{
	std::valarray<T> b(T(0), this->rows_);
	if ( this->cols_ < v.size() )
	{
		b.resize(0);
		throw;
        return b;
	}
	for (std::size_t i=0; i<this->rows_; i++)
	{
		std::valarray<T> tmp = this->data_[std::slice(i * this->cols_, v.size(), 1)] * v;
		b[i] = tmp.sum();
	}
	
	return b;
}





/* ---------- operator* ---------- */
/// compute    C = A*B
template<class T> Matrix<T> Matrix<T>::operator*(const Matrix<T>& B) const
{
	Matrix<T> C(this->rows_, B.col_dim());
	if ( this->cols_ < B.row_dim() )
	{
		C.clear();
		throw;
        return C;
	}
	for (std::size_t i=0; i<this->rows_; i++)
	{
		for (std::size_t j=0; j<B.col_dim(); j++)
		{
			std::valarray<T> tmp = this->data_[std::slice(i * this->cols_, B.row_dim(), 1)] * B[std::slice(j, B.row_dim(), B.col_dim())];
			C(i,j) = tmp.sum();
		}
	}
	return C;
}





/* ---------- norm ---------- */
/// compute the norm of the values of the matrix
template<class T> T Matrix<T>::norm() const
{
	return std::sqrt( ( this->data_ * this->data_ ).sum() );
}





#endif /* _MATRIX_ */
