/*
  fortran_io
  
  Fortran i/o routines written in C++
  
  Rado Faletic
  Department of Physics
  Faculty of Science
  Australian National University  ACT  0200
  Australia
  
  Rado.Faletic@anu.edu.au
  20th June 2004
*/


#ifndef _FORTRAN_IO_
#define _FORTRAN_IO_


/* ---------- standard header files ---------- */
#include <algorithm>
#include <fstream>
#include <string>
#include <typeinfo>
#include <valarray>
#include <vector>
/* ---------- user header files ---------- */
#include "file.h"
#include "front-end.h"
/* ---------------------------------- */


/* ------------------------------------------- */
/* ---------- function declarations ---------- */
/* ------------------------------------------- */


namespace Fortran
{
  //
  // functions typically follow the syntax
  //
  //   for integers:
  //     ifunc_name(iofile, value, byte_swapping, format);
  //
  //   for real:
  //     ffunc_name(iofile, value, byte_swapping, format, precision);
  //

  enum io_type
    {
      io_null,
      io_int,
      io_short,
      io_long,
      io_size_t,
      io_float,
      io_double,
      io_long_double
    };

  io_type integer_type();

  io_type real_type(const dataprecision& = Single);

  template<class T> bool is_integer(const T&);

  template<class T> bool is_real(const T&);

  template<class T> void read_byte_info(std::ifstream&, T&, const bool& = false, const dataformat& = Unformatted);

  template<class T> void write_byte_info(std::ofstream&, const T& = T(0), const bool& = false, const dataformat& = Unformatted);

  template<class T> void iread_item(std::ifstream&, T&, const bool& = false, const dataformat& = Formatted);

  void bread_item(std::ifstream&, bool&, const bool& = false, const dataformat& = Formatted);

  template<class T> void iwrite_item(std::ofstream&, const T& = T(0), const bool& = false, const dataformat& = Formatted);

  void bwrite_item(std::ofstream&, const bool& = false, const bool& = false, const dataformat& = Formatted);

  void iread_blank(std::ifstream&, const dataformat& = Formatted);

  template<class T> void iread_line(std::ifstream&, T&, bool&, const dataformat& = Formatted );

  template<class T> void iwrite_line(std::ofstream&, const T& = T(0), const bool& = false, const dataformat& = Formatted);

  template<class T> void iread_vector(std::ifstream&, T*,  T*, const bool& = false, const dataformat& = Formatted);

  template<class T> void iwrite_vector(std::ofstream&, T*, T*, const bool& = false, const dataformat& = Formatted);

  template<class T> void iread_vector(std::ifstream&, const size_t&, std::valarray<T>&, bool&, const dataformat& = Formatted);
  template<class T> void iread_vector(std::ifstream&, const size_t&, std::vector<T>&, bool&, const dataformat& = Formatted);

  template<class T> void iwrite_vector(std::ofstream&, const std::valarray<T>&, const bool& = false, const dataformat& = Formatted);

  template<class T> void peek(std::ifstream&, T&, const bool& = false, const dataformat& = Formatted);

  template<class T> bool fcheckset(std::ifstream&, const T&, bool&, const dataformat&, dataprecision&);

  template<class T> void fread_item(std::ifstream&, T&, const bool& = false, const dataformat& = Formatted, const dataprecision& = Single);

  template<class T> void fwrite_item(std::ofstream&, const T&, const bool& = false, const dataformat& = Formatted, const dataprecision& = Single);

  void fread_blank(std::ifstream&, const dataformat& = Formatted, const dataprecision& = Single);

  template<class T> void fread_line(std::ifstream&, T&, bool&, const dataformat&, dataprecision&);

  template<class T> void fwrite_line(std::ofstream&, const T&, const bool& = false, const dataformat& = Formatted, const dataprecision& = Single);

  template<class T> void fread_vector(std::ifstream&, T*, T*, const bool& = false, const dataformat& = Formatted, const dataprecision& = Single);

  template<class T> void fwrite_vector(std::ofstream&, T*, T*, const bool& = false, const dataformat& = Formatted, const dataprecision& = Single);

  template<class T> void fread_vector(std::ifstream&, const size_t&, std::valarray<T>&, bool&, const dataformat&, dataprecision&);
  template<class T> void fread_vector(std::ifstream&, const size_t&, std::vector<T>&, bool&, const dataformat&, dataprecision&);

  template<class T> void fwrite_vector(std::ofstream&, const std::valarray<T>&, const bool& = false, const dataformat& = Formatted, const dataprecision& = Single);

  template<class T> void fretrieve_vector(const std::string&, std::valarray<T>&, bool&, const dataformat& = Formatted);

  template<class T> void fdump_vector(const std::string&, const std::valarray<T>&, const bool& = false, const dataformat& = Formatted);
}


/* ------------------------------------------ */
/* ---------- function definitions ---------- */
/* ------------------------------------------ */


/* ---------- integer_type ---------- */
/* what type of integer is 4 bytes
   in size?                           */
Fortran::io_type
Fortran::integer_type()
{
  if ( sizeof(int) == 4 )
    {
      return io_int;
    }
  else if ( sizeof(short) == 4 )
    {
      return io_short;
    }
  else if ( sizeof(long) == 4 )
    {
      return io_long;
    }
  else if ( sizeof(size_t) == 4 )
    {
      return io_size_t;
    }
  return io_null;
}
/* ---------------------------------- */

/* ---------- real_type ---------- */
/* what type of float is 4 bytes
   (or 8 bytes for double
   precision) in size?             */
Fortran::io_type
Fortran::real_type(const dataprecision& precision)
{
  unsigned int fsize = ( precision == Double ) ? 8 : 4;
  if ( sizeof(float) == fsize )
    {
      return io_float;
    }
  else if ( sizeof(double) == fsize )
    {
      return io_double;
    }
  else if ( sizeof(long double) == fsize )
    {
      return io_long_double;
    }
  return io_null;
}
/* ------------------------------- */

/* ---------- is_integer ---------- */
/* is `value' an integer?           */
template<class T> bool
Fortran::is_integer(const T& value)
{
  if ( typeid(short) == typeid(T) ||
       typeid(int) == typeid(T) ||
       typeid(long) == typeid(T) ||
       typeid(unsigned short) == typeid(T) ||
       typeid(unsigned int) == typeid(T) ||
       typeid(unsigned long) == typeid(T) ||
       typeid(size_t) == typeid(T) ||
       typeid(std::streamsize) == typeid(T) )
    {
      return true;
    }
  return false;
}
/* -------------------------------- */

/* ---------- is_real ---------- */
/* is `value' a float?           */
template<class T> bool
Fortran::is_real(const T& value)
{
  if ( typeid(float) == typeid(T) ||
       typeid(double) == typeid(T) ||
       typeid(long double) == typeid(T) )
    {
      return true;
    }
  return false;
}
/* ------------------------------ */

/* ---------- read_byte_info ---------- */
/* read 4 bytes of information in
   binary, generally used for byte
   information about a block of data    */
template<class T> void
Fortran::read_byte_info(std::ifstream& file,
			T& value,
			const bool& byte_swapping,
			const dataformat& format)
{
  value = T(0);
  if ( format != Unformatted )
    {
      return;
    }
  else if ( !file )
    {
      debug("Fortran::read_byte_info", "no file to read from");
      throw; return;
    }

  short stemp;
  int itemp;
  long ltemp;

  // read the byte_information
  switch(Fortran::integer_type())
    {
    case io_short:
      file.read((char*)& stemp, 4);
      if ( byte_swapping )
	{
	  byte_swap(stemp);
	}
      value = stemp;
      break;
    case io_int:
      file.read((char*)& itemp, 4);
      if ( byte_swapping )
	{
	  byte_swap(itemp);
	}
      value = itemp;
      break;
    case io_long:
      file.read((char*)& ltemp, 4);
      if ( byte_swapping )
	{
	  byte_swap(ltemp);
	}
      value = ltemp;
      break;
    }
}
/* ------------------------------------ */

/* ---------- write_byte_info ---------- */
/* write 4 bytes of information in
   binary, generally used for byte
   information about a block of data     */
template<class T> void
Fortran::write_byte_info(std::ofstream& file,
			 const T& value,
			 const bool& byte_swapping,
			 const dataformat& format)
{
  if ( format != Unformatted )
    {
      return;
    }
  else if ( !file )
    {
      debug("Fortran::write_byte_info","no file to write to");
      throw; return;
    }

  short stemp;
  int itemp;
  long ltemp;

  switch(Fortran::integer_type())
    {
    case io_short:
      stemp = short(value);
      if ( byte_swapping )
	{
	  byte_swap(stemp);
	}
      file.write((char*)& stemp, 4);
      break;
    case io_int:
      itemp = int(value);
      if ( byte_swapping )
	{
	  byte_swap(itemp);
	}
      file.write((char*)& itemp, 4);
      break;
    case io_long:
      ltemp = long(value);
      if ( byte_swapping )
	{
	  byte_swap(ltemp);
	}
      file.write((char*)& ltemp, 4);
      break;
    }

}
/* ------------------------------------- */

/* ---------- iread_item ---------- */
/* read a single integer value      */
template<class T> void
Fortran::iread_item(std::ifstream& file,
		    T& value,
		    const bool& byte_swapping,
		    const dataformat& format)
{
  value = T(0);
  if ( !file )
    {
      debug("Fortran::iread_item","no file to read from");
      throw; return;
    }
  switch(format)
    {
    case Unformatted: case Binary:
      Fortran::read_byte_info(file, value, byte_swapping, Unformatted);
      break;
    case Formatted:
      file >> value;
      break;
    }
}
/* -------------------------------- */

/* ---------- bread_item ---------- */
/* read a bool value                */
void
Fortran::bread_item(std::ifstream& file,
		    bool& value,
		    const bool& byte_swapping,
		    const dataformat& format)
{
  value = false;
  if ( !file )
    {
      debug("Fortran::bread_item","no file to read from");
      throw; return;
    }
  short ivalue = 0;
  switch(format)
    {
    case Unformatted: case Binary:
      Fortran::read_byte_info(file, ivalue, byte_swapping, Unformatted);
      break;
    case Formatted:
      file >> ivalue;
      break;
    }
  if ( ivalue != 0 )
    {
      value = true;
    }
}
/* -------------------------------- */

/* ---------- iwrite_item ---------- */
/* write a single integer value      */
template<class T> void
Fortran::iwrite_item(std::ofstream& file,
		     const T& value,
		     const bool& byte_swapping,
		     const dataformat& format)
{
  if ( !file )
    {
      debug("Fortran::iwrite_item","no file to write to");
      throw; return;
    }
  switch(format)
    {
    case Unformatted: case Binary:
      Fortran::write_byte_info(file, value, byte_swapping, Unformatted);
      break;
    case Formatted:
      file << value << " ";
      break;
    }
}
/* --------------------------------- */

/* ---------- bwrite_item ---------- */
/* write a bool value                */
void
Fortran::bwrite_item(std::ofstream& file,
		     const bool& value,
		     const bool& byte_swapping,
		     const dataformat& format)
{
  if ( !file )
    {
      debug("Fortran::bwrite_item","no file to write to");
      throw; return;
    }
  short ivalue = ( value ) ? 1 : 0;
  switch(format)
    {
    case Unformatted: case Binary:
      Fortran::write_byte_info(file, ivalue, byte_swapping, Unformatted);
      break;
    case Formatted:
      file << ivalue << " ";
      break;
    }
}
/* --------------------------------- */

/* ---------- iread_blank ---------- */
/* read an integer, but do not record
   the value                         */
void Fortran::iread_blank(std::ifstream& file,
			  const dataformat& format)
{
  if ( !file )
    {
      debug("Fortran::iread_blank","no file to read from");
      throw; return;
    }
  short tmp;
  Fortran::iread_item(file, tmp, false, format);
}
/* --------------------------------- */

/* ---------- iread_line ---------- */
/* read a line than only contains a
   single integer                   */
template<class T> void
Fortran::iread_line(std::ifstream& file, T& value, bool& byte_swapping, const dataformat& format)
{
  value = T(0);

  short i1, i2;

  read_byte_info(file, i1, byte_swapping, format);
  if ( format == Unformatted && i1 != 4 )
    {
      i1 = 4;
      byte_swapping = !byte_swapping;
    }
  iread_item(file, value, byte_swapping, format);
  read_byte_info(file, i2, byte_swapping, format);
  if ( i1 != i2 )
    {
      value = T(0);
      byte_swapping = !byte_swapping;
    }
}
/* -------------------------------- */

/* ---------- iwrite_line ---------- */
/* write a line that only contains a
   single integer                    */
template<class T> void
Fortran::iwrite_line(std::ofstream& file,
		     const T& value,
		     const bool& byte_swapping,
		     const dataformat& format)
{
  if ( !file )
    {
      debug("Fortran::iwrite_line","no file to write to");
      throw; return;
    }
  Fortran::write_byte_info(file, 4, byte_swapping, format);
  Fortran::iwrite_item(file, value, byte_swapping, format);
  Fortran::write_byte_info(file, 4, byte_swapping, format);
  if ( format == Formatted )
    {
      file << std::endl;
    }
}
/* --------------------------------- */

/* ---------- iread_vector ---------- */
/* read an integer vector             */
template<class T> void
Fortran::iread_vector(std::ifstream& file,
		      T* in,
		      T* out,
		      const bool& byte_swapping,
		      const dataformat& format)
{
  if ( !file )
    {
      debug("Fortran::iread_vector","no file to read from");
      throw; return;
    }
  while ( in != out )
    {
      T tmp;
      Fortran::iread_item(file, tmp, byte_swapping, format);
      *in++ = tmp;
    }
}
/* ---------------------------------- */

/* ---------- iwrite_vector ---------- */
/* write an integer vector             */
template<class T> void
Fortran::iwrite_vector(std::ofstream& file,
		       T* in,
		       T* out,
		       const bool& byte_swapping,
		       const dataformat& format)
{
  if ( !file )
    {
      debug("Fortran::iwrite_vector","no file to write to");
      throw; return;
    }
  while ( in != out )
    {
      T tmp = *in++;
      Fortran::iwrite_item(file, tmp, byte_swapping, format);
    }
}
/* ----------------------------------- */

/* ---------- iread_vector ---------- */
/* automatically read an integer
   vector from an entire line         */
template<class T> void
Fortran::iread_vector(std::ifstream& file,
		      const size_t& vec_size,
		      std::valarray<T>& vec,
		      bool& byte_swapping,
		      const dataformat& format)
{
  vec.resize(0);
  if ( !vec_size )
    {
      return;
    }
  else if ( !file )
    {
      debug("Fortran::iread_vector","no file to read from");
      throw; return;
    }

  size_t isize;
  Fortran::read_byte_info(file, isize, byte_swapping, format);
  if ( format == Unformatted && isize != 4*vec_size )
    {
      isize = 4*vec_size;
      byte_swapping = !byte_swapping;
    }

  vec.resize(vec_size);

  for (size_t i=0; i<vec.size(); i++)
    {
      Fortran::iread_item(file, vec[i], byte_swapping, format);
    }

  size_t tsize;
  Fortran::read_byte_info(file, tsize, byte_swapping, format);
  if ( tsize != isize )
    {
      vec.resize(0);
      byte_swapping = !byte_swapping;
    }
}
/* ---------------------------------- */

/* ---------- iread_vector ---------- */
template<class T> void
Fortran::iread_vector(std::ifstream& file,
		      const size_t& vec_size,
		      std::vector<T>& vec,
		      bool& byte_swapping,
		      const dataformat& format)
{
  std::valarray<T> rvec;
  Fortran::iread_vector(file, vec_size, rvec, byte_swapping, format);
  vec.resize(rvec.size());
  std::copy(&rvec[0], &rvec[rvec.size()], vec.begin());
}
/* ---------------------------------- */

/* ---------- iwrite_vector ---------- */
/* automatically write an integer
   vector to an entire line            */
template<class T> void
Fortran::iwrite_vector(std::ofstream& file,
		       const std::valarray<T>& vec,
		       const bool& byte_swapping,
		       const dataformat& format)
{
  if ( !file || !vec.size() )
    {
      debug("Fortran::iwrite_vector","no file to write to");
      throw; return;
    }
  Fortran::write_byte_info(file, vec.size()*4, byte_swapping, format);
  for (size_t i=0; i<vec.size(); i++)
    {
      Fortran::iwrite_item(file, vec[i], byte_swapping, format);
    }
  Fortran::write_byte_info(file, vec.size()*4, byte_swapping, format);
  if ( format == Formatted )
    {
      file << std::endl;
    }
}
/* ----------------------------------- */

/* ---------- peek ---------- */
/* take a look at the next
   value, and compare it with
   `value', then put back the
   gotten characters          */
template<class T> void
Fortran::peek(std::ifstream& file, T& value, const bool& byte_swapping, const dataformat& format)
{
  value = T(0);
  if ( format != Unformatted )
    {
      return;
    }
  else if ( !file )
    {
      debug("Fortran::peek","no file to read from");
      throw; return;
    }
  Fortran::iread_item(file, value, byte_swapping, format);
  file.unget();
  file.unget();
  file.unget();
  file.unget();
}
/* -------------------------- */

/* ---------- fcheckset ---------- */
/* check that you have set the
   correct endieness and float
   precision for reading your file */
template<class T> bool
Fortran::fcheckset(std::ifstream& file, const T& value, bool& byte_swapping, const dataformat& format, dataprecision& precision)
{
  if ( format == Formatted )
    {
      return true;
    }
  else if ( !file )
    {
      debug("Fortran::fcheckset","no file to read from");
      throw; return false;
    }
  T temp;

  Fortran::peek(file, temp, byte_swapping, Unformatted);

  if ( temp != value )
    {
      if ( precision == Single && temp == 2 * value )
	{
	  precision = Double;
	  return true;
	}
      else if ( precision == Double && 2 * temp == value )
	{
	  precision = Single;
	  return true;
	}
      byte_swapping = !byte_swapping;
      Fortran::peek(file, temp, byte_swapping, Unformatted);
      if ( temp != value )
	{
	  if ( precision == Single && temp == 2 * value )
	    {
	      precision = Double;
	      return true;
	    }
	  else if ( precision == Double && 2 * temp == value )
	    {
	      precision = Single;
	      return true;
	    }
	  byte_swapping = !byte_swapping;
	  return false;
	}
    }
  return true;
}
/* ------------------------------- */

/* ---------- fread_item ---------- */
/* read a float                     */
template<class T> void
Fortran::fread_item(std::ifstream& file,
		    T& value,
		    const bool& byte_swapping,
		    const dataformat& format,
		    const dataprecision& precision)
{
  value = T(0);
  if ( !file )
    {
      debug("Fortran::fread_item","no file to read from");
      throw; return;
    }
  std::streamsize fsize = ( precision == Double ) ? 8 : 4;

  float ftemp;
  double dtemp;
  long double ltemp;

  switch(format)
    {
    case Unformatted: case Binary:
      switch(Fortran::real_type(precision))
	{
	case io_float:
	  file.read((char*)& ftemp, fsize);
	  if ( byte_swapping )
	    {
	      byte_swap(ftemp);
	    }
	  value = ftemp;
	  break;
	case io_double:
	  file.read((char*)& dtemp, fsize);
	  if ( byte_swapping )
	    {
	      byte_swap(dtemp);
	    }
	  value = dtemp;
	  break;
	case io_long_double:
	  file.read((char*)& ltemp, fsize);
	  if ( byte_swapping )
	    {
	      byte_swap(ltemp);
	    }
	  value = ltemp;
	  break;
	}
      break;
    case Formatted:
      file >> value;
      break;
    }
}
/* -------------------------------- */

/* ---------- fwrite_item ---------- */
/* write a float                     */
template<class T> void
Fortran::fwrite_item(std::ofstream& file,
		     const T& value,
		     const bool& byte_swapping,
		     const dataformat& format,
		     const dataprecision& precision)
{
  if ( !file )
    {
      debug("Fortran::fwrite_item","no file to write to");
      throw; return;
    }
  std::streamsize fsize = ( precision == Double ) ? 8 : 4;

  float ftemp;
  double dtemp;
  long double ltemp;

  switch(format)
    {
    case Unformatted: case Binary:
      switch(Fortran::real_type(precision))
	{
	case io_float:
	  ftemp = value;
	  if ( byte_swapping )
	    {
	      byte_swap(ftemp);
	    }
	  file.write((char*)& ftemp, fsize);
	  break;
	case io_double:
	  dtemp = value;
	  if ( byte_swapping )
	    {
	      byte_swap(dtemp);
	    }
	  file.write((char*)& dtemp, fsize);
	  break;
	case io_long_double:
	  ltemp = value;
	  if ( byte_swapping )
	    {
	      byte_swap(ltemp);
	    }
	  file.write((char*)& ltemp, fsize);
	  break;
	}
      break;
    case Formatted:
      if ( value == T(long(value)) )
	{
	  file.setf(std::ios_base::fmtflags(0), std::ios_base::floatfield);
	}
      else
	{
	  file.precision(8);
	  file.setf(std::ios_base::scientific, std::ios_base::floatfield);
	}
      file << value << ' ';
      break;
    }
}
/* --------------------------------- */

/* ---------- fread_blank ---------- */
/* read a float, but do not store the
   value                             */
void
Fortran::fread_blank(std::ifstream& file,
		     const dataformat& format,
		     const dataprecision& precision)
{
  if ( !file )
    {
      debug("Fortran::fread_blank","no file to read from");
      throw; return;
    }
  float tmp;
  fread_item(file, tmp, false, format, precision);
}
/* --------------------------------- */

/* ---------- fread_line ---------- */
/* read an entire line, that
   consists of a single float       */
template<class T> void
Fortran::fread_line(std::ifstream& file,
		    T& value,
		    bool& byte_swapping,
		    const dataformat& format,
		    dataprecision& precision)
{
  value = T(0);
  if ( !file )
    {
      debug("Fortran::fread_line","no file to read from");
      throw; return;
    }

  // check that the file attributes work
  std::streamsize fsize = ( precision == Double ) ? 8 : 4;
  if ( !Fortran::fcheckset(file, fsize, byte_swapping, format, precision) )
    {
      debug("Fortran::fread_line","bad checkset");
      throw; return;
    }
  fsize = ( precision == Double ) ? 8 : 4;
  std::streamsize htemp;

  // read in the byte header
  Fortran::read_byte_info(file, htemp, byte_swapping, format);

  // read in the data value
  Fortran::fread_item(file, value, byte_swapping, format, precision);

  // read in the byte tail
  Fortran::read_byte_info(file, fsize, byte_swapping, format);
  if ( htemp != fsize )
    {
      value = T(0);
    }
}
/* -------------------------------- */

/* ---------- fwrite_line ---------- */
/* write an entire line, that
   consists of a single float        */
template<class T> void
Fortran::fwrite_line(std::ofstream& file,
		     const T& value,
		     const bool& byte_swapping,
		     const dataformat& format,
		     const dataprecision& precision)
{
  if ( !file )
    {
      debug("Fortran::fwrite_line","no file to write to");
      throw; return;
    }
  std::streamsize fsize = ( precision == Double ) ? 8 : 4;

  Fortran::write_byte_info(file, fsize, byte_swapping, format);
  Fortran::fwrite_item(file, value, byte_swapping, format, precision);
  Fortran::write_byte_info(file, fsize, byte_swapping, format);
  if ( format == Formatted )
    {
      file << std::endl;
    }
}
/* --------------------------------- */

/* ---------- fread_vector ---------- */
/* read a float vector                */
template<class T> void
Fortran::fread_vector(std::ifstream& file,
		      T* in,
		      T* out,
		      const bool& byte_swapping ,
		      const dataformat& format,
		      const dataprecision& precision)
{
  if ( !file )
    {
      debug("Fortran::fread_vector","no file to read from");
      throw; return;
    }
  while ( in != out )
    {
      T tmp;
      Fortran::fread_item(file, tmp, byte_swapping, format, precision);
      *in++ = tmp;
    }
}
/* ---------------------------------- */

/* ---------- fwrite_vector ---------- */
/* write a float vector                */
template<class T> void
Fortran::fwrite_vector(std::ofstream& file,
		       T* in,
		       T* out,
		       const bool& byte_swapping,
		       const dataformat& format,
		       const dataprecision& precision)
{
  if ( !file )
    {
      debug("Fortran::iwrite_vector","no file to write to");
      throw; return;
    }
  while ( in != out )
    {
      T tmp = *in++;
      Fortran::fwrite_item(file, tmp, byte_swapping, format, precision);
    }
}
/* ----------------------------------- */

/* ---------- fread_vector ---------- */
/* automatically read a float vector
   from a single line                 */
template<class T> void
Fortran::fread_vector(std::ifstream& file,
		      const size_t& vec_size,
		      std::valarray<T>& vec,
		      bool& byte_swapping,
		      const dataformat& format,
		      dataprecision& precision)
{
  vec.resize(0);
  if ( !vec_size )
    {
      return;
    }
  else if ( !file )
    {
      debug("Fortran::fread_vector","no file to read from");
      throw; return;
    }

  size_t fsize = ( precision == Double ) ? 8 : 4;
  fsize *= vec_size;
  if ( !Fortran::fcheckset(file, fsize, byte_swapping, format, precision) )
    {
      debug("Fortran::fread_vector","bad checkset");
      throw; return;
    }

  Fortran::read_byte_info(file, fsize, byte_swapping, format);

  vec.resize(vec_size);

  for (size_t i=0; i<vec.size(); i++)
    {
      Fortran::fread_item(file, vec[i], byte_swapping, format, precision);
    }

  size_t tsize;
  Fortran::read_byte_info(file, tsize, byte_swapping, format);
  if ( tsize != fsize )
    {
      vec.resize(0);
    }
}
/* ---------------------------------- */

/* ---------- fread_vector ---------- */
template<class T> void
Fortran::fread_vector(std::ifstream& file,
		      const size_t& vec_size,
		      std::vector<T>& vec,
		      bool& byte_swapping,
		      const dataformat& format,
		      dataprecision& precision)
{
  std::valarray<T> rvec;
  Fortran::fread_vector(file, vec_size, rvec, byte_swapping, format, precision);
  vec.resize(rvec.size());
  std::copy(&rvec[0], &rvec[rvec.size()], vec.begin());
}
/* ---------------------------------- */

/* ---------- fwrite_vector ---------- */
/* automatically write a float vector
   to an entire line                   */
template<class T> void
Fortran::fwrite_vector(std::ofstream& file,
		       const std::valarray<T>& vec,
		       const bool& byte_swapping,
		       const dataformat& format,
		       const dataprecision& precision)
{
  if ( !vec.size() )
    {
      return;
    }
  else if ( !file )
    {
      debug("Fortran::fwrite_vector","no file to write to");
      throw; return;
    }
  size_t fsize = ( precision == Double ) ? 8 : 4;
  fsize *= vec.size();
  Fortran::write_byte_info(file, fsize, byte_swapping, format);
  for (size_t i=0; i<vec.size(); i++)
    {
      Fortran::fwrite_item(file, vec[i], byte_swapping, format, precision);
    }
  Fortran::write_byte_info(file, fsize, byte_swapping, format);
  if ( format == Formatted )
    {
      file << std::endl;
    }
}
/* ----------------------------------- */

/* ---------- fretrieve_vector ---------- */
/* read a file that contains only a single
   float vector                           */
template<class T> void
Fortran::fretrieve_vector(const std::string& filename,
			  std::valarray<T>& vec,
			  bool& byte_swapping,
			  const dataformat& format)
{
  std::ifstream file;
  switch(format)
    {
    case Unformatted: case Binary:
      file.open(filename.c_str(),std::ios_base::binary);
      break;
    default:
      file.open(filename.c_str());
      break;
    }

  size_t vec_size;
  Fortran::iread_line(file, vec_size, byte_swapping, format);
  vec.resize(vec_size);

  dataprecision precision = Single;
  size_t fsize = 4;
  if ( !Fortran::fcheckset(file, vec_size*fsize, byte_swapping, format, precision) )
    {
      debug("Fortran::fretrieve_vector","bad checkset");
      vec.resize(0);
      throw; return;
    }
  fsize = ( precision == Double ) ? 8 : 4;
  Fortran::fread_vector(file, vec_size, vec, byte_swapping, format, precision);

  file.close();
}
/* -------------------------------------- */

/* ---------- fdump_vector ---------- */
/* write a file that consists of a
   single float vector                */
template<class T> void
Fortran::fdump_vector(const std::string& filename,
		      const std::valarray<T>& vec,
		      const bool& byte_swapping,
		      const dataformat& format)
{

  std::ofstream file;
  switch(format)
    {
    case Unformatted: case Binary:
      file.open(filename.c_str(),std::ios_base::binary);
      break;
    default:
      file.open(filename.c_str());
      break;
    }

  dataprecision precision = Double;
  if ( sizeof(T) <= 4 )
    {
      precision = Single;
    }
  Fortran::iwrite_line(file, vec.size(), byte_swapping, format);
  Fortran::fwrite_vector(file, vec, byte_swapping, format, precision);

  file.close();
}
/* ---------------------------------- */


#endif /* _FORTRAN_IO_ */
