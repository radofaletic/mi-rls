/*
  text-front
  
  a text based interface
  
  Rado Faletic
  Department of Physics
  Faculty of Science
  Australian National University  ACT  0200
  Australia
  
  Rado.Faletic@anu.edu.au
  19th May 2004
*/


#ifndef _TEXT_FRONT_
#define _TEXT_FRONT_

/* ---------- standard header files ---------- */
#include <algorithm>
#include <iostream>
#include <string>
#include <sstream>
#include <valarray>
#include <vector>
/* ---------- control flags ---------- */
#ifdef USE_MESSAGES
#include <fstream>
extern std::ofstream messages;
#endif /* USE_MESSAGES */
/* ------- user header files -------- */
#include "conversions.h"
/* ---------------------------------- */


/* ---------------------------------- */


/* ------------------------------------------- */
/* ---------- function declarations ---------- */
/* ------------------------------------------- */


void writestuff(const std::string&);
void writestufff(const std::string&);

void message(const std::string&);

void debug(const std::string&, const std::string&);
void debug(const std::string&);
void debugn(const std::string&);
void debugn(const size_t&, const std::string&);
void debugnr(const std::string&);
void debugnr(const size_t&, const std::string&);

void counter(const std::string&, const size_t&);
void counter(const std::string&, const size_t&, const size_t&);

bool yesno(const std::string&, const bool& = false);

template<class T> T question(const std::string&, const T&);
template<class T> T question(const std::string&, const std::vector<T>&, const T&);
template<class T> T question(const std::string&, const std::valarray<T>&, const T&);


/* ------------------------------------------ */
/* ---------- function definitions ---------- */
/* ------------------------------------------ */


/* ---------- writestuff ---------- */
void
writestuff(const std::string& stuff)
{
  std::cout << stuff;
#ifdef USE_MESSAGES
  messages << stuff;
#endif /* USE_MESSAGES */
}
/* -------------------------------- */

/* ---------- writestufff ---------- */
void
writestufff(const std::string& stuff)
{
  writestuff(stuff);

  std::cout << std::endl;
#ifdef USE_MESSAGES
  messages << std::endl;
#endif /* USE_MESSAGES */
}
/* --------------------------------- */

/* ---------- message ---------- */
void
message(const std::string& message)
{
  writestufff(message);
}
/* ----------------------------- */

/* ---------- debug ---------- */
void inline
debug(const std::string& function, const std::string& message)
{
#if DEBUG > 0

  std::clog << "DEBUG::"+function+":\n\t"+message << std::endl;
#ifdef USE_MESSAGES
  messages << "DEBUG::"+function+":\n\t"+message << std::endl;
#endif /* USE_MESSAGES */

#endif /* DEBUG */
}
/* --------------------------- */

/* ---------- debug ---------- */
void inline
debug(const std::string& message)
{
#if DEBUG > 0

  if ( message != "" && message != "\n" )
    {
      std::clog << "\t";
#ifdef USE_MESSAGES
	  messages << "\t";
#endif /* USE_MESSAGES */
    }
  std::clog << message << std::endl;
#ifdef USE_MESSAGES
  messages << message << std::endl;
#endif /* USE_MESSAGES */

#endif /* DEBUG */
}
/* --------------------------- */

/* ---------- debugn ---------- */
void inline
debugn(const std::string& message)
{
#if DEBUG > 0

  std::string s = message;
  std::string::size_type pb = 0;
  while ( pb != s.size() )
    {
      if ( s[pb] == '\n' )
	{
	  s.replace(pb,1,"\n\t");
	}
      pb++;
    }

  std::clog << s;
#ifdef USE_MESSAGES
  messages << s;
#endif /* USE_MESSAGES */

#endif /* DEBUG */
}
/* ---------------------------- */

/* ---------- debugn ---------- */
void inline
debugn(const size_t& num, const std::string& message)
{
#if DEBUG > 0

  std::string my_message = "";
  for (size_t i=0; i<num; i++)
    {
      my_message += message;
    }
  debugn(my_message);

#endif /* DEBUG */
}
/* ---------------------------- */

/* ---------- debugnr ---------- */
void inline
debugnr(const std::string& message)
{
#if DEBUG > 0

  std::string s = message;
  std::string::size_type pb = 0;
  while ( pb != s.size() )
    {
      if ( s[pb] == '\n' )
	{
	  s.replace(pb,1,"\n\t");
	}
      pb++;
    }

  std::clog << s;

#endif /* DEBUG */
}
/* ----------------------------- */

/* ---------- debugnr ---------- */
void inline
debugnr(const size_t& num, const std::string& message)
{
#if DEBUG > 0

  std::string my_message = "";
  for (size_t i=0; i<num; i++)
    {
      my_message += message;
    }
  debugnr(my_message);

#endif /* DEBUG */
}
/* ----------------------------- */

/* ---------- counter ---------- */
void inline
counter(const std::string& s, const size_t& val)
{
  std::cerr << s << "[";
  std::cerr.width(7);
  std::cerr << val;
  std::cerr << "]" << "\b\b\b\b\b\b\b\b\b";
  for (size_t i=0; i<s.size(); i++)
    {
      std::cerr << "\b";
    }
}
/* ----------------------------- */

/* ---------- counter ---------- */
void inline
counter(const std::string& s, const size_t& sz, const size_t& val)
{
  static size_t counter_percentage = 0;
  if ( val == 0 || val == 1 )
    {
      counter_percentage = 99999;
    }
  size_t np = ( 100 * val ) / sz;
  if ( np != counter_percentage )
    {
      counter_percentage = np;
      std::cerr << s << " ";
      std::cerr.width(3);
      std::cerr << counter_percentage;
      std::cerr << "%" << "\b\b\b\b\b";
      for (size_t i=0; i<s.size(); i++)
	{
	  std::cerr << "\b";
	}
    }
}
/* ----------------------------- */

/* ---------- yesno ---------- */
bool
yesno(const std::string& question, const bool& init_value)
{
  writestuff(question+": [");
  if ( init_value )
    {
      writestuff("Yes] ");
    }
  else
    {
      writestuff("No] ");
    }

  std::string response = "";
  std::getline(std::cin, response);
  if ( response == std::string("Y") ||
       response == std::string("y") ||
       response == std::string("YES") ||
       response == std::string("Yes") ||
       response == std::string("yes") ||
       response == std::string("T") ||
       response == std::string("t") ||
       response == std::string("TRUE") ||
       response == std::string("True") ||
       response == std::string("true") )
    {
#ifdef USE_MESSAGES
      writestufff("true");
#endif /* USE_MESSAGES */
      return true;
    }
  else if ( response == std::string("N") ||
	    response == std::string("n") ||
	    response == std::string("NO") ||
	    response == std::string("No") || 
	    response == std::string("no") ||
	    response == std::string("F") ||
	    response == std::string("f") ||
	    response == std::string("FALSE") ||
	    response == std::string("False") || 
	    response == std::string("false") )
    {
#ifdef USE_MESSAGES
      writestufff("false");
#endif /* USE_MESSAGES */
      return false;
    }
#ifdef USE_MESSAGES
  std::string tf = ( init_value ) ? "true" : "false";
  writestufff(tf);
#endif /* USE_MESSAGES */
  return init_value;
}
/* --------------------------- */

/* ---------- question ---------- */
template<class T> T
question(const std::string& question, const T& init_value)
{
  writestuff(question+": ["+ntos(init_value)+"] ");

  std::string response = "";
  std::getline(std::cin,response);
  if ( response != std::string("") &&
       response != std::string(" ") &&
       response != std::string("\n") &&
       response != std::string("\0") &&
       response != std::string("\t") )
    {
      T ret;
      std::istringstream word(response);
      word >> ret;
#ifdef USE_MESSAGES
      writestufff(ntos(ret));
#endif /* USE_MESSAGES */
      return ret;
    }
#ifdef USE_MESSAGES
  writestufff(ntos(init_value));
#endif /* USE_MESSAGES */
  return init_value;
}
/* ------------------------------ */

/* ---------- question ---------- */
template<class T> T
question(const std::string& question, const std::vector<T>& list, const T& init_value)
{
  for (size_t i=0; i<list.size(); i++)
    {
      if ( i < 9 )
	{
	  writestuff("  ");
	}
      else if ( i < 99 )
	{
	  writestuff(" ");
	}
      writestufff(ntos(i+1)+". "+ntos(list[i]));
    }
  writestuff(question+": ["+ntos(init_value)+"] ");

  std::string response = "";
  std::getline(std::cin,response);
  for (size_t i=0; i<list.size(); i++)
    {
      if ( response == ntos(i+1) || response == std::string(list[i]) )
	{
#ifdef USE_MESSAGES
	  writestufff(ntos(list[i]));
#endif /* USE_MESSAGES */
	  return list[i];
	}
    }
#ifdef USE_MESSAGES
  writestufff(ntos(init_value));
#endif /* USE_MESSAGES */
  return init_value;
}
/* ------------------------------ */

/* ---------- question ---------- */
template<class T> T
question(const std::string& tquestion, const std::valarray<T>& list, const T& init_value)
{
  std::vector<T> vlist(list.size());
  std::copy(&list[0], &list[list.size()], &vlist[0]);
  return question(tquestion, vlist, init_value);
}
/* ------------------------------ */


#endif /* _TEXT_FRONT_ */
