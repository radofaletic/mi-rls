/**
 text-front
 
 a text based interface
 
 Rado Faletic
 Department of Physics
 Faculty of Science
 Australian National University  ACT  0200
 Australia
 
 Rado.Faletic@anu.edu.au
 19th May 2004
 19th April 2022, updated to C++20
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





/* ---------- function declarations ---------- */

void writestuff(const std::string&);
void writestufff(const std::string&);

void message(const std::string&);

void debug(const std::string&, const std::string&);
void debug(const std::string&);
void debugn(const std::string&);
void debugn(const std::size_t&, const std::string&);
void debugnr(const std::string&);
void debugnr(const std::size_t&, const std::string&);

void counter(const std::string&, const std::size_t&);
void counter(const std::string&, const std::size_t&, const std::size_t&);

bool yesno(const std::string&, const bool& = false);

template<class T> T question(const std::string&, const T&);
template<class T> T question(const std::string&, const std::vector<T>&, const T&);
template<class T> T question(const std::string&, const std::valarray<T>&, const T&);





/* ---------- function definitions ---------- */





/* ---------- writestuff ---------- */
/// write string to std::cout
void writestuff(const std::string& stuff)
{
	std::cout << stuff;
#ifdef USE_MESSAGES
	messages << stuff;
#endif /* USE_MESSAGES */
}





/* ---------- writestufff ---------- */
/// write string to std::cout and add an endline character
void writestufff(const std::string& stuff)
{
	writestuff(stuff);
	std::cout << std::endl;
#ifdef USE_MESSAGES
	messages << std::endl;
#endif /* USE_MESSAGES */
}





/* ---------- message ---------- */
/// an alias for writestufff
void message(const std::string& message)
{
	writestufff(message);
}





/* ---------- debug ---------- */
/// write debugging information to std::clog
void inline debug(const std::string& function, const std::string& message)
{
#if DEBUG > 0
	std::clog << "DEBUG::" + function + ":\n\t" + message << std::endl;
#ifdef USE_MESSAGES
	messages << "DEBUG::" + function + ":\n\t" + message << std::endl;
#endif /* USE_MESSAGES */
#endif /* DEBUG */
}





/* ---------- debug ---------- */
/// write debugging information to std::clog
void inline debug(const std::string& message)
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





/* ---------- debugn ---------- */
void inline debugn(const std::string& message)
{
#if DEBUG > 0
	std::string s = message;
	std::string::size_type pb = 0;
	while ( pb != s.size() )
	{
		if ( s[pb] == '\n' )
		{
			s.replace(pb, 1, "\n\t");
		}
		pb++;
	}
	std::clog << s;
#ifdef USE_MESSAGES
	messages << s;
#endif /* USE_MESSAGES */
#endif /* DEBUG */
}





/* ---------- debugn ---------- */
void inline debugn(const std::size_t& num, const std::string& message)
{
#if DEBUG > 0
	std::string my_message = "";
	for (std::size_t i=0; i<num; i++)
	{
		my_message += message;
	}
	debugn(my_message);
#endif /* DEBUG */
}





/* ---------- debugnr ---------- */
void inline debugnr(const std::string& message)
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





/* ---------- debugnr ---------- */
void inline debugnr(const std::size_t& num, const std::string& message)
{
#if DEBUG > 0
	std::string my_message = "";
	for (std::size_t i=0; i<num; i++)
	{
		my_message += message;
	}
	debugnr(my_message);
#endif /* DEBUG */
}





/* ---------- counter ---------- */
/// write a counter to std::cerr
void inline counter(const std::string& s, const std::size_t& val)
{
	std::cerr << s << "[";
	std::cerr.width(7);
	std::cerr << val;
	std::cerr << "]" << "\b\b\b\b\b\b\b\b\b";
	for (std::size_t i=0; i<s.size(); i++)
	{
		std::cerr << "\b";
	}
}





/* ---------- counter ---------- */
/// write progressive percentage to std::cerr
void inline counter(const std::string& s, const std::size_t& sz, const std::size_t& val)
{
	static std::size_t counter_percentage = 0;
	if ( val == 0 || val == 1 )
	{
		counter_percentage = 99999;
	}
    std::size_t np = ( 100 * val ) / sz;
	if ( np != counter_percentage )
	{
		counter_percentage = np;
		std::cerr << s << " ";
		std::cerr.width(3);
		std::cerr << counter_percentage;
		std::cerr << "%" << "\b\b\b\b\b";
		for (std::size_t i=0; i<s.size(); i++)
		{
			std::cerr << "\b";
		}
	}
}





/* ---------- yesno ---------- */
/// get a "yes" or "no" response from std::cin
bool yesno(const std::string& question, const bool& init_value)
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





/* ---------- question ---------- */
/// get a response from a question via std::cin
template<class T> T question(const std::string& question, const T& init_value)
{
	writestuff(question + ": [" + std::to_string(init_value) + "] ");
	
	std::string response = "";
	std::getline(std::cin, response);
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
		writestufff(std::to_string(ret));
#endif /* USE_MESSAGES */
		return ret;
	}
#ifdef USE_MESSAGES
	writestufff(std::to_string(init_value));
#endif /* USE_MESSAGES */
	return init_value;
}





/* ---------- question ---------- */
/// get a response from a question via std::cin
std::string question(const std::string& question, const std::string& init_value)
{
    writestuff(question + ": [" + init_value + "] ");
    
    std::string response = "";
    std::getline(std::cin, response);
    if ( response != std::string("") &&
        response != std::string(" ") &&
        response != std::string("\n") &&
        response != std::string("\0") &&
        response != std::string("\t") )
    {
        std::string ret;
        std::istringstream word(response);
        word >> ret;
#ifdef USE_MESSAGES
        writestufff(ret);
#endif /* USE_MESSAGES */
        return ret;
    }
#ifdef USE_MESSAGES
    writestufff(init_value);
#endif /* USE_MESSAGES */
    return init_value;
}





/* ---------- question ---------- */
/// get a response from a question with specified options via std::cin
template<class T> T question(const std::string& question, const std::vector<T>& list, const T& init_value)
{
	for (std::size_t i=0; i<list.size(); i++)
	{
		if ( i < 9 )
		{
			writestuff("  ");
		}
		else if ( i < 99 )
		{
			writestuff(" ");
		}
		writestufff(std::to_string(i + 1) + ". " + std::to_string(list[i]));
	}
	writestuff(question + ": [" + std::to_string(init_value) + "] ");
	
	std::string response = "";
	std::getline(std::cin,response);
	for (std::size_t i=0; i<list.size(); i++)
	{
		if ( response == std::to_string(i + 1) || response == std::string(list[i]) )
		{
#ifdef USE_MESSAGES
			writestufff(std::to_string(list[i]));
#endif /* USE_MESSAGES */
			return list[i];
		}
	}
#ifdef USE_MESSAGES
	writestufff(std::to_string(init_value));
#endif /* USE_MESSAGES */
	return init_value;
}





/* ---------- question ---------- */
/// get a response from a question with specified options via std::cin
std::string question(const std::string& question, const std::vector<std::string>& list, const std::string& init_value)
{
    for (std::size_t i=0; i<list.size(); i++)
    {
        if ( i < 9 )
        {
            writestuff("  ");
        }
        else if ( i < 99 )
        {
            writestuff(" ");
        }
        writestufff(std::to_string(i + 1) + ". " + list[i]);
    }
    writestuff(question + ": [" + init_value + "] ");
    
    std::string response = "";
    std::getline(std::cin,response);
    for (std::size_t i=0; i<list.size(); i++)
    {
        if ( response == std::to_string(i + 1) || response == std::string(list[i]) )
        {
#ifdef USE_MESSAGES
            writestufff(list[i]);
#endif /* USE_MESSAGES */
            return list[i];
        }
    }
#ifdef USE_MESSAGES
    writestufff(init_value);
#endif /* USE_MESSAGES */
    return init_value;
}





/* ---------- question ---------- */
template<class T> T question(const std::string& tquestion, const std::valarray<T>& list, const T& init_value)
{
	std::vector<T> vlist(list.size());
	std::copy(&list[0], &list[list.size()], &vlist[0]);
	return question(tquestion, vlist, init_value);
}





#endif /* _TEXT_FRONT_ */
