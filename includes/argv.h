/*
 argv.h
 
 Functions for getting any command line inputs and options.
 
 Rado Faletic
 Department of Physics
 Faculty of Science
 Australian National University  ACT  0200
 Australia
 
 Rado.Faletic@anu.edu.au
 8th September 2004
 */


#include <string>
#include <vector>


class args
{
private:
	std::string var_;
	std::string val_;
public:
	args() : var_(""), val_("") { };
	std::string var() const { return this->var_; };
	std::string& var() { return this->var_; };
	bool var(const std::string& cmp) const { return ( this->var_ == cmp ); };
	std::string val() const { return this->val_; };
	std::string& val() { return this->val_; };
	bool val(const std::string& cmp) const { return ( this->val_ == cmp ); };
};


std::vector<args> get_args(const int& argc, char* argv[])
{
	std::vector<args> output(argc-1);
	for (size_t i=0; i<output.size(); i++)
	{
		std::string istring(argv[i+1]);
		if ( istring.substr(0,2) != std::string("--") )
		{
			continue;
		}
		std::string::size_type f1 = istring.find("--") + 2;
		std::string::size_type f2 = istring.find("=");
		if ( f2 <= f1 )
		{
			continue;
		}
		output[i].var() = istring.substr(f1, f2-f1);
		f2++;
		output[i].val() = istring.substr(f2, istring.size()-f2);
	}
	for (int i=output.size()-1; i>=0; i--)
	{
		if ( !output[i].var().size() || !output[i].val().size() )
		{
			output.erase(output.begin()+i);
		}
	}
	return output;
}
