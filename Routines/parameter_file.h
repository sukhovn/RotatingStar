#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <iomanip> 

#ifndef PARAMETER_FILE_H
#define PARAMETER_FILE_H

class Parameter_File{
	public:
		int num;
		int sz;
		std::vector<std::string> fields;
		std::vector<unsigned int> places;
		std::vector<double> values;
		std::string pfile_name;
		
	public:
		//This Parameter_File constructor takes a string with parameter description eg
		//Parameter name: parameter_place\n
		//Where parameter_place is an index value in values array
		Parameter_File(const std::string &input);

		//The same constructor, but takes a pointer to the string with its size
		Parameter_File(const char *input, int input_size);

		//This Parameter_File constructor takes a file with parameter description eg
		//Parameter name: parameter_place
		//Where parameter_place is an index value in values array
		Parameter_File(const char *file_name);

		//This function prints out name from parameter description
		std::string name(void) const{return pfile_name;}

		//This function reads parameter files
		//Parameter name: value
		//And places the value in the corresponding spot in values array 
		void read_file(const std::string &file_name);

		//This function prints out parameters in parameter file
		void print_file(const std::string &file_name) const;

		//This function returns values array size
		int size(void) const{return sz;}

		//This function gives out a constant value out of value array
		const double operator[](unsigned int i) const{
			if(i < sz) return values[i];
			std::cout << "Index out of range\n";
			throw("Index out of range");
		}

		//This function provides acces to the value array
		double& operator[](unsigned int i){
			if(i < sz) return values[i];
			std::cout << "Index out of range\n";
			throw("Index out of range");
		}
};

#endif