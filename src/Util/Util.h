//
// Util - Common data structures and functions
//
#ifndef UTIL_H
#define UTIL_H

#include <vector>
#include <list>
#include <string>
#include <istream>
#include <fstream>
#include <sstream>
#include <stdint.h>
#include <iostream>

typedef std::vector<std::string> StringVector;

//
// Functions
//
std::string stripDirectories(const std::string& filename);
std::ifstream::pos_type getFilesize(const std::string& filename);

// Wrapper function for opening a reader of compressed or uncompressed file
std::istream* createReader(const std::string& filename, std::ios_base::openmode mode = std::ios_base::in);
std::ostream* createWriter(const std::string& filename, 
                           std::ios_base::openmode mode = std::ios_base::out);

void assertFileOpen(std::ifstream& fh, const std::string& fn);
void assertFileOpen(std::ofstream& fh, const std::string& fn);

// Number string conversions

template <typename T>
std::string NumberToString ( T Number )
{
	std::stringstream ss;
	ss << Number;
	return ss.str();
}

template <typename T>
T StringToNumber ( const std::string &Text )//Text not by const reference so that the function can be used with a
{                               //character array as argument
	std::stringstream ss(Text);
	T result;
	return ss >> result ? result : 0;
}


#endif
