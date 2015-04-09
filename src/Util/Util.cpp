/*
 * Util.cpp
 *
 *  Created on: Mar 20, 2014
 *      Author: zd1
 */

#include <iostream>
#include <fstream>
#include <istream>
#include <sstream>
#include <algorithm>
#include <set>
#include <map>
#include <vector>
#include "Util.h"

// Reverse complement a sequence
std::string reverseComplement(const std::string& seq)
{
    std::string out(seq.length(), 'A');
    size_t last_pos = seq.length() - 1;
    for(int i = last_pos; i >= 0; --i)
    {
        out[last_pos - i] = complement(seq[i]);
    }
    return out;
}

std::istream* createReader(const std::string& filename)
{
	std::ifstream* pReader = new std::ifstream(filename.c_str(), std::ifstream::in);
	return pReader;
}

std::ostream* createWriter(const std::string& filename)
{
	std::ofstream* pWriter = new std::ofstream(filename.c_str(), std::ifstream::out);
	return pWriter;
}

// Returns the size of the file. Code from stackoverflow.
std::ifstream::pos_type getFilesize(const std::string& filename)
{
    std::ifstream in(filename.c_str(), std::ifstream::in | std::ifstream::binary);
    in.seekg(0, std::ifstream::end);
    return in.tellg();
}
