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


std::vector<range>::iterator searchRange(std::vector<range>::iterator startit, std::vector<range>::iterator endit, range r){

	typedef std::vector<range> set_type;
	set_type::iterator it_l = std::lower_bound(startit, endit, r, compare());
	set_type::iterator it_u = std::upper_bound(startit, endit, r, compare());

//	std::cout << "      lower: " << it_l->first << ' ' << it_l->second << '\n';
// 	std::cout << "      upper: " << it_u->first << ' ' << it_u->second << '\n';
// 	std::cout << "      it_l == end: " << (it_l == endit) << "\n";
//	std::cout << "      it_u == end: " << (it_u == endit) << "\n";
// 	std::cout << "      it_l == it_u: " << (it_l == it_u) << "\n";

	if (it_l == it_u) { // if not found, return the end
		return endit;
	}else{
		return it_l; // if found, return the lower bound.
	}
}


std::map< std::string, std::vector<range> > readBedAsVector(std::string filename){
    
	std::map< std::string, std::vector<range> > bed;
   	std::ifstream inputFile(filename);
    std::string line;
    
    while (getline(inputFile, line))
    {
    	std::istringstream ss(line);
        std::string chrm;
        int s,e;
        ss >> chrm >> s >> e;
        bed[chrm].push_back({s,e});
    }
    return bed;
}



std::map< std::string, std::set<range> > readBedAsMap(std::string filename){

	std::map< std::string, std::set<range> > bed;
   	std::ifstream inputFile(filename);
    std::string line;

    while (getline(inputFile, line))
    {
    	std::istringstream ss(line);
        std::string chrm;
        int s,e;
        ss >> chrm >> s >> e;
        bed[chrm].insert({s,e});
    }

    return bed;
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






