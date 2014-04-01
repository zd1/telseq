/*
 * Util.h
 *
 *  Created on: Mar 20, 2014
 *      Author: zd1
 */

#ifndef UTIL_H_
#define UTIL_H_


#include <iostream>
#include <fstream>
#include <istream>
#include <sstream>
#include <set>

typedef std::pair<unsigned int, unsigned int> range;
typedef std::set<range> set_type;


//
// template
//
template <typename T> std::string NumberToString (T Number ){
     std::ostringstream ss;
     ss << Number;
     return ss.str();
}


struct compare {
	bool operator()(range l, range r) const {
		return l.second <= r.first;
	}
};

std::ifstream::pos_type getFilesize(const std::string& filename);
std::istream* createReader(const std::string& filename);
std::ostream* createWriter(const std::string& filename);

//std::set<range>::iterator searchRange(std::set<range>::iterator startit, std::set<range>::iterator endit, range r);
//std::map< std::string, std::set<range> > readBedAsMap(std::string filename);
std::vector<range>::iterator searchRange(std::vector<range>::iterator startit, std::vector<range>::iterator endit, range r);
std::map< std::string, std::vector<range> > readBedAsVector(std::string filename);

inline static std::string refID2Name(int refid)
{	
    if(refid < 22){
		return NumberToString(++refid);    	
    }else if (refid == 22){
    	return "X";
    }else if (refid == 23){
    	return "Y";
    }else{
    	return "-1";
    }
}


#endif /* UTIL_H_ */
