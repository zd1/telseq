/*
 * merge.cpp
 *
 *  Created on: 3 Sep 2013
 *      Author: zd1
 */

#include <iostream>
#include <istream>
#include <iterator>
#include <fstream>
#include <algorithm>
#include "merge.h"
#include "Timer.h"
#include "Util.h"
#include "dirent.h"
#include "assert.h"
#include "scan.h"

//
// Getopt
//
#define SUBPROGRAM "merge"

static const char *MERGE_VERSION_MESSAGE =
SUBPROGRAM " Version " PACKAGE_VERSION "\n"
"Written by Zhihao Ding.\n"
"\n"
"Copyright 2013 Wellcome Trust Sanger Institute\n";

static const char *MERGE_USAGE_MESSAGE =
"Usage: " PACKAGE_NAME " " SUBPROGRAM " [OPTION]\n"
"Merge the results generated from the previous scan step. \n"
"   -o, --output_dir=STR     output directory for results\n"
"   --help                   display this help and exit\n"

"\nReport bugs to " PACKAGE_BUGREPORT "\n\n";

namespace opt
{
    static std::string outputdir = ".";
}

static const char* shortopts = "o:";

enum { OPT_HELP = 1, OPT_VERSION };

static const struct option longopts[] = {
    { "help",               no_argument,       NULL, OPT_HELP },
    { "version",            no_argument,       NULL, OPT_VERSION },
    { "output-dir",		optional_argument, NULL, 'o' },
    { NULL, 0, NULL, 0 }
};

int mergeMain(int argc, char** argv)
{
	Timer* pTimer = new Timer("merge results");
	parseMergeOptions(argc, argv);
	mergeBam();
    delete pTimer;
    return 0;
}

int mergeBam()
{
	StringVector bamscans;
	DIR *dir;
	struct dirent *ent;
	if ((dir = opendir(opt::outputdir.c_str())) != NULL) {
	  /* print all the files and directories within directory */
	  while ((ent = readdir (dir)) != NULL) {
		  std::string rf = ent->d_name;
		  size_t p1 = rf.find(ScanParameters::SCAN_FILE_SUFFIX);
		  if(p1!=std::string::npos){
			  bamscans.push_back(rf);
		  }
	  }
	  closedir(dir);
	} else {
	  /* could not open directory */
		std::cerr << "could not open specified output directory" << opt::outputdir << "\n";
		exit(EXIT_FAILURE);
	}

	std::cout << "Found " << bamscans.size() << " result files \n";

	std::vector<ScanResults> allresults;
	for(size_t i=0; i < bamscans.size();i++){
		ScanResults sr = parseResultFile(bamscans[i]);
		allresults.push_back(sr);
	}

	std::cout << "Parsed " << bamscans.size() << " result files \n";
	std::string filepath = opt::outputdir+"/summary";
	printsummary(allresults, filepath);
	std::cout << "Wrote summary result file " << bamscans.size() << "\n";
	return 0;
}

int printsummary(std::vector<ScanResults> &results, std::string filepath){

	std::ostream* pWriter = createWriter(filepath);
	StringVector headers;
	headers.push_back(ScanParameters::LABEL_SAMPLE);
	headers.push_back(ScanParameters::LABEL_TOTAL);
	headers.push_back(ScanParameters::LABEL_MAPPED);
	headers.push_back(ScanParameters::LABEL_DUP);
	headers.push_back(ScanParameters::LABEL_LEN);

	for(int i=0;i < ScanParameters::TEL_MOTIF_N;i++){
		std::string h = ScanParameters::LABEL_TEL + NumberToString(i);
		headers.push_back(h);
	}
	for(int i=0;i<ScanParameters::GC_BIN_N;i++){
		std::string h = ScanParameters::LABEL_GC + NumberToString(i);
		headers.push_back(h);
	}

	for(size_t h=0; h<headers.size();h++){
		*pWriter << headers[h] << ScanParameters::FIELD_SEP;
	}
	*pWriter << "\n";

	for(size_t i=0; i< results.size();i++){
		ScanResults result=results[i];
		*pWriter << result.sample << ScanParameters::FIELD_SEP;
		*pWriter << result.numTotal << ScanParameters::FIELD_SEP;
		*pWriter << result.numMapped << ScanParameters::FIELD_SEP;
		*pWriter << result.numDuplicates << ScanParameters::FIELD_SEP;
		*pWriter << result.telLenEstimate << ScanParameters::FIELD_SEP;

		for (std::size_t j = 0, max = result.telcounts.size(); j != max; ++j){
			*pWriter << result.telcounts[j] << ScanParameters::FIELD_SEP;
		}
		for (std::size_t k = 0, max = result.gccounts.size(); k != max; ++k){
			*pWriter << result.gccounts[k] << ScanParameters::FIELD_SEP;
		}
		*pWriter << "\n";
	}

	delete pWriter;
	return 0;
}

// read result file
ScanResults parseResultFile(std::string filename)
{

    size_t filesize = getFilesize(filename);
    if(filesize == 0)
    {
        std::cerr << SUBPROGRAM ": result file is empty\n";
        exit(EXIT_FAILURE);
    }

    std::istream* pReader = createReader(filename);
    std::string line;
    ScanResults results;

    while(getline(*pReader, line))
    {
    	if(line.empty()){
    		continue;
    	}

    	StringVector fields = split(line, ':');
    	assert(fields.size()==2);
    	std::string label=fields[0];
    	std::string val = fields[1];

    	bool found = false;

    	if(label == ScanParameters::LABEL_SAMPLE){
    		results.sample = val;
    		found = true;
    	}else if(label == ScanParameters::LABEL_BAM){
    		results.bam = val;
    		found = true;
    	}else if(label == ScanParameters::LABEL_LEN){
    		results.telLenEstimate = StringToNumber<double>(val);
    		found = true;
    	}else if(label == ScanParameters::LABEL_TOTAL){
    		results.numTotal = StringToNumber<int>(val);
    		found = true;
    	}else if(label == ScanParameters::LABEL_MAPPED){
    		results.numMapped = StringToNumber<int>(val);
    		found = true;
    	}else if(label == ScanParameters::LABEL_DUP){
    		results.numDuplicates = StringToNumber<int>(val);
    		found = true;
    	}else if(label.find(ScanParameters::LABEL_TEL,0)){
    		std::string idx_str = label.substr (ScanParameters::LABEL_TEL.size());
    		int idx = StringToNumber<int>(idx_str);
    		results.telcounts[idx]=StringToNumber<int>(val);
    		found = true;
    	}else if(label.find(ScanParameters::LABEL_GC,0)){
    		std::string idx_str = label.substr (ScanParameters::LABEL_GC.size());
    		int idx = StringToNumber<int>(idx_str);
    		results.gccounts[idx]= StringToNumber<int>(val);
    		found = true;
    	}

    	if(!found){
    		std::cerr << "field not found. field={"<< label << "}\n";
    		exit(EXIT_FAILURE);
    	}

    }

    delete pReader;
    return results;
}


//
// Handle command line arguments
//
void parseMergeOptions(int argc, char** argv)
{

    bool die = false;
    for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;)
    {
        std::istringstream arg(optarg != NULL ? optarg : "");
        switch (c)
        {   case 'o':
        		arg >> opt::outputdir; break;
            case OPT_HELP:
                std::cout << MERGE_USAGE_MESSAGE;
                exit(EXIT_SUCCESS);
            case OPT_VERSION:
                std::cout << MERGE_VERSION_MESSAGE;
                exit(EXIT_SUCCESS);
        }
    }


    if (argc - optind > 1)
    {
        std::cerr << SUBPROGRAM ": too many arguments\n";
        die = true;
    }

    if (die)
    {
        std::cout << "\n" << MERGE_USAGE_MESSAGE;
        exit(EXIT_FAILURE);
    }

}



