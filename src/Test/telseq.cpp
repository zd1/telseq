/*
 * Telseq.cpp
 *
 *  Created on: 22 Aug 2013
 *      Author: zd1
 */

#include <string>
#include <iostream>
#include "scan.h"
#include "merge.h"

#define PROGRAM_BIN "telseq"
#define AUTHOR "Zhihao Ding"

static const char *TELSEQ_VERSION_MESSAGE =
"TelSeq Version " PACKAGE_VERSION "\n"
"Written by " AUTHOR ".\n"
"Copyright 2013 Wellcome Trust Sanger Institute\n" ;

static const char *TELSEQ_USAGE_MESSAGE =
"Program: " PACKAGE_NAME "\n"
"Version: " PACKAGE_VERSION "\n"
"Contact: " AUTHOR " [" PACKAGE_BUGREPORT "]\n"
"Usage: " PROGRAM_BIN " [OPTION] <in.1.bam> <in.2.bam> <...> \n\n"
"Scan BAM and estimate telomere length. \n"
"   <in.bam>                 one or more BAM files to be analyzed. This can also be from a pipe, which should contain 1 BAM path per row.\n"
"   -f, --bamlist=STR        a file that contains a list of file paths of BAMs. It should have only one column, \n"
"                            with each row a BAM file path. -f has higher priority than <in.bam>. When specified, <in.bam> are ignored.\n"
"   -o, --output_dir=STR     output directory for results\n"
"   -H                       remove header line (by default output header line)\n"
"   -k                       threshold of the amount of TTAGGG/CCCTAA repeats in read for a read to be considered telomeric. default = 7.\n"
"   --help                   display this help and exit\n"

"\nReport bugs to " PACKAGE_BUGREPORT "\n\n";

int main(int argc, char** argv)
{
    if(argc <= 1)
    {
        std::cout << TELSEQ_USAGE_MESSAGE;
        return 0;
    }
    else
    {
        std::string command(argv[1]);
        if(command == "help" || command == "--help")
        {
            std::cout << TELSEQ_USAGE_MESSAGE;
            return 0;
        }
        else if(command == "version" || command == "--version")
        {
            std::cout << TELSEQ_VERSION_MESSAGE;
            return 0;
        }

        if(command == "scan")
        {
        	scanMain(argc - 1, argv + 1);
        }

        else if(command == "merge")
        {
        	mergeMain(argc - 1, argv + 1);
        }
        else
        {
            std::cerr << "Unrecognized command: " << command << "\n";
            return 1;
        }
    }

    return 0;
}
