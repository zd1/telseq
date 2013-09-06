/*
 * Bamscan.cpp
 *
 *  Created on: 22 Aug 2013
 *      Author: zd1
 */


#include <iostream>
#include <fstream>


#define PACKAGE_NAME "TelSeq"
#define PACKAGE_VERSION "0.0.1"
#define PACKAGE_BUGREPORT "zd1@sanger.ac.uk"

//
// Getopt
//
#define SUBPROGRAM "SCAN"
static const char *BAMSCAN_VERSION_MESSAGE =
SUBPROGRAM " Version " PACKAGE_VERSION "\n"
"Written by Zhihao Ding.\n"
"\n"
"Copyright 2013 Wellcome Trust Sanger Institute\n";

static const char *BAMSCAN_USAGE_MESSAGE =
"Usage: " PACKAGE_NAME " " SUBPROGRAM " [OPTION] ...\n"
"Scan BAMs \n"
"\n"
"      --help                           display this help and exit\n"
"      -v, --verbose                    display verbose output\n"
"\nInput/Output options:\n"
"      -o, --out=FILE                   write the reads to FILE (default: stdout)\n"
"\nReport bugs to " PACKAGE_BUGREPORT "\n\n";

enum QualityScaling
{
    QS_UNDEFINED,
    QS_NONE,
    QS_SANGER,
    QS_PHRED64
};

namespace opt
{
    static unsigned int verbose;
    static std::string outFile;
    static unsigned int qualityTrim = 0;
    static unsigned int hardClip = 0;
    static int qualityFilter = -1;
    static unsigned int peMode = 0;
    static double sampleFreq = 1.0f;

    static bool bDiscardAmbiguous = true;
    static QualityScaling qualityScale = QS_SANGER;

    static bool bFilterGC = false;
    static bool bDustFilter = false;
    static double dustThreshold = 4.0f;
    static std::string suffix;
    static double minGC = 0.0f;
    static double maxGC = 1.0;
    static bool bIlluminaScaling = false;
    static bool bDisablePrimerCheck = false;
    static std::string orphanFile;
    static std::string adapterF;  // adapter sequence forward
    static std::string adapterR; // adapter sequence reverse
}

static const char* shortopts = "o:q:m:h:p:r:c:s:f:vi";

enum { OPT_HELP = 1, OPT_VERSION, OPT_PERMUTE, OPT_QSCALE, OPT_MINGC, OPT_MAXGC,
       OPT_DUST, OPT_DUST_THRESHOLD, OPT_SUFFIX, OPT_PHRED64, OPT_OUTPUTORPHANS, OPT_DISABLE_PRIMER };
enum {no_argument=0, required_argument=1};

static const struct option longopts[] = {
    { "verbose",                no_argument,       NULL, 'v' },
    { "out",                    required_argument, NULL, 'o' },
    { "quality-trim",           required_argument, NULL, 'q' },
    { "quality-filter",         required_argument, NULL, 'f' },
    { "pe-mode",                required_argument, NULL, 'p' },
    { "hard-clip",              required_argument, NULL, 'h' },
    { "min-length",             required_argument, NULL, 'm' },
    { "sample",                 required_argument, NULL, 's' },
    { "remove-adapter-fwd",     required_argument, NULL, 'r' },
    { "remove-adapter-rev",     required_argument, NULL, 'c' },
    { "dust",                   no_argument,       NULL, OPT_DUST},
    { "dust-threshold",         required_argument, NULL, OPT_DUST_THRESHOLD },
    { "suffix",                 required_argument, NULL, OPT_SUFFIX },
    { "phred64",                no_argument,       NULL, OPT_PHRED64 },
    { "pe-orphans",             required_argument, NULL, OPT_OUTPUTORPHANS },
    { "min-gc",                 required_argument, NULL, OPT_MINGC},
    { "max-gc",                 required_argument, NULL, OPT_MAXGC},
    { "help",                   no_argument,       NULL, OPT_HELP },
    { "version",                no_argument,       NULL, OPT_VERSION },
    { "permute-ambiguous",      no_argument,       NULL, OPT_PERMUTE },
    { "no-primer-check",        no_argument,       NULL, OPT_DISABLE_PRIMER },
    { NULL, 0, NULL, 0 }
};

//
// Main
//
int preprocessMain(int argc, char** argv)
{
	parsePreprocessOptions(argc, argv);
    return 0;
}


double calcGC(const std::string& seq)
{
    double num_gc = 0.0f;
    double num_total = 0.0f;
    for(size_t i = 0; i < seq.size(); ++i)
    {
        if(seq[i] == 'C' || seq[i] == 'G')
            ++num_gc;
        ++num_total;
    }
    return num_gc / num_total;
}

//
// Handle command line arguments
//
void parsePreprocessOptions(int argc, char** argv)
{
    bool die = false;
    for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;)
    {
        std::istringstream arg(optarg != NULL ? optarg : "");
        switch (c)
        {
            case 'o': arg >> opt::outFile; break;
            case 'q': arg >> opt::qualityTrim; break;
            case 'f': arg >> opt::qualityFilter; break;
            case 'i': arg >> opt::bIlluminaScaling; break;
            case 'h': arg >> opt::hardClip; break;
            case 'p': arg >> opt::peMode; break;
            case 'r': arg >> opt::adapterF; break;
            case 'c': arg >> opt::adapterR; break;
            case 's': arg >> opt::sampleFreq; break;
            case '?': die = true; break;
            case 'v': opt::verbose++; break;
            case OPT_DUST_THRESHOLD: arg >> opt::dustThreshold; opt::bDustFilter = true; break;
            case OPT_SUFFIX: arg >> opt::suffix; break;
            case OPT_MINGC: arg >> opt::minGC; opt::bFilterGC = true; break;
            case OPT_OUTPUTORPHANS: arg >> opt::orphanFile; break;
            case OPT_MAXGC: arg >> opt::maxGC; opt::bFilterGC = true; break;
            case OPT_PHRED64: opt::qualityScale = QS_PHRED64; break;
            case OPT_PERMUTE: opt::bDiscardAmbiguous = false; break;
            case OPT_DUST: opt::bDustFilter = true; break;
            case OPT_DISABLE_PRIMER: opt::bDisablePrimerCheck = true; break;
            case OPT_HELP:
                std::cout << BAMSCAN_USAGE_MESSAGE;
                exit(EXIT_SUCCESS);
            case OPT_VERSION:
                std::cout << BAMSCAN_VERSION_MESSAGE;
                exit(EXIT_SUCCESS);
        }
    }

    if (argc - optind < 1)
    {
        std::cerr << SUBPROGRAM ": missing arguments\n";
        die = true;
    }

    if (die)
    {
        std::cout << "\n" << BAMSCAN_USAGE_MESSAGE;
        exit(EXIT_FAILURE);
    }

    if(opt::peMode > 2)
    {
        std::cerr << SUBPROGRAM ": error pe-mode must be 0,1 or 2 (found: " << opt::peMode << ")\n";
        exit(EXIT_FAILURE);
    }

    if(opt::adapterF.empty() != opt::adapterR.empty())
    {
        std::cerr << SUBPROGRAM ": Forward and Reverse sequence is necessary to perform adapter removal.\n";
        exit(EXIT_FAILURE);
    }

    if(!opt::outFile.empty() && opt::outFile == opt::orphanFile)
    {
        std::cerr << SUBPROGRAM ": Output file and orphan file must be different\n";
        exit(EXIT_FAILURE);
    }
}


