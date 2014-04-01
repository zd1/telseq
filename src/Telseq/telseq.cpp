
#include <iostream>
#include <istream>
#include <iterator>
#include <fstream>
#include <algorithm>
#include <vector>
#include <istream>
#include <sstream>
#include <set>
#include <map>
#include <unistd.h>

#include "telseq.h"
#include "Timer.h"
#include "math.h"
#include "api/BamReader.h"
#include "api/BamWriter.h"
#include "Util.h"
#include "prettyprint.h"

//
// Getopt
//

#define PROGRAM_BIN "telseq"
#define AUTHOR "Zhihao Ding"

static const char *TELSEQ_VERSION_MESSAGE =
"TelSeq Version " PACKAGE_VERSION "\n"
"Written by " AUTHOR ".\n"
"Copyright 2013 Wellcome Trust Sanger Institute\n" ;

static const char *TELSEQ_USAGE_MESSAGE =
"Program: " PACKAGE_NAME "\n"
"Version: " PACKAGE_VERSION "\n"
"Contact: " AUTHOR " [" PACKAGE_BUGREPORT "]\n\n"
"Usage: " PROGRAM_BIN " [OPTION] <in.1.bam> <in.2.bam> <...> \n"
"Scan BAM and estimate telomere length. \n"
"   <in.bam>                 one or more BAM files to be analysed. File names can also be passed from a pipe, \n "
"                            with each row containing 1 BAM path.\n"
"   -f, --bamlist=STR        a file that contains a list of file paths of BAMs. It should contain only one column, \n"
"                            with each row a BAM file path. -f has higher priority than <in.bam>. When specified, \n"
"                            <in.bam> are ignored.\n"
"   -o, --output_dir=STR     output file for results. Ignored when input is from stdin, in which case output will be stdout. \n"
"   -H                       remove header line, which is printed by default.\n"
"   -h                       print the header line only. The text can be used to attach to result files, useful\n"
"                            when the headers of the result files are suppressed. \n"
"   -m                       merge read groups by taking weighted average, by the total number of reads in each read group,\n"
"                            across all read groups for a sample. Default is to output each readgroup separately.\n"
"   -u                       ignore read groups. Treat all reads in BAM as if they were from a same read group.\n"
"   -k                       threshold of the amount of TTAGGG/CCCTAA repeats in read for a read to be considered telomeric. default = 7.\n"
"\nTesting functions\n------------\n"
"   -z                       use user specified pattern for searching.\n"
"   -e, --exomebed=STR       specifiy exome regions in the BED format. These regions will be excluded when \n"
"                            calculating estimates, thus an estimate of the absolute length can be obtained, \n"
"                            same as is for the whole genome sequence BAMs.\n"
"   --help                   display this help and exit\n"

"\nReport bugs to " PACKAGE_BUGREPORT "\n\n";

namespace opt
{
    static StringVector bamlist;
    static std::string outputfile = "";
    static std::string exomebedfile = "";
    static std::map< std::string, std::vector<range> > exomebed;
    static bool writerheader = true;
    static bool mergerg = false;
    static bool ignorerg = false;
    static int tel_k= ScanParameters::TEL_MOTIF_CUTOFF;
    static std::string unknown = "UNKNOWN";
    static std::string PATTERN=ScanParameters::PATTERN;
    static std::string PATTERN_REV=ScanParameters::PATTERN_REVCOMP;

}

static const char* shortopts = "f:o:k:z:e:Hhvmu";

enum { OPT_HELP = 1, OPT_VERSION };

static const struct option longopts[] = {
	{ "bamlist",		optional_argument, NULL, 'f' },
    { "output-dir",		optional_argument, NULL, 'o' },
    { "exomebed",		optional_argument, NULL, 'e' },
    { "help",               no_argument,       NULL, OPT_HELP },
    { "version",            no_argument,       NULL, OPT_VERSION },
    { NULL, 0, NULL, 0 }
};



int scanBam()
{

	std::vector< std::map<std::string, ScanResults> > resultlist;
//	std::ostream* pWriter;
//	pWriter = &std::cout;
	bool isExome = opt::exomebedfile.size()==0? false: true;
	
	std::cout << opt::bamlist << "\n";
	
    for(std::size_t i=0; i<opt::bamlist.size(); i++) {

        // storing results for each read group (RG tag). use
        // read group ID as key.
        std::map<std::string, ScanResults> resultmap;

        // store where the overlap was last found in the case of exome seq
    	std::map<std::string, std::vector<range>::iterator> lastfound;
    	std::vector<range>::iterator searchhint;
        
        std::cerr << "Start analysing BAM " << opt::bamlist[i] << "\n";

        // Open the bam files for reading/writing
        BamTools::BamReader* pBamReader = new BamTools::BamReader;

        pBamReader->Open(opt::bamlist[i]);

        // get bam headers
        const BamTools::SamHeader header = pBamReader ->GetHeader();

//        for(BamTools::SamSequenceConstIterator it = header.Sequences.Begin();
//        						it != header.Sequences.End();++it){
//        	std::cout << "Assembly ID:" << it->AssemblyID << ", Name:" << it->Name << std::endl;
//        }
//        exit(0);

        bool rggroups=false;
        
        if(opt::ignorerg){ // ignore read groups
        	std::cerr << "Treat all reads in BAM as if they were from a same sample" << std::endl;
        	ScanResults results;
        	results.sample = opt::unknown;
        	resultmap[opt::unknown]=results;
        }else{


			std::map <std::string, std::string> readgroups;
			std::map <std::string, std::string> readlibs;

			rggroups = header.HasReadGroups();

			if(rggroups){
				for(BamTools::SamReadGroupConstIterator it = header.ReadGroups.Begin();
						it != header.ReadGroups.End();++it){
					readgroups[it->ID]= it->Sample;
					if(it->HasLibrary()){
						readlibs[it->ID] = it -> Library;
					}else{
						readlibs[it->ID] = opt::unknown;
					}
				}
				std::cerr<<"Specified BAM has "<< readgroups.size()<< " read groups" << std::endl;

				for(std::map<std::string, std::string>::iterator it = readgroups.begin(); it != readgroups.end(); ++it){
					ScanResults results;
					std::string rgid = it -> first;
					results.sample = it -> second;
					results.lib = readlibs[rgid];
					resultmap[rgid]=results; //results are identified by RG tag.
				}

			}else{
				std::cerr << "Warning: can't find RG tag in the BAM header" << std::endl;
				std::cerr << "Warning: treat all reads in BAM as if they were from a same sample" << std::endl;
				ScanResults results;
				results.sample = opt::unknown;
				results.lib = opt::unknown;
				resultmap[opt::unknown]=results;
			}
        }

        BamTools::BamAlignment record1;
        bool done = false;
        
        int nprocessed=0; // number of reads analyzed
        int ntotal=0; // number of reads scanned in bam (we skip some reads, see below)
        while(!done)
        {
            ntotal ++;
            done = !pBamReader -> GetNextAlignment(record1);
            std::string tag = opt::unknown;
            if(rggroups){
                
                // skip reads that do not have read group tag
            	if(record1.HasTag("RG")){
					record1.GetTag("RG", tag);
	//            	std::cerr << c << " reads:{" << record1.QueryBases << "} tag:{" << tag << "}\n";
				}else{
					std::cerr << "can't find RG tag for read at position {" << record1.RefID << ":" << record1.Position << "}" << std::endl;
					std::cerr << "skip this read" << std::endl;
					continue;
				}
            }
            
            // skip reads with readgroup not defined in BAM header
            if(resultmap.find(tag) == resultmap.end()){
				std::cerr << "RG tag {" << tag << "} for read at position ";
				std::cerr << "{" << record1.RefID << ":" << record1.Position << "} doesn't exist in BAM header.";
				continue;
            }

            // for exome, exclude reads mapped to the exome regions.
            if(isExome){
				range rg;
				rg.first = record1.Position;
				rg.second = record1.Position + record1.Length;
				std::string chrm =  refID2Name(record1.RefID);

				if(chrm != "-1"){ // check if overlap exome when the read is mapped to chr1-22, X, Y
					// std::cerr << "read: " << chrm << " " << rg << "\n" << std::endl;
					std::map<std::string, std::vector<range> >::iterator chrmit = opt::exomebed.find(chrm);

					if(chrmit == opt::exomebed.end())
					{
						std::cerr<<"chromosome or reference sequence: " << chrm << " is not present in the specified exome bed file." <<std::endl;
						std::cerr<<"please check sequence name encoding, i.e. for chromosome one, is it chr1 or 1" << std::endl;
						resultmap[tag].n_exreadsChrUnmatched +=1;
					}else{

						std::vector<range>::iterator itend = opt::exomebed[chrm].end();
						std::map<std::string, std::vector<range>::iterator>::iterator lastfoundchrmit = lastfound.find(chrm);
						if(lastfoundchrmit == lastfound.end()){ // first entry to this chrm
							lastfound[chrm] = chrmit->second.begin();// start from begining
						}

						// set the hint to where the previous found is
						searchhint = lastfound[chrm];
						std::vector<range>::iterator itsearch = searchRange(searchhint, itend, rg);
						// if found
						if(itsearch != itend){// if found
							searchhint = itsearch;
							resultmap[tag].n_exreadsExcluded +=1;
							lastfound[chrm] = searchhint; // update search hint
							continue;
						}
					}

				}
            }

            resultmap[tag].numTotal +=1;

            if(!record1.IsMapped())
            {
            	resultmap[tag].numMapped += 1;
            }
            if(record1.IsDuplicate()){
            	resultmap[tag].numDuplicates +=1;
            }

            double gc = calcGC(record1.QueryBases);
            int ptn_count = countMotif(record1.QueryBases, opt::PATTERN, opt::PATTERN_REV);

            resultmap[tag].telcounts[ptn_count]+=1;

            if(gc >= ScanParameters::GC_LOWERBOUND && gc <= ScanParameters::GC_UPPERBOUND){
            	// get index for GC bin.
            	int idx = floor((gc-ScanParameters::GC_LOWERBOUND)/ScanParameters::GC_BINSIZE);
            	assert(idx >=0 && idx <= ScanParameters::GC_BIN_N-1);
//            	std::cerr << c << " GC:{"<< gc << "} telcounts:{"<< ptn_count <<"} GC idx{" << idx << "}\n";
            	if(idx > ScanParameters::GC_BIN_N-1){
            		std::cerr << nprocessed << " GC:{"<< gc << "} telcounts:{"<< ptn_count <<"} GC bin index out of bound:" << idx << "\n";
            		exit(EXIT_FAILURE);
            	}
            	resultmap[tag].gccounts[idx]+=1;
            }
            
            if(resultmap[tag].n_exreadsChrUnmatched > 1000){
            	std::cerr<<"too many reads found with unmatched chromosome ID between BAM and exome BED. \n" << std::endl;
            }
            
            nprocessed++;

            if( nprocessed%10000000 == 0){
            	std::cerr << "[scan] processed " << nprocessed << " reads \n" ;
                
            }
            
        }
        
		pBamReader->Close();
        delete pBamReader;
        
        std::cerr << "[scan] total reads in BAM scanned " << ntotal << std::endl;
        std::cerr << "Completed scanning BAM\n";
	    
        resultlist.push_back(resultmap);
    }

    outputresults(resultlist);
    if(isExome){
    	printlog(resultlist);
    }
    return 0;
}

void printlog(std::vector< std::map<std::string, ScanResults> > resultlist){

	for(size_t i =0; i< resultlist.size(); i++){
		auto rmap = resultlist[i];
		for(std::map<std::string, ScanResults>::iterator it= rmap.begin();
				it != rmap.end(); ++it){
			
			std::string rg = it ->first;
			ScanResults result = it -> second;

			std::cerr << "BAM:" << rg << std::endl;
			std::cerr << "	chr ID unmatched reads: " << result.n_exreadsChrUnmatched << std::endl;
			std::cerr << "	exome reads excluded: " << result.n_exreadsExcluded << std::endl;
		}	
	}

}


void printout(std::string rg, ScanResults result, std::ostream* pWriter){

	*pWriter << rg << ScanParameters::FIELD_SEP;
	*pWriter << result.lib << ScanParameters::FIELD_SEP;
	*pWriter << result.sample << ScanParameters::FIELD_SEP;
	*pWriter << result.numTotal << ScanParameters::FIELD_SEP;
	*pWriter << result.numMapped << ScanParameters::FIELD_SEP;
	*pWriter << result.numDuplicates << ScanParameters::FIELD_SEP;

	result.telLenEstimate = calcTelLength(result);
	if(result.telLenEstimate==-1){
		std::cerr << "Telomere length estimate unknown. No read was found with telomeric GC composition.\n";
		*pWriter << opt::unknown << ScanParameters::FIELD_SEP;
	}else if(result.telLenEstimate>1000000){
		std::cerr << "Telomere length estimate unknown. Possibly due to not enough representation of genome.\n";
		*pWriter << opt::unknown << ScanParameters::FIELD_SEP;
	}else if(result.telLenEstimate==0){
		std::cerr << "Telomere length estimate unknown. No read contains more than " << opt::tel_k << " telomere repeats.\n";
		*pWriter << opt::unknown << ScanParameters::FIELD_SEP;
	}
	else{
		*pWriter << result.telLenEstimate << ScanParameters::FIELD_SEP;
	}

	for (std::size_t j = 0, max = result.telcounts.size(); j != max; ++j){
		*pWriter << result.telcounts[j] << ScanParameters::FIELD_SEP;
	}
	for (std::size_t k = 0, max = result.gccounts.size(); k != max; ++k){
		*pWriter << result.gccounts[k] << ScanParameters::FIELD_SEP;
	}
	*pWriter << "\n";
}


int outputresults(std::vector< std::map<std::string, ScanResults> > resultlist){

	std::ostream* pWriter;
	bool tostdout = opt::outputfile.empty() ? true:false;
	if(tostdout){
		pWriter = &std::cout;
	}else{
		pWriter = createWriter(opt::outputfile);
	}

	if(opt::writerheader){
		Headers hd;
		for(size_t h=0; h<hd.headers.size();h++){
			*pWriter << hd.headers[h] << ScanParameters::FIELD_SEP;
		}
		*pWriter << "\n";
	}

	ScanResults mergedrs;
	std::string grpnames = "";

	for(size_t i=0; i < resultlist.size();++i){

		std::map<std::string, ScanResults> resultmap = resultlist[i];

		// if merge read groups, take weighted average for all measures
		bool domg = opt::mergerg && resultmap.size() > 1? true:false;

		for(std::map<std::string, ScanResults>::iterator it= resultmap.begin();
				it != resultmap.end(); ++it){

			std::string rg = it ->first;
			ScanResults result = it -> second;

			if(domg){
				if(grpnames.size()==0){
					grpnames += rg;
				}else{
					grpnames += "|"+rg;
				}

				mergedrs.sample = result.sample;
				mergedrs.numTotal += result.numTotal;
				mergedrs.numMapped += result.numMapped * result.numTotal;
				mergedrs.numDuplicates += result.numDuplicates * result.numTotal;
				mergedrs.telLenEstimate += calcTelLength(result)* result.numTotal;

				for (std::size_t j = 0, max = result.telcounts.size(); j != max; ++j){
					mergedrs.telcounts[j] += result.telcounts[j]* result.numTotal;
				}
				for (std::size_t k = 0, max = result.gccounts.size(); k != max; ++k){
					mergedrs.gccounts[k] += result.gccounts[k]* result.numTotal;
				}
				continue;
			}else{
				printout(rg, result, pWriter);
			}
		}

		//in this case calculate weighted average
		if(domg){

			mergedrs.numMapped /= mergedrs.numTotal;
			mergedrs.numDuplicates /= mergedrs.numTotal;
			mergedrs.telLenEstimate /= mergedrs.numTotal;

			for (std::size_t j = 0, max = mergedrs.telcounts.size(); j != max; ++j){
				mergedrs.telcounts[j] /= mergedrs.numTotal;
			}
			for (std::size_t k = 0, max = mergedrs.gccounts.size(); k != max; ++k){
				mergedrs.gccounts[k] /= mergedrs.numTotal;
			}

			mergedrs.numTotal  /= resultmap.size();

			printout(grpnames, mergedrs, pWriter);
		};
	}

	if(!tostdout){
		delete pWriter;
	}
	return 0;
}

double calcTelLength(ScanResults results){

	float acc = 0;
	for(std::size_t i=opt::tel_k, max = results.telcounts.size(); i !=max; ++i){
		acc += results.telcounts[i];
	}

	float gc_tel = 0;
	for(std::size_t i=0, max = results.gccounts.size(); i !=max; ++i){
		float gc1=ScanParameters::GC_LOWERBOUND + ScanParameters::GC_BINSIZE*i;
		float gc2=ScanParameters::GC_LOWERBOUND + ScanParameters::GC_BINSIZE*(i+1);
		if(gc1 >= ScanParameters::GC_TELOMERIC_LOWERBOUND && gc2 <= ScanParameters::GC_TELOMERIC_UPPERBOUND ){
			gc_tel += results.gccounts[i];
		}
	}

	if(gc_tel == 0){
		return -1;
	}

	return (acc/gc_tel)*float(ScanParameters::GENOME_LENGTH_AT_TEL_GC)/46000;
}



int countMotif(std::string &read, std::string pattern, std::string pattern_revcomp){

	int motifcount1 = 0;
	int motifcount2 = 0;

	size_t p1 = read.find(pattern, 0);
	while(p1 != std::string::npos)
	{
	    p1 = read.find(pattern,p1+pattern.size());
	    motifcount1 += 1;
	}

	size_t p2 = read.find(pattern_revcomp, 0);
	while(p2 != std::string::npos)
	{
	    p2 = read.find(pattern_revcomp,p2+pattern_revcomp.size());
	    motifcount2 += 1;
	}

	return motifcount1 > motifcount2? motifcount1:motifcount2;

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
void parseScanOptions(int argc, char** argv)
{

	std::string bamlistfile =  "";
	std::string rev = "";

	Headers hd;
    for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;)
    {
        std::istringstream arg(optarg != NULL ? optarg : "");
        switch (c)
        {
            case 'f':
            	arg >> bamlistfile; break;
            case 'o':
            	arg >> opt::outputfile; break;
            case 'H':
            	opt::writerheader=false; break;
            case 'm':
                opt::mergerg = true; break;
            case 'u':
                opt::ignorerg = true; break;
            case 'h':
        		for(size_t h=0; h<hd.headers.size();h++){
        			std::cout << hd.headers[h] << ScanParameters::FIELD_SEP;
        		}
        		std::cout << "\n";
        		exit(EXIT_SUCCESS);
            case 'k':
            	arg >> opt::tel_k;
            	if(opt::tel_k < 1 or opt::tel_k > ScanParameters::TEL_MOTIF_N-1){
            		std::cerr << "k out of bound. k must be an integer from 1 to " <<  ScanParameters::TEL_MOTIF_N-1 << "\n";
            		exit(EXIT_FAILURE);
            	}
            	break;
            case 'z':
            	arg >> opt::PATTERN;
            	for(size_t i=0;i<opt::PATTERN.size();i++){
            		rev += complement(opt::PATTERN[i]);
            	}
            	opt::PATTERN_REV =  rev;
            	break;

            case 'e':
            	arg >> opt::exomebedfile;
            	opt::exomebed = readBedAsVector(opt::exomebedfile);
            	std::cout << "loaded "<< opt::exomebed.size() << " exome regions \n"<< std::endl;
//            	std::cout << opt::exomebed << "\n";
            	break;

            case OPT_HELP:
                std::cout << TELSEQ_USAGE_MESSAGE;
                exit(EXIT_SUCCESS);
            case OPT_VERSION:
                std::cout << TELSEQ_VERSION_MESSAGE;
                exit(EXIT_SUCCESS);
        }
    }

    // deal with cases of API usage:
    // telseq a.bam b.bam c.bam ...
    // | telseq
    // telseq -f

    if (argc - optind < 1) // no argument specified
    {
    	// check if it is from pipe
		if(!isatty(fileno(stdin))){
			std::string line;
			while (std::getline(std::cin, line))
			{
				if(line.empty()){
					continue;
				}
				opt::bamlist.push_back(line);
			}
		}else if(bamlistfile.empty() ){ // check if not from a pipe, -f must be spceified
//    		std::cerr << SUBPROGRAM ": No BAM specified. Please specify BAM either directly, by using -f option or piping BAM file path.\n";
    		std::cout << TELSEQ_USAGE_MESSAGE;
    		exit(EXIT_FAILURE);
    	}
    }

    else if (argc - optind >= 1) // if arguments are specified
    {
    	// -f has higher priority, when specified, ignore arguments.
    	if(bamlistfile.empty()){
    		for(int i = optind; i < argc; ++i ){
    		    opt::bamlist.push_back(argv[i]);
    		}
    	}
    }

    // read in bamlist
    if(!bamlistfile.empty()){
        size_t filesize = getFilesize(bamlistfile);
        if(filesize == 0)
        {
            std::cerr << PROGRAM_BIN ": BAMLIST file specified by -f is empty\n";
            exit(EXIT_FAILURE);
        }

        std::istream* pReader = createReader(bamlistfile);
        std::string line;

        while(getline(*pReader, line))
        {
        	if(line.empty()){
        		continue;
        	}
            opt::bamlist.push_back(line);
        }

        size_t bamsize = opt::bamlist.size();
        if(bamsize == 0 ){
            std::cerr << PROGRAM_BIN ": Could find any sample in BAMLIST file specified.\n";
            exit(EXIT_FAILURE);
        }
        delete pReader;
    }

    // check output
//    for(size_t i = 0; i < opt::bamlist.size(); ++i ){
//    	std::cerr << "opt::bamlist:" << opt::bamlist[i] << "\n";
//    }
//    std::cerr << "opt::bamlistfile:" << bamlistfile << "\n";
//    std::cerr << "opt::writerheader:" << opt::writerheader << "\n";
//    std::cerr << "opt::tel_k:" << opt::tel_k << "\n";
//    std::cerr << "opt::outputdir:" << opt::outputfile << "\n";

}



int main(int argc, char** argv)
{
	Timer* pTimer = new Timer("scan BAM");
	parseScanOptions(argc, argv);
    scanBam();
    delete pTimer;
    return 0;
}


