
#include <iostream>
#include <map>
#include "Util.h"

// Strip the leadering directories from a filename
std::string stripDirectories(const std::string& filename)
{
    size_t lastDirPos = filename.find_last_of('/');
    
    if(lastDirPos == std::string::npos)
        return filename; // no directories
    else
        return filename.substr(lastDirPos + 1);
}

// Returns the size of the file. Code from stackoverflow.
std::ifstream::pos_type getFilesize(const std::string& filename)
{
    std::ifstream in(filename.c_str(), std::ifstream::in | std::ifstream::binary);
    in.seekg(0, std::ifstream::end);
    return in.tellg();
}

// Open a file that may or may not be gzipped for reading
// The caller is responsible for freeing the handle
std::istream* createReader(const std::string& filename, std::ios_base::openmode mode)
{
	std::ifstream* pReader = new std::ifstream(filename.c_str(), mode);
	assertFileOpen(*pReader, filename);
	return pReader;
}

// Open a file that may or may not be gzipped for writing
// The caller is responsible for freeing the handle
std::ostream* createWriter(const std::string& filename,
                           std::ios_base::openmode mode)
{
	std::ofstream* pWriter = new std::ofstream(filename.c_str(), mode);
	assertFileOpen(*pWriter, filename);
	return pWriter;
}

// Ensure a filehandle is open
void assertFileOpen(std::ifstream& fh, const std::string& fn)
{
    if(!fh.is_open())
    {
        std::cerr << "Error: could not open " << fn << " for read\n";
        exit(EXIT_FAILURE);
    }
}

// Ensure a filehandle is open
void assertFileOpen(std::ofstream& fh, const std::string& fn)
{
    if(!fh.is_open())
    {
        std::cerr << "Error: could not open " << fn << " for write\n";
        exit(EXIT_FAILURE);
    }
}
