/*
 * merge.h
 *
 *  Created on: 3 Sep 2013
 *      Author: zd1
 */

#ifndef MERGE_H_
#define MERGE_H_

#include <getopt.h>
#include "config.h"
#include "scan.h"

int mergeMain(int argc, char** argv);
int mergeBam();
ScanResults parseResultFile(std::string filename);
void parseMergeOptions(int argc, char** argv);
int printsummary(std::vector<ScanResults> &results, std::string filepath);

#endif /* MERGE_H_ */
