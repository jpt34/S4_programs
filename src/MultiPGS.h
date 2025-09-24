#ifndef MULTIPGS_H
#define MULTIPGS_H

#include <string>
#include <vector>

using namespace std;

void MultiPGSregression(const string& summarystats_file, const vector<string>& pgm_file, const string& bgenix_reffile_prefix, const string& bgen_ref_file_prefix, const string& include_file, const string& chrX_include_file, const string& geneticmap_file, double genetic_dist);

void PrintCorrelationMatrix(const string& bgenix_reffile, const string& ref_file, const vector<string>& pgm_file, const string& include_file, string geneticmap_file, double genetic_distance, const string& outputfile);

void MultiPGSregressionFromCorrelationFiles(const string& summarystats_file, const vector<string>& pgm_file, const string& corr_file_prefix);

#endif // MULTIPGS_H
