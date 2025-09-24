#ifndef PGSREGRESSION_H
#define PGSREGRESSION_H

#include <string>
#include <vector>
#include <map>

#include <Eigen/Dense>

#include <jlibrary.h>
#include <Posclass.h>

using namespace std;
using namespace Eigen;

vector<vector<double> > GetPRSstats(const string& testphenotypefile, const vector<string>& testfile, const string& bgenix_prefix, const vector<string>& chrx_files, const string& analysis_group, map<Posclass, vector<double>>& prsmap, const string& weightings_file, const vector<vector<double>>& prs_values={});

void PrintPRSstats(const string& testphenotypefile, const vector<string>& testfile, const string& bgenix_prefix, const vector<string>& chrx_files, const string& analysis_group, const vector<string> &prs_file, int cv_num);

#endif // PGSREGRESSION_H
