#ifndef LOWSIGNIFICANCECORRECTION_H
#define LOWSIGNIFICANCECORRECTION_H

#include <string>
#include <vector>

#include <Eigen/Dense>
#include <jlibrary.h>
#include <VectorclassNew4/vectorclass.h>
#include <Posclass.h>

using namespace std;
using namespace Eigen;


struct ParametersLS
{
    double threshold;
    double maxr;
    double r2impute;
    bool correct_r2;
    double base_num;
    bool compute_uncorrected = false;
    int ancestry = 0;
};

struct SumStatsLS{
    string rsnum;
    Posclass posclass;
    double eaf;
    double genetic_position;
    vector<double> w;   // weights from PRS weights file
    double b;   // summary statistics coefficient estimate
    double s;   // summary statistics standard error estimate
    double p;
    VectorXd ref_dosages; // reference dosages for SNP
};

struct PGMweights{
    string rsnum;
    Posclass posclass;
    double effect;
    double eaf;
};

vector<vector<PGMweights>> CorrectionForLSFromSummaryStatistics(const string& summarystats_file, const string& bgenix_reffile, const string& ref_file, const map<Posclass, vector<double>>& prsmap, const string& include_file, string geneticmap_file, double genetic_distance, int chr, const ParametersLS& parameters);

vector<SumStatsLS> GetSelectedLSstatistics(const map<Posclass, vector<double>>& prsmap, const map<Posclass, long>& position_map, const string& summarystats_file, const string& ref_file, const string& include_file, const string& geneticmap_file, int chr);

void GetLowSignificancePGS(const string& summarystats_file, const string& pgm_file, const string& hybrid_file, const string& bgen_ref_file_prefix, const string& bgenix_reffile_prefix, const string& include_file, const string& geneticmap_file, double genetic_dist, const string& outputprefix, const ParametersLS& parameters);


#endif // LOWSIGNIFICANCECORRECTION_H
