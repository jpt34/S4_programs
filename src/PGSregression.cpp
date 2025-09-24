#include "PGSregression.h"
#include "routines.h"

#include <iostream>
#include <string>
#include <regex>
#include <Eigen/Dense>


#include "bgen/parser.hpp"
#include "bgen/IndexQuery.hpp"

#include <jlibrary.h>


// Get statistics for PRS values

vector<vector<double>> GetPRSstats(const string& testphenotypefile, const vector<string>& testfile, const string& bgenix_prefix, const vector<string>& chrx_files, const string& analysis_group, map<Posclass, vector<double>>& prsmap, const string& weightings_file, const vector<vector<double>>& prs_values){

    auto phenotypedata = GetPhenotypeData(testphenotypefile, analysis_group);

    // generate vector of prs values
    // if prsmap is empty than generate values from the bgen file, otherwise use the prsmap values
    // uses a lambda to initialise the vector

    const vector<VectorXd> prs = [&]{
        if(prs_values.empty()){
            return vector<VectorXd>{GetPrincipalComponentValuesFromSeveralBgenFiles(testfile, prsmap, {0,"", false}, bgenix_prefix, chrx_files)};
        }
        else{
            vector<VectorXd> result;

            for(const auto& v:prs_values){
                result.push_back((VectorXd(Map<const VectorXd>(v.data(), v.size()))));
            }

            return result;
        }
    }();

    // if weightings file is specified then write PRS values to that file

    if(weightings_file!=""){
        zofstream output(weightings_file);

        output << "prs\n";

        forc(i, prs[0]){
            output << prs[0][i] << '\n';
        }
    }

    auto prsstats = (phenotypedata.status.size()!=0?GetPRSchi2(phenotypedata, prs):GetRegressionPRS(phenotypedata, prs));

    forc(j, prs){
        VectorXd prs_filter(phenotypedata.study.size());

        int filternum = 0;

        forc(i, phenotypedata.include){
            if(phenotypedata.include[i]==1){
                prs_filter[filternum++] = prs[j][i];
            }
        }

        double auc = (phenotypedata.status.size()!=0?GetStratRankSumStat(phenotypedata.status, phenotypedata.study, prs_filter):GetPearsonCorrelation(phenotypedata.regval, prs_filter));

        prsstats[j].push_back(auc);

        double prs_se = sqrt((prs_filter.array() - prs_filter.mean()).square().sum()/(prs_filter.size()-1));

        prsstats[j].push_back(prs_se);
    }

    return prsstats;
}

// print statistics for PRS, either from analysis weights and bgen file, or from prs doses

void PrintPRSstats(const string& testphenotypefile, const vector<string>& testfile, const string& bgenix_prefix, const vector<string>& chrx_files, const string& analysis_group, const vector<string>& prs_file, int cv_num)
{
    map<Posclass, vector<double>> prsmap;

    string first_file = (!prs_file.empty()?prs_file[0]:"");

    zifstream prs_input(first_file);

    string rsnum;

    Posclass posclass;

    double oddsratio, eaf;

    // ismapfile determines whether the file is just dosages when the number of columns will be one

    bool ismapfile = [&]{
        zifstream input(first_file);
        SplitString splitstr;
        input >> splitstr;
        return splitstr.size() != 1;
    }();

    vector<vector<double>> prs;

    if(ismapfile){
        while(prs_input >> rsnum >> posclass >> oddsratio >> eaf){

            prsmap[posclass] = {oddsratio, eaf};
        }
    }
    else{

        for(const auto& v:prs_file){

            zifstream input(v);

            input >> skipline;

            double temp;
// Push an empty vector to prs
            prs.push_back({});
// Read the values from the file until the end of the file and push each into the last vector in prs
            while(input >> temp){
                prs.back().push_back(temp);
            }
        }
    }

    auto prsstats = GetPRSstats(testphenotypefile, testfile, bgenix_prefix, chrx_files, analysis_group, prsmap, "", prs);

    for(auto& v:prsstats){

        if(v.size()>5){
            Printtabline(cout, v[0], v[1]*v[5], v[2], v[4], v[5], v[0]*v[5], (v[0]-1.96*v[1])*v[5], (v[0]+1.96*v[1])*v[5], v[3]);
        }
        else{
            Printtabline(cout, v[0], v[1]*v[4], v[2], v[3], v[4], v[0]*v[4], (v[0]-1.96*v[1])*v[4], (v[0]+1.96*v[1])*v[4]);
        }
    }
}
