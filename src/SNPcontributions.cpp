#include "SNPcontributions.h"
#include "routines.h"

#include <iostream>
#include <string>
#include <Eigen/Dense>
#include <Eigen/Sparse>

#include "bgen/parser.hpp"
#include "bgen/IndexQuery.hpp"

#include <jlibrary.h>

using namespace std;
using namespace Eigen;
using  SpMat = SparseMatrix<double>;


struct PGSstats{
    string rsnum;
    Posclass posclass;
    double genetic_position;
    double w;   // weights from PRS weights file
    double s;   // standard error estimate for statistics
    VectorXd ref_dosages; // reference dosages for SNP
};

// Get SNP statisitics of snps that are in the reference file
// Column 2 chromosome
// Column 3 position
// Column 4 baseline
// Column 5 effect
// Column 8 odds ratio
// Column 9 standard error

vector<PGSstats> GetSelectedSNPstatistics(const map<Posclass, vector<double>>& prsmap, const map<Posclass, long>& position_map, const string& ref_file, const string& include_file, const string& geneticmap_file)
{
    vector<PGSstats> results;

    BgenParser bgenParser(ref_file);

    vector<int> include;

    if(include_file != ""){
        include = FileRead<vector<int>>(include_file);
    }

    string origsnpstr, snpstr;

    Imputeclass imputeclass;

    struct FlipCheck{
        string ref;
        string alt;
    };

    auto genmap = ReadGeneticMap(geneticmap_file);

    for(const auto& [posclass,coeff]: prsmap){
        try {
            if(auto itor2 = position_map.find(posclass);itor2!=position_map.end()){

                bgenParser.JumpToPosition(itor2->second);

                if(include_file != ""){
                    bgenParser.ReadImputeclass(imputeclass, include);
                }
                else{
                    bgenParser.ReadAllImputeclass(imputeclass);
                }

                if(isfinite(imputeclass.r2) && imputeclass.r2 > 0){
                    double mean = 2.0*imputeclass.eaf;

                    try{
                        results.push_back(PGSstats());

                        results.back().w = coeff[0];

                        results.back().posclass = posclass;

                        results.back().s = sqrt(2.0*imputeclass.eaf*(1-imputeclass.eaf));

                        results.back().genetic_position = CalcGeneticPosition(genmap, posclass.position);

                        results.back().ref_dosages = VectorXd::Zero(imputeclass.size());

                        forc(i, imputeclass){
                            if(imputeclass[i]>-0.5){
                                results.back().ref_dosages[i] = imputeclass[i] - mean;
                            }
                            else{
                                results.back().ref_dosages[i]=0;
                            }
                        }

                        results.back().ref_dosages = results.back().ref_dosages/results.back().ref_dosages.norm();

                        int flip2 = MatchAlleles(imputeclass, FlipCheck{posclass.ref, posclass.alt});

                        if(flip2 == -1){
                            results.back().w *= -1;
                        }
                    }
                    catch(...){
                    }
                }
            }
        } catch (...) {

        }
    }

    return results;
}


void CalcSNPcontributions(string pgm_file, string bgenix_reffile_prefix, string bgen_ref_file_prefix, string include_file, string geneticmap_file, double genetic_dist)
{
    auto prsmap = GeneratePRSmap(pgm_file);

    set<int> chrset;

    for(auto& [v,w]: prsmap){
        chrset.insert(v.chr);
    }

    double quadprod = 0;

    for(auto& i:chrset){
        string bgenix_reffile = Paste(bgenix_reffile_prefix,i,".txt");

        string bgen_ref_file = Paste(bgen_ref_file_prefix,i,".bgen");

        vector<string> rsids=[&]{
            if(bgenix_reffile_prefix!=""){
                auto unique_ids=ConvertPositionsToIds(Paste(bgenix_reffile_prefix,i,".txt"), prsmap, false);

                return vector<string>(unique_ids.begin(), unique_ids.end());
            }
            else{
                return vector<string>();
            }
        }();

        auto position_map = GetPositionsFromIds(rsids, Paste(bgen_ref_file_prefix,i,".bgen"));

        auto summary_stats = GetSelectedSNPstatistics(prsmap, position_map, bgen_ref_file, include_file, geneticmap_file);

        int ndim = summary_stats.size();

        using Tlet = Eigen::Triplet<double>;

        vector<Tlet> triplet;

        forl(i, ndim){
            triplet.push_back(Tlet(i,i,1));
            forl(j, i){
                if(abs(summary_stats[i].genetic_position - summary_stats[j].genetic_position) < genetic_dist){
                    double correlation = summary_stats[i].ref_dosages.dot(summary_stats[j].ref_dosages);

    //                cout << j << '\t' << i << '\t' << correlation << '\n';

                    triplet.push_back(Tlet(i,j,correlation));

                    triplet.push_back(Tlet(j,i,correlation));
                }
            }
        }

        SpMat corrmat(ndim, ndim);

        corrmat.setFromTriplets(triplet.begin(), triplet.end());

    // get b and s values

        VectorXd s(ndim), w(ndim);

        forl(i, ndim){
            s[i] = summary_stats[i].s;
            w[i] = summary_stats[i].w;

            if(!isfinite(s[i])){
                cout << summary_stats[i].rsnum << ' ' << summary_stats[i].posclass << '\n';
            }
        }


        forl(i, ndim){
            if(s[i]<1e-5 || !isfinite(s[i]) || !isfinite(w[i])){
                cout << i << ' ' << s[i] << ' ' << w[i] << '\n';
            }
        }

        VectorXd x = w.cwiseProduct(s);

        VectorXd product = corrmat * x;

        VectorXd diff = 2.0*product.cwiseProduct(x) - x.cwiseAbs2();

        quadprod += (x.transpose() * corrmat * x).value();

        forc(i , x){
            cout << summary_stats[i].posclass << '\t' << diff[i] << '\t' << Sqr(x[i]) << '\n';
        }
    }
}
