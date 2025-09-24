#include "LowSignificanceCorrection.h"
#include "routines.h"

#include <iostream>
#include <string>
#include <regex>
#include <Eigen/Dense>

#include "bgen/parser.hpp"
#include "bgen/IndexQuery.hpp"

#include <jlibrary.h>

using namespace std;
using namespace Eigen;

// Use summary statistics to predict the accuracy of PRS weights
// prs_file is the same as the file used to generate PRS weights

vector<vector<PGMweights>> CorrectionForLSFromSummaryStatistics(const string& summarystats_file, const string& bgenix_reffile, const string& ref_file, const map<Posclass, vector<double>>& prsmap, const string& include_file, string geneticmap_file, double genetic_distance, int chr, const ParametersLS& parameters)
{
    vector<string> rsids=[&]{
        if(bgenix_reffile!=""){ // if the bgenix reference file is not empty
            auto unique_ids=ConvertPositionsToIds(bgenix_reffile, prsmap, false); // convert the positions to ids using the bgenix reference file

            return vector<string>(unique_ids.begin(), unique_ids.end()); // return the ids as a vector of strings
        }
        else{
            return vector<string>(); // otherwise return an empty vector
        }
    }();

    auto position_map = GetPositionsFromIds(rsids, ref_file); // get a map of positions from ids using the reference file

    auto summary_stats = GetSelectedLSstatistics(prsmap, position_map, summarystats_file, ref_file, include_file, geneticmap_file, chr); // get the summary statistics for the selected SNPs

    vector<SumStatsLS> correction_snps;

    vector<SumStatsLS> hybrid_snps;

    for(auto& v: summary_stats){
        if(v.w.size()>0){
            correction_snps.push_back(v);
        }
        else{
            hybrid_snps.push_back(v);
        }
    }

    vector<PGMweights> correction;

    vector<PGMweights> correction_se;

    vector<PGMweights> hybrid;

    vector<PGMweights> hybrid_se;

//    vector<PGMweights> orig;

//    vector<PGMweights> orig_se;

    PGMweights tempcorrection, tempcorrection_se;

    for(auto& v: hybrid_snps){
        double r2_v = 1.0/(v.eaf*(1-v.eaf)*Sqr(v.s)*parameters.base_num);

        if(r2_v>1.0) r2_v = 1.0;

        double correction_effect = 0;

        bool pass_corr = r2_v > parameters.r2impute;

        for(auto &w: correction_snps){
            if(abs(v.genetic_position - w.genetic_position) < genetic_distance){
                double correlation = v.ref_dosages.dot(w.ref_dosages);

                double r2_w = 1.0/(w.eaf*(1-w.eaf)*Sqr(w.s)*parameters.base_num);

                if(r2_w>1.0) r2_w = 1.0;


                if(parameters.correct_r2){
                    double corr_correct = sqrt(min(r2_v, r2_w)/max(r2_v, r2_w));

                    correlation *= corr_correct;
                }

                correction_effect -= correlation*w.w[0]/w.s;

                if(correlation > parameters.maxr){
                    pass_corr = false;
                }
            }
        }

        tempcorrection.eaf = tempcorrection_se.eaf = v.eaf;

        tempcorrection.posclass = tempcorrection_se.posclass = v.posclass;

        tempcorrection.rsnum = tempcorrection_se.rsnum = v.rsnum;

        if(parameters.compute_uncorrected){
            tempcorrection.effect = v.b;

            tempcorrection_se.effect = v.b / (Sqr(v.s) * parameters.base_num);
        }
        else{
            tempcorrection.effect = correction_effect * v.s;

            tempcorrection_se.effect = correction_effect / (v.s * parameters.base_num);
        }

        if(pass_corr){
            correction.push_back(tempcorrection);
            correction_se.push_back(tempcorrection_se);
        }
    }

    return {correction, correction_se, hybrid, hybrid_se};
}


// Get SNP statisitics of snps that are in the reference file
// Column 2 chromosome
// Column 3 position
// Column 4 baseline
// Column 5 effect
// Column 8 odds ratio
// Column 9 standard error

vector<SumStatsLS> GetSelectedLSstatistics(const map<Posclass, vector<double>>& prsmap, const map<Posclass, long>& position_map, const string& summarystats_file, const string& ref_file, const string& include_file, const string& geneticmap_file, int chr)
{
    vector<SumStatsLS> results;

    BgenParser bgenParser(ref_file);

    vector<int> include;

    if(include_file != ""){
        include = FileRead<vector<int>>(include_file);
    }

    string origsnpstr, snpstr;
    SplitStringView splitstr;

    Imputeclass imputeclass;

    struct FlipCheck{
        string ref;
        string alt;
    };

    auto genmap = ReadGeneticMap(geneticmap_file);

    zifstream input(summarystats_file);

    while(input >> splitstr){
        try {
            Posclass posclass(chr, stoi(string(splitstr[1])), splitstr[3], splitstr[2]);

            if(auto itor_prsmap = GetPosItor(prsmap, posclass);itor_prsmap!=prsmap.end()){     // snp found in map
                if(auto itor2 = GetPosItor(position_map, posclass);itor2!=position_map.end()){

                    bgenParser.JumpToPosition(itor2->second);

                    if(include_file != ""){
                        bgenParser.ReadImputeclass(imputeclass, include);
                    }
                    else{
                        bgenParser.ReadAllImputeclass(imputeclass);
                    }

                    if(isfinite(imputeclass.r2) && imputeclass.r2 > 0){
                        int flip = MatchAlleles(imputeclass, FlipCheck{posclass.ref, posclass.alt});    // are the alleles in the same order as the test statistics

                        double mean = 2.0*imputeclass.eaf;

                        try{
                            double bvalue = stod(string(splitstr[6]));

                            double svalue = stod(string(splitstr[7]));

                            if(isfinite(bvalue) && isfinite(svalue) && svalue>0){
                                results.push_back(SumStatsLS());

                                results.back().w = itor_prsmap->second;

                                results.back().posclass = posclass;

                                results.back().b = bvalue;

                                results.back().s = svalue;

                                results.back().rsnum = string(splitstr[0]);

                                results.back().p = stod(string(splitstr[8]));

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

                                results.back().ref_dosages = flip*results.back().ref_dosages/results.back().ref_dosages.norm();

                                int flip2 = MatchAlleles(FlipCheck{posclass.ref, posclass.alt}, FlipCheck{itor_prsmap->first.ref, itor_prsmap->first.alt});

                                if(flip2 == -1){
                                    for(auto& value:results.back().w)
                                    value *= -1;    // won't matter if frequency as not used here
                                }

                                results.back().eaf = stod(string(splitstr[4]));
                            }
                        }
                        catch(...){
                        }
                    }
                }
            }
        } catch (...) {
        }
    }

    return results;
}


void GetLowSignificancePGS(const string& summarystats_file, const string& pgm_file, const string& hybrid_file, const string& bgen_ref_file_prefix, const string& bgenix_reffile_prefix, const string& include_file, const string& geneticmap_file, double genetic_dist, const string& outputprefix, const ParametersLS& parameters)
{
//    map<Posclass, vector<double>> prsmap;  // create a map of position classes and vectors of doubles

    int nrows = pgm_file.size(); // get the number of rows from the pgm file

    auto prsmap = GeneratePRSmap(pgm_file);

//    forc(i, pgm_file){ // loop over the pgm file
//        auto tempmap = GeneratePRSmap(pgm_file[i]); // generate a PRS map from each file

//        for(auto& [w1,w2]:tempmap){ // loop over the PRS map
//            if(prsmap[w1].empty()){ // if the prsmap for the position class is empty
//                prsmap[w1].resize(nrows,0.0); // resize it to nrows and fill with zeros
//            }

//            prsmap[w1][i] = w2[0]; // assign the first value of w2 to the prsmap
//        }
//    }

    {
        zifstream input(hybrid_file);

        string rsnum, ref, alt;

        int chr, position;

        while(input >> rsnum >> chr >> position >> ref >> alt){

            Posclass temppos(chr, position, ref, alt);

            if(prsmap.count(temppos)==0){
                prsmap[temppos] ={};
            }
        }
    }

    set<int> chrset; // create a set to store the chromosomes

    for(auto& [v,w]: prsmap){ // loop over the prsmap
        chrset.insert(v.chr); // insert the chromosome of each position class into the set
    }

    regex pattern(R"(chr\d+)"); // create a regex pattern for chromosome numbers

    VectorXd numerator_sum = VectorXd::Zero(nrows); // create a vector for numerator sum and initialize with zeros

    MatrixXd corr_sum = MatrixXd::Zero(nrows, nrows); // create a matrix for corr sum and initialize with zeros


    zofstream output(Paste(outputprefix, ".txt"));

    zofstream output_se(Paste(outputprefix, "_se.txt"));

//    zofstream output_hybrid(Paste(outputprefix, "_hybrid.txt"));

//    zofstream output_hybrid_se(Paste(outputprefix, "_hybrid_se.txt"));

    for(auto& i:chrset){ // loop over the chromosomes
        string replacement = Paste("chr",i); // create a replacement string with chr and i

        string bgenix_reffile = regex_replace(bgenix_reffile_prefix, pattern, replacement); // replace the bgenix reference file prefix with the replacement

        string bgen_ref_file = regex_replace(bgen_ref_file_prefix, pattern, replacement); // replace the bgen reference file prefix with the replacement

        string geneticmap_file_chr = regex_replace(geneticmap_file, pattern, replacement); // replace the genetic map file with the replacement

        string summarystats_chr = regex_replace(summarystats_file, pattern, replacement);

        auto corrections = CorrectionForLSFromSummaryStatistics(summarystats_chr, bgenix_reffile, bgen_ref_file, prsmap, include_file, geneticmap_file_chr, genetic_dist, i, parameters); // estimate the numerator and corr from summary statistics

        for(auto& v:corrections[0]){
            Printtabline(output, v.rsnum, v.posclass, v.effect, v.eaf);
        }

        for(auto& v:corrections[1]){
            Printtabline(output_se, v.rsnum, v.posclass, v.effect, v.eaf);
        }
    }
}
