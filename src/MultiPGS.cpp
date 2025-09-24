#include "MultiPGS.h"
#include "routines.h"

#include <iostream>
#include <string>
#include <regex>
#include <Eigen/Dense>
#include <Eigen/Sparse>

#include "bgen/parser.hpp"
#include "bgen/IndexQuery.hpp"

#include <jlibrary.h>

using namespace std;
using namespace Eigen;
using  SpMat = SparseMatrix<double>;


// Use summary statistics to predict the accuracy of PRS weights
// prs_file is the same as the file used to generate PRS weights


tuple<VectorXd, MatrixXd> EstimateMultiPGSFromSummaryStatistics(const string& summarystats_file, const string& bgenix_reffile, const string& ref_file, const map<Posclass, vector<double>>& prsmap, const string& include_file, string geneticmap_file, double genetic_distance)
{    
    int nrows = prsmap.begin()->second.size(); // get the number of rows from the prsmap

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

    auto summary_stats = GetSelectedSNPstatistics(prsmap, position_map, summarystats_file, ref_file, include_file, geneticmap_file); // get the summary statistics for the selected SNPs

    int ndim = summary_stats.size(); // get the number of dimensions from the summary statistics

    using Tlet = Eigen::Triplet<double>; // define a type alias for Eigen triplets

    vector<Tlet> triplet; // create a vector of triplets

    forl(i, ndim){ // loop over the dimensions
        triplet.push_back(Tlet(i,i,1)); // push back a triplet with diagonal value 1
        forl(j, i){ // loop over the lower triangle
            if(abs(summary_stats[i].genetic_position - summary_stats[j].genetic_position) < genetic_distance){ // if the genetic distance between two SNPs is less than the threshold
                double correlation = summary_stats[i].ref_dosages.dot(summary_stats[j].ref_dosages); // calculate the correlation between their reference dosages

                triplet.push_back(Tlet(i,j,correlation)); // push back a triplet with lower triangular value correlation

                triplet.push_back(Tlet(j,i,correlation)); // push back a triplet with upper triangular value correlation
            }
        }
    }

    SpMat corrmat(ndim, ndim); // create a sparse matrix for the correlation matrix

    corrmat.setFromTriplets(triplet.begin(), triplet.end()); // set the sparse matrix from the triplets


// get b and s values

    VectorXd b(ndim); // create a vector for b values

    MatrixXd w(nrows, ndim); // create a matrix for w values

    forc(i, b){ // loop over b values
        b[i] = summary_stats[i].b/summary_stats[i].s; // calculate b[i] as b/s

        w.col(i) = Map<VectorXd>(summary_stats[i].w.data(), nrows)/summary_stats[i].s; // map w column i as w/s

        if(!isfinite(summary_stats[i].s)){ // if s is not finite
            cout << summary_stats[i].rsnum << ' ' << summary_stats[i].posclass << '\n'; // print rsnum and posclass
        }
    }

    MatrixXd corr = w * corrmat * w.transpose(); // calculate corr as w * corrmat * w'

    VectorXd numerator = w * b; // calculate numerator as w * b

    return {numerator, corr}; // return a tuple of numerator and corr
}


tuple<VectorXd, MatrixXd> EstimateMultiPGSFromSummaryStatisticsAndCorrelationFile(const string& summarystats_file, const string& corr_file, const map<Posclass, vector<double>>& prsmap)
{
    int nrows = prsmap.begin()->second.size(); // get the number of rows from the prsmap

    zifstream corr_input(corr_file);

    map<Posclass, int> corr_indices;

    using Tlet = Eigen::Triplet<double>; // define a type alias for Eigen triplets

    SplitStringView splitstr;

    int ind1, ind2;

    double corrval;

    while(corr_input >> splitstr){
        if(splitstr.size()==5){
            corr_indices[Posclass{ConvertChrStr(string(splitstr[0])), stoi(string(splitstr[1])), splitstr[2], splitstr[3]}] = stoi(string(splitstr[4]));
        }
        else if(splitstr.size()==3){
            ind1 = stoi(string(splitstr[0]));
            ind2 = stoi(string(splitstr[1]));
            corrval = stod(string(splitstr[2]));
            break;
        }

//        if(corr_triplets.size()%1000==0){
//            cout << "Read " << corr_triplets.size() << " correlations\n";
//            cout << "line is " << splitstr.str() << "\n";
//        }
    }

    vector<SumStatsPRS> summary_stats;


    string origsnpstr, snpstr;

    struct FlipCheck{
        string ref;
        string alt;
    };


    zifstream input(summarystats_file);

    map<int, int> index_conversions;

    while(input >> splitstr){
        try {
            Posclass posclass(ConvertChrStr(string(splitstr[1])), stoi(string(splitstr[2])), splitstr[4], splitstr[3]);

            if(auto itor_prsmap = GetPosItor(prsmap, posclass);itor_prsmap!=prsmap.end()){     // snp found in map
                if(auto itor_corrmatrix = GetPosItor(corr_indices, posclass);itor_corrmatrix!=corr_indices.end()){

                    FlipCheck ref_flip{itor_corrmatrix->first.ref, itor_corrmatrix->first.alt};

                    int flip_sumstats = MatchAlleles(ref_flip, FlipCheck{posclass.ref, posclass.alt});    // are the alleles in the same order as the test statistics

                    try{
                        double bvalue = flip_sumstats * stod(string(splitstr[7]));

                        double svalue = stod(string(splitstr[8]));

                        if(isfinite(bvalue) && isfinite(svalue) && svalue>0){
                            summary_stats.push_back(SumStatsPRS());

                            summary_stats.back().w = itor_prsmap->second;

                            summary_stats.back().posclass = posclass;

                            summary_stats.back().b = bvalue;

                            summary_stats.back().s = svalue;

                            index_conversions[itor_corrmatrix->second] = int(summary_stats.size())-1;

                            int flip_prs = MatchAlleles(ref_flip, FlipCheck{itor_prsmap->first.ref, itor_prsmap->first.alt});

                            if(flip_prs == -1){
                                for(auto& value:summary_stats.back().w)
                                    value *= -1;    // won't matter if frequency as not used here
                            }
                        }
                    }
                    catch(...){
                    }
                }
            }
        } catch (...) {
        }
    }

//    cout << "Summary statistics read\n";

    int ndim = summary_stats.size(); // get the number of dimensions from the summary statistics

    SpMat corrmat(ndim, ndim); // create a sparse matrix for the correlation matrix

    // This might read values faster now that there are 3 columns

    {

        vector<Tlet> converted_triplets;

        do{
            if(index_conversions.count(ind2) && index_conversions.count(ind1)){
                converted_triplets.push_back(Tlet(index_conversions[ind2], index_conversions[ind1], corrval));

                if(ind2!=ind1){
                    converted_triplets.push_back(Tlet(index_conversions[ind1], index_conversions[ind2], corrval));
                }
            }
        }while(corr_input >> ind1 >> ind2 >> corrval);


        // Check for duplicates and then delete set

        {
            set<pair<int, int>> check_duplicates;

            for(auto& v: converted_triplets){
                if(check_duplicates.insert({v.col(), v.row()}).second == false){
                    cerr << "Duplicate for elements {" << v.col() << "," << v.row() << "}\n";
                }
            }
        }

        if(converted_triplets.size()==0){
            return {VectorXd::Zero(nrows), MatrixXd::Zero(nrows, nrows)};
        }

        corrmat.setFromTriplets(converted_triplets.begin(), converted_triplets.end()); // set the sparse matrix from the triplets
    }

//    cout << "corrmat created\n";

// get b and s values

    VectorXd b(ndim); // create a vector for b values

    MatrixXd w(nrows, ndim); // create a matrix for w values

    forc(i, b){ // loop over b values
        b[i] = summary_stats[i].b/summary_stats[i].s; // calculate b[i] as b/s

        w.col(i) = Map<VectorXd>(summary_stats[i].w.data(), nrows)/summary_stats[i].s; // map w column i as w/s

        if(!isfinite(summary_stats[i].s)){ // if s is not finite
            cout << summary_stats[i].rsnum << ' ' << summary_stats[i].posclass << '\n'; // print rsnum and posclass
        }
    }

//    cout << "Before matrix calculation\n";

    MatrixXd corr = w * corrmat * w.transpose(); // calculate corr as w * corrmat * w'

    VectorXd numerator = w * b; // calculate numerator as w * b

//    cout << "Matrix calculted\n";

    return {numerator, corr}; // return a tuple of numerator and corr
}


void MultiPGSregression(const string& summarystats_file, const vector<string>& pgm_file, const string& bgenix_reffile_prefix, const string& bgen_ref_file_prefix, const string& include_file, const string& chrX_include_file, const string& geneticmap_file, double genetic_dist)
{
    map<Posclass, vector<double>> prsmap;  // create a map of position classes and vectors of doubles

    int nrows = pgm_file.size(); // get the number of rows from the pgm file

    forc(i, pgm_file){ // loop over the pgm file
        auto tempmap = GeneratePRSmap(pgm_file[i]); // generate a PRS map from each file

        for(auto& [w1,w2]:tempmap){ // loop over the PRS map
            if(prsmap[w1].empty()){ // if the prsmap for the position class is empty
                prsmap[w1].resize(nrows,0.0); // resize it to nrows and fill with zeros
            }

            prsmap[w1][i] = w2[0]; // assign the first value of w2 to the prsmap
        }
    }

    set<int> chrset; // create a set to store the chromosomes

    for(auto& [v,w]: prsmap){ // loop over the prsmap
        chrset.insert(v.chr); // insert the chromosome of each position class into the set
    }

    regex pattern(R"(chr\d+)"); // create a regex pattern for chromosome numbers

    VectorXd numerator_sum = VectorXd::Zero(nrows); // create a vector for numerator sum and initialize with zeros

    MatrixXd corr_sum = MatrixXd::Zero(nrows,  nrows); // create a matrix for corr sum and initialize with zeros

    for(auto& i:chrset){ // loop over the chromosomes
        string replacement = Paste("chr",i); // create a replacement string with chr and i

        string bgenix_reffile = regex_replace(bgenix_reffile_prefix, pattern, replacement); // replace the bgenix reference file prefix with the replacement

        string bgen_ref_file = regex_replace(bgen_ref_file_prefix, pattern, replacement); // replace the bgen reference file prefix with the replacement

        string geneticmap_file_chr = regex_replace(geneticmap_file, pattern, replacement); // replace the genetic map file with the replacement

        string used_include_file = ((i==23 && chrX_include_file!="")?chrX_include_file:include_file);

        auto [numerator, corr] = EstimateMultiPGSFromSummaryStatistics(summarystats_file, bgenix_reffile, bgen_ref_file, prsmap, used_include_file, geneticmap_file_chr, genetic_dist); // estimate the numerator and corr from summary statistics

        numerator_sum += numerator; // add the numerator to the numerator sum

        corr_sum += corr; // add the corr to the corr sum

        if(nrows>1){
            cout << "chr" << i << "\n"; // print the chromosome number

            cout << numerator << "\n"; // print the numerator vector

            cout << corr << endl; // print the corr matrix and flush the output stream
        }
        else{
            cout << "chr" << i << '\t' << numerator[0] << '\t' << corr(0,0) << endl; // print the values on the same line if only one PGS
        }
    }

    MatrixXd variance = corr_sum.inverse(); // calculate the variance matrix as the inverse of the corr sum matrix

    VectorXd coeff = variance*numerator_sum;

    cout << "All results\n";

    if(nrows>1){
        cout << numerator_sum << "\n\n"; // print the numerator sum

        cout << corr_sum << "\n"; // print the corr sum

        cout << "Estimated beta\n";

        cout << coeff << "\n";   // print the product of the variance matrix and the numerator sum vector

        cout << "Estimated variance matrix\n";

        cout << variance << '\n';   // print the variance matrix
    }
    else{
        cout << numerator_sum[0] << '\t' << corr_sum(0,0) << '\n';
    }

    cout << "Coefficients, standard errors and chi2\n";

    for(auto & v:coeff){
        cout << v << '\t';      // print each coeff value
    }

    forl(i, nrows){
        cout << sqrt(variance(i,i)) << '\t';   // print the square root of the diagonal elements
    }

    forl(i, nrows){
        cout << Sqr(coeff[i])/variance(i,i) << (i==nrows-1?'\n':'\t');      // Get the chi2 values for the regression
    }
}


void MultiPGSregressionFromCorrelationFiles(const string& summarystats_file, const vector<string>& pgm_file, const string& corr_file_prefix)
{
    map<Posclass, vector<double>> prsmap;  // create a map of position classes and vectors of doubles

    int nrows = pgm_file.size(); // get the number of rows from the pgm file

    forc(i, pgm_file){ // loop over the pgm file
        auto tempmap = GeneratePRSmap(pgm_file[i]); // generate a PRS map from each file

        for(auto& [w1,w2]:tempmap){ // loop over the PRS map
            if(prsmap[w1].empty()){ // if the prsmap for the position class is empty
                prsmap[w1].resize(nrows,0.0); // resize it to nrows and fill with zeros
            }

            prsmap[w1][i] = w2[0]; // assign the first value of w2 to the prsmap
        }
    }

    set<int> chrset; // create a set to store the chromosomes

    for(auto& [v,w]: prsmap){ // loop over the prsmap
        chrset.insert(v.chr); // insert the chromosome of each position class into the set
    }

    regex pattern(R"(chr\d+)"); // create a regex pattern for chromosome numbers

    VectorXd numerator_sum = VectorXd::Zero(nrows); // create a vector for numerator sum and initialize with zeros

    MatrixXd corr_sum = MatrixXd::Zero(nrows,  nrows); // create a matrix for corr sum and initialize with zeros

//    cout << "PRS files read\n";

    for(auto& i:chrset){ // loop over the chromosomes
        string replacement = Paste("chr",i); // create a replacement string with chr and i

        string corr_file = regex_replace(corr_file_prefix, pattern, replacement); // replace the bgenix reference file prefix with the replacement

        auto [numerator, corr] = EstimateMultiPGSFromSummaryStatisticsAndCorrelationFile(summarystats_file, corr_file, prsmap); // estimate the numerator and corr from summary statistics

        numerator_sum += numerator; // add the numerator to the numerator sum

        corr_sum += corr; // add the corr to the corr sum

        if(nrows>1){
            cout << "chr" << i << "\n"; // print the chromosome number

            cout << numerator << "\n"; // print the numerator vector

            cout << corr << endl; // print the corr matrix and flush the output stream
        }
        else{
            cout << "chr" << i << '\t' << numerator[0] << '\t' << corr(0,0) << endl; // print the values on the same line if only one PGS
        }
    }

    MatrixXd variance = corr_sum.inverse(); // calculate the variance matrix as the inverse of the corr sum matrix

    VectorXd coeff = variance*numerator_sum;

    cout << "All results\n";

    if(nrows>1){
        cout << numerator_sum << "\n\n"; // print the numerator sum

        cout << corr_sum << "\n"; // print the corr sum

        cout << "Estimated beta\n";

        cout << coeff << "\n";   // print the product of the variance matrix and the numerator sum vector

        cout << "Estimated variance matrix\n";

        cout << variance << '\n';   // print the variance matrix
    }
    else{
        cout << numerator_sum[0] << '\t' << corr_sum(0,0) << '\n';
    }

    cout << "Coefficients, standard errors and chi2\n";

    for(auto & v:coeff){
        cout << v << '\t';      // print each coeff value
    }

    forl(i, nrows){
        cout << sqrt(variance(i,i)) << '\t';   // print the square root of the diagonal elements
    }

    forl(i, nrows){
        cout << Sqr(coeff[i])/variance(i,i) << (i==nrows-1?'\n':'\t');      // Get the chi2 values for the regression
    }
}

// Get SNP statisitics of snps that are in the reference file
// Column 2 chromosome
// Column 3 position
// Column 4 baseline
// Column 5 effect
// Column 8 odds ratio
// Column 9 standard error

vector<SumStatsPRS> GetSelectedDosages(const map<Posclass, vector<double>>& prsmap, const map<Posclass, long>& position_map, const string& ref_file, const string& include_file, const string& geneticmap_file)
{
    vector<SumStatsPRS> results;

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

    for(const auto& [posclass,w]: prsmap){
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
                    int flip = MatchAlleles(imputeclass, FlipCheck{posclass.ref, posclass.alt});    // are the alleles in the same order as the test statistics

                    double mean = 2.0*imputeclass.eaf;

                    results.push_back(SumStatsPRS());

                    results.back().posclass = posclass;

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
                }
            }
        }catch (...) {
        }
    }

    return results;
}


void PrintCorrelationMatrix(const string& bgenix_reffile, const string& ref_file, const vector<string>& pgm_file, const string& include_file, string geneticmap_file, double genetic_distance, const string& outputfile)
{
    map<Posclass, vector<double>> prsmap;  // create a map of position classes and vectors of doubles

    int nrows = pgm_file.size(); // get the number of rows from the pgm file

    forc(i, pgm_file){ // loop over the pgm file
        auto tempmap = GeneratePRSmap(pgm_file[i]); // generate a PRS map from each file

        for(auto& [w1,w2]:tempmap){ // loop over the PRS map
            if(prsmap[w1].empty()){ // if the prsmap for the position class is empty
                prsmap[w1].resize(nrows,0.0); // resize it to nrows and fill with zeros
            }

            prsmap[w1][i] = w2[0]; // assign the first value of w2 to the prsmap
        }
    }

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

    auto snp_dosages = GetSelectedDosages(prsmap, position_map, ref_file, include_file, geneticmap_file); // get the summary statistics for the selected SNPs

    int ndim = snp_dosages.size(); // get the number of dimensions from the summary statistics

    using Tlet = Eigen::Triplet<double>; // define a type alias for Eigen triplets

    vector<Tlet> triplet; // create a vector of triplets

    forl(i, ndim){ // loop over the dimensions
        triplet.push_back(Tlet(i,i,1)); // push back a triplet with diagonal value 1
        forl(j, i){ // loop over the lower triangle
            if(abs(snp_dosages[i].genetic_position - snp_dosages[j].genetic_position) < genetic_distance){ // if the genetic distance between two SNPs is less than the threshold
                double correlation = snp_dosages[i].ref_dosages.dot(snp_dosages[j].ref_dosages); // calculate the correlation between their reference dosages

                triplet.push_back(Tlet(i,j,correlation)); // push back a triplet with lower triangular value correlation

//                triplet.push_back(Tlet(j,i,correlation)); // push back a triplet with upper triangular value correlation
            }
        }
    }

    zofstream output(outputfile);

    forc(i, snp_dosages){
        output << snp_dosages[i].posclass << '\t' << i << '\n';
    }

    for(auto& v: triplet){
        output << v.col() << '\t' << v.row() << '\t' << v.value() << '\n';
    }
}
