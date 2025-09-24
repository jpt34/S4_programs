#include <iostream>
#include <random>
#include <experimental/string_view>
#include <boost/program_options.hpp>
#include <boost/math/distributions/inverse_gaussian.hpp>
#include <Eigen/Dense>
#include <Eigen/Sparse>

#include <MachineTimer.h>
#include <jlibrary.h>
#include <Posclass.h>

#include "bgen/parser.hpp"
#include "bgen/IndexQuery.hpp"
#include "Mlogit.h"
#include "routines.h"
#include "SNPcontributions.h"
#include "MultiPGS.h"
#include "PGSregression.h"
#include "LowSignificanceCorrection.h"


using namespace std;
using namespace Eigen;

namespace po = boost::program_options;
using  SpMat = SparseMatrix<double>;


// Input parameters to run program

struct Parameters
{
    double maf;                 // minimum maf
    double threshold;
    double test_threshold;
    double minr;
    double maxr;
    double maxdistance;
    double r2impute;
    string chromosome;
    bool correct_r2;
    double base_num;
    bool use_multiple_r2;
    double zchange;
    string includefile;
    string method;
    string refsnpfile;
    string genmapfile;
    double sparse;
    double max_individual_r;
    bool novarqc;
    double corr_damp;
    bool logresults;
    vector<double> qc;
    double indiv_change;
    double cojo_threshold;
    bool use_effect_threshold;
    double subset_threshold;
    string hybrid_file;
    double hybrid_threshold;
    int ancestry = 0;
    vector<double> bounds;
    double minindr;
    string excludesnps;
    double regval = 0.0;
    bool functional_scores=false;

    double CalcMaxPsi(double zscore) const
    {
        if(regval > 1e-10) return 1.0/regval;

        if(bounds.empty()) return 1e50;

        double result = bounds[0];

        for(int i=2;i<int(bounds.size());i+=2){
            if(abs(zscore) < bounds[i-1]){
                result = bounds[i];
                break;
            }
        }

        return result;
    };

    bool IsCojo() const
    {
        return (method=="cojo" || method=="reduced");
    }
};


struct SumStatsParameters
{
    double genetic_distance;
    double corr_damp;
};


// Get positions in bgen file of SNPs of interest

tuple<map<string, string>, map<Posclass, long>> GetPositions(const string& inputfile, const vector<string>& analysis_files, const Parameters& parameters)
{
// Get SNPs of interest

    struct SNPdetails
    {
        int position;
    };

// read positions and convert

    map<Posclass, string> snpconv;

    if(parameters.refsnpfile!=""){
        zifstream refinput(parameters.refsnpfile);

        SplitString splitstr;

        Posclass posclass;

        while(refinput >> splitstr){
            if(splitstr.size()>=7){
                try {
                    posclass = Posclass(0, stoi(splitstr[3]), splitstr[5], splitstr[6]);

                    snpconv[posclass] = splitstr[1];
                } catch (...) {
//                    cerr << "Error reading line containing\n" << splitstr.str() << "\n";
                }
            }
        }
    }

    set<PosclassUnordered> hybrid_snps;

    if(parameters.hybrid_file!=""){
        zifstream hybrid_input(parameters.hybrid_file);

        SplitString splitstr;

        PosclassUnordered posclass;

        while(hybrid_input >> splitstr){
            if(splitstr.size()>=5){
                try {
                    posclass = PosclassUnordered(ConvertChrStr(splitstr[1]), stoi(splitstr[2]), splitstr[3], splitstr[4]);

                    hybrid_snps.insert(posclass);
                } catch (...) {
                    //                    cerr << "Error reading line containing\n" << splitstr.str() << "\n";
                }
            }
        }
    }

    BgenParser bgenParser(inputfile);

    map<Posclass,string> snp_list;

    set<string> ids_set;

    bool exclude_amb = (string::npos != parameters.excludesnps.find("amb"));

    bool exclude_indel = (string::npos != parameters.excludesnps.find("indel"));

    auto IsAmb = [&](const string& ref, const string& alt){
        const static std::map<std::string, int> conv{{"A", 1}, {"T",1}, {"C",2}, {"G",2}};
        try{
            if(conv.at(ref)==conv.at(alt)){
                return true;
            }
            else{
                return false;
            }
        }
        catch(...){
            return false;
        }
    };

    for(auto v:analysis_files){
        zifstream input(v);

        SplitString splitstr;

        bool firstline = true;

        while(input >> splitstr){
            try{
                if(splitstr.size()<9) throw 1;

                if(firstline){
                    firstline = false;

                    if(!std::all_of(splitstr[1].begin(), splitstr[1].end(), ::isdigit)) continue;
                }

                double threshold = (parameters.method=="cojo"?parameters.cojo_threshold:parameters.threshold);

                int in_hybrid = 0;

                if(parameters.hybrid_file!=""){
                    in_hybrid = hybrid_snps.count(PosclassUnordered(ConvertChrStr(parameters.chromosome), stoi(splitstr[1]),splitstr[2], splitstr[3]));

//                    if(!in_hybrid && (splitstr[2].size()!=1 || splitstr[3].size()!=1)){
//                        in_hybrid = hybrid_snps.count(Posclass(ConvertChrStr(parameters.chromosome), stoi(splitstr[1]),splitstr[3], splitstr[2]));
//                    }
                }

                bool excludesnp = (exclude_amb && IsAmb(splitstr[2], splitstr[3])) || (exclude_indel && (splitstr[2].size()!=1 || splitstr[3].size()!=1));

                if(((!excludesnp && stold(splitstr[8])< threshold && abs(stod(splitstr[4])-0.5) < 0.5-parameters.maf) || in_hybrid) && stod(splitstr[5])>=parameters.r2impute){
                    double eaf = stod(splitstr[4]);

                    double fscore = parameters.functional_scores && splitstr.size()>9?
                                         stod(splitstr[9]):
                                         1.0;

                    double r2_ss = fscore*(parameters.IsCojo()?1.0:1/(eaf*(1-eaf)*Sqr(stod(splitstr[7]))*parameters.base_num));

                    if(stold(splitstr[8])/r2_ss < threshold || (in_hybrid && stold(splitstr[8])/r2_ss < parameters.hybrid_threshold)){
                        if(parameters.refsnpfile==""){
                            ids_set.insert(splitstr[0]);
                        }
                        else{
                            auto itor_snpconv = GetPosItor(snpconv, Posclass(0,stoi(splitstr[1]),splitstr[2], splitstr[3]));

                            if(itor_snpconv!=snpconv.end()){
//                                ids_set.insert(snpconv[Posclass(0,stoi(splitstr[1]),splitstr[2], splitstr[3])]);
                                ids_set.insert(itor_snpconv->second);
                                snp_list[Posclass(0,stoi(splitstr[1]),splitstr[3], splitstr[2])] = splitstr[0];
                            }
                        }
                    }
                }
            }
            catch(...){
                cerr << "Error reading line containing:\n" << splitstr.str() << "\n";
            }
        }

// only choose variants for ancestry analysis based on first file

        if(parameters.ancestry){
            break;
        }
    }

    vector<string> ids(ids_set.begin(), ids_set.end());

    genfile::bgen::SqliteIndexQuery query(inputfile +".bgi");

    query.include_rsids(ids);

    query.initialise();

    tuple<map<string, string>, map<Posclass, long>> result;

    for(size_t i=0;i<query.number_of_variants();i++){

        auto fileposition = query.locate_variant(i);

        bgenParser.JumpToPosition(fileposition.first);

        std::string chromosome ;
        uint32_t position ;
        std::string rsid ;
        std::vector< std::string > alleles ;
        std::vector< std::vector< double > > probs ;

        bgenParser.read_variant( &chromosome, &position, &rsid, &alleles );

        get<1>(result)[Posclass(ConvertChrStr(chromosome), position, alleles[0], alleles[1])] = fileposition.first;

        if(parameters.refsnpfile==""){
            get<0>(result)[rsid] = rsid;
        }
        else{
            get<0>(result)[snp_list[Posclass(0, position, alleles[0], alleles[1])]] = rsid;
        }

    }

    return result;
}


// contains statistics about SNPs together with their names, positions etc.

struct SNPdetails
{
    string rsnum;

    string chromosome;

    int position;

    double genpos;

    string effect;

    string baseline;

    double eaf;         // Effect allele frequency

    double r2;

    double ref_r2;

    double ss_r2;

    long double minpvalue = 1.0;        // minpvalue is long double so as to be able to use very small p values

    double functional_score = 1.0;

    vector<tuple<double,double>> stats;     // this contains beta estimates and their standard errors for the SNP.
                                            // It is a vector so that extra histotypes can be analysed.

    // Work out the best individual statistics if more that one histotype

    double CalcBestIndividualStats() const{
        double maxstat=0;

        for(auto& v:stats){
            double absstat=abs(get<0>(v)/get<1>(v));

            if(absstat>maxstat) maxstat=absstat;
        }

        return maxstat;
    }
};


struct SNPvector
{
    SNPdetails snpinfo;

    VectorXd snpvec;
};


auto First_abs = [](auto x){
    return abs(get<0>(x[0])/get<1>(x[0]));
};

auto Max_abs = [](auto x){
    auto v = *max_element(begin(x), end(x),[](auto x, auto y){return abs(get<0>(x)/get<1>(x))<abs(get<0>(y)/get<1>(y));});

    return abs(get<0>(v)/get<1>(v));
};

auto Max_index = [](auto x){
    auto v = max_element(begin(x), end(x),[](auto x, auto y){return abs(get<0>(x)/get<1>(x))<abs(get<0>(y)/get<1>(y));});

    return (v - begin(x)) + 1;
};


auto CalcBestStats = [](auto x){
    VectorXd beststat = VectorXd::Zero(get<0>(x[0]).size());

    VectorXd bestindex = VectorXd::Zero(get<0>(x[0]).size());

    int num = 0;

    for(auto& v: x){
        VectorXd z = get<0>(v).cwiseQuotient(get<1>(v)).cwiseAbs();

        ++num;

        forc(i, beststat){
            if(z[i]>beststat[i]){
                bestindex[i] = num;
            }
        }

        beststat = beststat.cwiseMax(z);
    }

    return make_tuple(beststat, bestindex);
};


struct SNPgrouping
{
    vector<SNPvector> snps;

    VectorXd best_statistics;

    VectorXd best_index;

    MatrixXd corr;

    MatrixXd corr_triangular;

    // return correlation matrix when adding a new variable contain correlations between it and the other variables

    MatrixXd CalcCorr(const VectorXd& newcorr, int index) const
    {
        int nrows = newcorr.size() + (index == newcorr.size());

        MatrixXd z(nrows, nrows);

        if(index==newcorr.size()){
            z << corr, newcorr, newcorr.transpose(), 1;
        }else{
           z = corr;
           forl(i,newcorr.size()-1){
               int non_index = i + (i>=index);
               z(index, non_index) = z(non_index, index) = newcorr[non_index];
           }
        }

        return z;
    }

    // calculate conditional statistics based on summary statistics and the correlation matrix

    tuple<vector<tuple<VectorXd,VectorXd>>,MatrixXd, double, double> CalcStats(const SNPdetails& snpdetails, const VectorXd& newcorr, int index, bool use_only_first) const
    {
        MatrixXd z = CalcCorr(newcorr, index);

        int nrows = newcorr.size() + (index == newcorr.size());

        VectorXd b(nrows);

        VectorXd se(nrows);

        vector<tuple<VectorXd,VectorXd>> result;

        double maxchange = 0;           // maxchange is the maximum increase in z-score between a condtional statistic over the individual test statistic.
                                        // If this is above a threshold then the SNP shouldn't be added to the model as something is probably wrong with the correlation structure.

        double indiv_change  = 0;

        forl(i, snpdetails.stats.size()){
            forc(j, newcorr){
                b[j] = get<0>(snps[j].snpinfo.stats[i]);

                se[j] = get<1>(snps[j].snpinfo.stats[i]);
            }

            b[index] = get<0>(snpdetails.stats[i]);

            se[index] = get<1>(snpdetails.stats[i]);

            DiagonalMatrix<double, -1> semat = se.asDiagonal().inverse();

            MatrixXd bmat = semat * z * semat;

//            MatrixXd variance = bmat.inverse();

            MatrixXd variance = bmat.llt().solve(MatrixXd::Identity(nrows, nrows));

            VectorXd beta = variance*semat*semat*b;

            VectorXd se_result = variance.diagonal().cwiseSqrt();

            forc(j, beta){
                if((abs(beta[j])-abs(b[j]))/se[j] > maxchange) maxchange = (abs(beta[j])-abs(b[j]))/se[j];

//                if(abs(b[j])>maxchange & b[j]*beta[j] < 0) maxchange = abs(b[j]);
            }

            indiv_change = abs(beta[index])-abs(b[index])/se[index];

            result.emplace_back(beta, se_result);

            if(use_only_first) break;
        }

        return make_tuple(result, z, maxchange, indiv_change);
    }


    // calculate conditional statistics based on summary statistics and the correlation matrix

    tuple<vector<tuple<VectorXd,VectorXd>>,MatrixXd, double, double> CalcStatsNoVariance(const SNPdetails& snpdetails, const VectorXd& newcorr, int index, const Parameters& parameters) const
    {
        MatrixXd z = CalcCorr(newcorr, index);

        double zstat = get<0>(snpdetails.stats[0])/get<1>(snpdetails.stats[0]);

        z(index, index) = 1.0 + 1.0/parameters.CalcMaxPsi(zstat);

        int nrows = newcorr.size() + (index == newcorr.size());

        double r2 = 1/(snpdetails.eaf*(1-snpdetails.eaf)*Sqr(get<1>(snpdetails.stats[0]))*parameters.base_num);

        if(r2 > 1) r2 = 1;

        VectorXd b(nrows);

        VectorXd se(nrows);

        vector<tuple<VectorXd,VectorXd>> result;

        double maxchange = 0;           // maxchange is the maximum increase in z-score between a condtional statistic over the individual test statistic.
                                        // If this is above a threshold then the SNP shouldn't be added to the model as something is probably wrong with the correlation structure.

        LLT<MatrixXd> zllt(z);

        double indiv_change = 0;

        forl(i, snpdetails.stats.size()){
            forc(j, newcorr){
                b[j] = get<0>(snps[j].snpinfo.stats[i]);

                se[j] = get<1>(snps[j].snpinfo.stats[i]);
            }

            b[index] = get<0>(snpdetails.stats[i]);

            se[index] = get<1>(snpdetails.stats[i]);

            VectorXd zscores = b.cwiseQuotient(se);

            VectorXd beta = zllt.solve(zscores);

            forc(j, beta){
                if(abs(beta[j])-abs(zscores[j]) > maxchange) maxchange = (abs(beta[j])-abs(zscores[j]));

//                if(abs(beta[j])>maxchange & zscores[j]*beta[j] < 0) maxchange = abs(beta[j]);
            }

            result.emplace_back(beta.cwiseProduct(se), se);

            if(parameters.ancestry) break;  // only analyse the first value if data is from several ancestries
        }

        return make_tuple(result, z, maxchange, indiv_change);
    }

    // Only calculate multiple r2 for correlations about subset threshold
    // This is done to avoid issues with multiple r2 taking into account every SNP in a long list

    double CalculateSubsetR2(const VectorXd& newcorr, double subset_threshold)
    {
        vector<int> indices;

        forc(i, newcorr){
            if(abs(newcorr[i]) >= subset_threshold){
                indices.push_back(i);
            }
        }

        VectorXd subcorr = newcorr(indices);

        MatrixXd submat = corr(indices, indices);

        LLT<MatrixXd> llt(submat);

        MatrixXd transform = llt.matrixL();

        return transform.triangularView<Lower>().solve(subcorr).norm();
    }

    // calculate conditional statistics based on summary statistics, the correlation matrix and restricting analysis to correlations above subset_threshold

    tuple<vector<tuple<VectorXd,VectorXd>>,MatrixXd, double, double> CalcSubsetStatsNoVariance(const SNPdetails& snpdetails, const VectorXd& newcorr, const Parameters& parameters) const
    {
        vector<int> indices;

        forc(i, newcorr){
            if(abs(newcorr[i]) >= parameters.subset_threshold){
                indices.push_back(i);
            }
        }

        VectorXd subcorr = newcorr(indices);

        MatrixXd submat = corr(indices, indices);

        int nrows = subcorr.size() + 1;

        MatrixXd z(nrows, nrows);

        z << submat, subcorr, subcorr.transpose(), 1;

        VectorXd b(nrows);

        VectorXd se(nrows);

        vector<tuple<VectorXd,VectorXd>> result;

        double maxchange = 0;           // maxchange is the maximum increase in z-score between a condtional statistic over the individual test statistic.
                                        // If this is above a threshold then the SNP shouldn't be added to the model as something is probably wrong with the correlation structure.

        LLT<MatrixXd> zllt(z);

        int index = nrows-1;

        double indiv_change = 0;

        forl(i, snpdetails.stats.size()){
            forc(j, indices){
                b[j] = get<0>(snps[indices[j]].snpinfo.stats[i]);

                se[j] = get<1>(snps[indices[j]].snpinfo.stats[i]);
            }

            b[index] = get<0>(snpdetails.stats[i]);

            se[index] = get<1>(snpdetails.stats[i]);

            VectorXd zscores = b.cwiseQuotient(se);

            VectorXd beta = zllt.solve(zscores);

            forc(j, beta){
                if(abs(beta[j])-abs(zscores[j]) > maxchange) maxchange = (abs(beta[j])-abs(zscores[j]));
            }

            indiv_change = abs(beta[index])-abs(zscores[index]);

            result.emplace_back(beta.cwiseProduct(se), se);

            if(parameters.ancestry) break;
        }

        return make_tuple(result, z, maxchange, indiv_change);
    }


    // This function adds the SNP into the model if it posses a few tests
    // there are three methods:
    // s4 which is the main s4 algorithm to select SNPs to then apply to the shrinkage algorithm
    // cojo which does a conditional analysis based on the correlation stucture. This is mainly used for finding GWAS significant SNPs.
    // pca which is to select uncorrelated SNPs for a principal component analysis
    // zthreshold is only used if a cojo analysis is used.

    void UpdateIfSignificant(const SNPdetails& v, double zthreshold, const VectorXd& temp, const Parameters& parameters)
    {
        VectorXd newcorr(corr.rows());

        // ncomp is the expected variance for the test snp compared to what would be the expected variance from the imputed value
        // this is set to a maximum of one in case the standard error is smaller than what would be expected
        // ncomp is roughly equivalent to the imputation accuracy.

        double ncomp = min(1.0, 1.0/(v.ref_r2*v.eaf*(1-v.eaf)*Sqr(get<1>(v.stats[0]))*parameters.base_num));

        if(ncomp*v.ref_r2 < parameters.r2impute && parameters.correct_r2) return;

        forc(i,newcorr){

            // ncomp2 is the expected variance for each snp already included

            double ncomp2 = min(1.0, 1.0/(snps[i].snpinfo.ref_r2*snps[i].snpinfo.eaf*(1-snps[i].snpinfo.eaf)*Sqr(get<1>(snps[i].snpinfo.stats[0]))*parameters.base_num));

            // the correlation between the variables is reduced if one is less well imputed than the other

            if(parameters.correct_r2){
                newcorr[i] = parameters.corr_damp*snps[i].snpvec.dot(temp)*sqrt(min(ncomp, ncomp2)/max(ncomp, ncomp2));
            }
            else{
                newcorr[i] = snps[i].snpvec.dot(temp);
            }
        }

        double r2 = 1/(v.eaf*(1-v.eaf)*Sqr(get<1>(v.stats[0]))*parameters.base_num);

        if(r2 > 1) r2 = 1;

        bool in_hybrid = parameters.hybrid_file!="" && v.minpvalue/r2 > parameters.threshold;

        double r2_min = newcorr.cwiseAbs().maxCoeff();

        double r2_value;

        if(in_hybrid){
            r2_value = 0;
        }
        else{
            r2_value = ((parameters.use_multiple_r2 && r2_min < parameters.maxr) ? parameters.subset_threshold>0 ? CalculateSubsetR2(newcorr, parameters.subset_threshold)
                                                                                        : corr_triangular.triangularView<Lower>().solve(newcorr).norm()
                                                                                        : r2_min); // don't calculate if r2_min is above threshold
        }


//        double r2_compare = CalculateSubsetR2(newcorr, parameters.subset_threshold);

//        cout << r2_value << ' ' << r2_compare << ' ' << r2_value-r2_compare << '\n';

//        if(r2_value < r2_min){
//            cout << r2_value << ' ' << r2_min << ' ' << r2_value - r2_min << '\n';
//        }

        tuple<vector<tuple<VectorXd,VectorXd>>,MatrixXd,double, double> result;

//        cout << "Test1" << endl;

//        cout << parameters.use_multiple_r2 << endl;

//        cout << (parameters.use_multiple_r2 && r2_min < parameters.maxr) << endl;

//        cout << parameters.subset_threshold << endl;

//        cout << r2_min << endl;

//        cout << corr_triangular.triangularView<Lower>().solve(newcorr).norm() << endl;

//        cout << corr_triangular << endl << endl;

//        cout << newcorr << endl << endl;

//        cout << r2_value << '\t' << corr_triangular.triangularView<Lower>().solve(newcorr).norm() << '\t' << r2_min << '\t' << parameters.maxr << endl;

        if(r2_value < parameters.maxr){

            // if just choosing snps for a principal component analysis then add any snp below a correlation threshold

            if(parameters.method=="pca"){
                snps.push_back({v, temp});
                corr = CalcCorr(newcorr, newcorr.size());
                LLT<MatrixXd> llt(corr);
                corr_triangular = llt.matrixL();
                return;
            }

            if(parameters.IsCojo()){
                result = CalcStats(v, newcorr, newcorr.size(), parameters.ancestry);
            }
//            else if(parameters.subset_threshold > 0){
//                result = CalcSubsetStatsNoVariance(v, newcorr, parameters.subset_threshold);
//            }
            else{
                result = CalcStatsNoVariance(v, newcorr, newcorr.size(), parameters);
            }

            auto [beststat, bestindex] = CalcBestStats(get<0>(result));

//            auto result2 = CalcSubsetStatsNoVariance(v, newcorr, 0.05);

//            auto [beststat2, bestindex2] = CalcBestStats(get<0>(result2));

//            cout << get<2>(result) << ' ' << get<2>(result2) << '\n';

            // This is run if the maximum change in the test statisic is below a threshold (suggesting no qc problems)
            // and the method isn't cojo or the method is cojo and the test statistic is above the threshold.

//            double zcoeff = get<0>(get<0>(result)[0])[newcorr.size()]/get<1>(get<0>(result)[0])[newcorr.size()];

//            double zorig =

//            cout << get<2>(result) << '\n';

            if((!parameters.IsCojo() || beststat[newcorr.size()] > zthreshold) && get<2>(result)<parameters.zchange && get<3>(result)<parameters.indiv_change){
//                cout << "Passed " << snps.size() << "\n";
                snps.push_back({v, temp});
//                if(parameters.subset_threshold>0){
//                    best_statistics = (VectorXd(best_statistics.size()+1) << best_statistics, abs(get<0>(v.stats[0])/ get<1>(v.stats[0]))).finished();
//                    best_index = (VectorXd(best_index.size()+1) << best_index, 1).finished();
//                    corr = (MatrixXd(corr.rows()+1, corr.cols()+1) << corr, newcorr, newcorr.transpose(), 1).finished();
//                    return;
//                }
//                else{
                    best_statistics = beststat;
                    best_index = bestindex;
                    corr = get<1>(result);
                    LLT<MatrixXd> llt(corr);
                    corr_triangular = llt.matrixL();
                    return;
//                }
            }
            else{
//                cout << "Failed " << snps.size() << "\n";
            }
        }

        // if the method is cojo then test all statistics where one of the SNPs already in the model is replaced by the tested SNP.

        if(parameters.method=="cojo"){
            for(int i=newcorr.size()-1;i>=0;i--){

                VectorXd testmax = newcorr;

                testmax[i] = 0;

                double r2_value2 = [&]{
                    if(parameters.use_multiple_r2){
                        MatrixXd corr_i = corr;

                        forc(j, newcorr){
                            corr_i(i,j)=corr_i(j,i) = 0;
                        }

                        corr_i(i,i) = 1;

                        return sqrt(testmax.transpose()*corr_i.inverse()*testmax);
                    }
                    else{
                        return testmax.cwiseAbs().maxCoeff();
                    }
                }();

                if(r2_value2 < parameters.maxr){
                    result = CalcStats(v, newcorr, i, parameters.ancestry);

                    auto [beststat,bestindex] = CalcBestStats(get<0>(result));

                    if(beststat[i] > best_statistics[i] + 0.5*(i==0 || parameters.method=="sc") && get<2>(result)<parameters.zchange){
                        snps[i] = SNPvector{v, temp};
                        best_statistics = beststat;
                        best_index = bestindex;
                        corr = get<1>(result);
                        LLT<MatrixXd> llt(corr);
                        corr_triangular = llt.matrixL();
                        return;
                    }
                }
            }
        }
    }


    // works out conditional statistics which can be printed if necessary

    vector<tuple<VectorXd,VectorXd>> CalcConditionalStats(bool use_only_first) const
    {
        int nrows = corr.rows();

        VectorXd b(nrows);      // contains the beta estimates

        VectorXd se(nrows);     // contains the standard errors

        vector<tuple<VectorXd,VectorXd>> result;

        forl(i, snps[0].snpinfo.stats.size()){
            forl(j, nrows){
                b[j] = get<0>(snps[j].snpinfo.stats[i]);

                se[j] = get<1>(snps[j].snpinfo.stats[i]);
            }

            MatrixXd semat = se.asDiagonal().inverse();

            MatrixXd bmat = semat * corr * semat;

            MatrixXd variance = bmat.inverse();

            VectorXd beta = variance*semat*semat*b;

            VectorXd se_result = variance.diagonal().cwiseSqrt();

            result.emplace_back(beta, se_result);   // beta is the conditional beta estimates and se_results is the conditional standard error

            if(use_only_first) break;
        }

        return result;
    }
};


// Get SNP statisitics of snps that are in the reference file and passed the p-value threshold

map<string, SNPdetails> GetSNPstatistics(const map<string, string>& snpfilepositions, const vector<string>& analysis_files, const Parameters& parameters)
{
    map<string, SNPdetails> results;

    forc(i, analysis_files){
        string origsnpstr, snpstr;
        SplitString remainingline;
        zifstream input(analysis_files[i]);
        while(input >> origsnpstr){
            if(auto itor = snpfilepositions.find(origsnpstr);itor!=snpfilepositions.end()){     // snp found in map
                snpstr = itor->second;

                input >> remainingline;

                if(remainingline.size()<9){     // check line has enough fields
                    continue;
                }

                if(results[origsnpstr].stats.empty()){
                    results[origsnpstr].stats = vector<tuple<double,double>>(analysis_files.size(),{0.0,0.0});

                    tie(results[origsnpstr].rsnum, results[origsnpstr].position, results[origsnpstr].effect, results[origsnpstr].baseline, results[origsnpstr].eaf, results[origsnpstr].r2) = make_tuple(snpstr, stoi(remainingline[1]), remainingline[2], remainingline[3],stod(remainingline[4]), stod(remainingline[5]));
                }

                struct TempStr{
                    string ref;
                    string alt;
                };

                int factor = (remainingline[2] == results[origsnpstr].effect) + (remainingline[3] == results[origsnpstr].baseline) - 1;

                int factor2 = (remainingline[2] == results[origsnpstr].baseline) + (remainingline[3]==results[origsnpstr].effect) - 1;

                if(factor==0 || factor+factor2!=0 || stoi(remainingline[1]) != results[origsnpstr].position){
                    cerr << "Allele or position mismatch for " << origsnpstr << "\n";
                    exit(1);
                }

                results[origsnpstr].stats[i] = {factor*stod(remainingline[6]),stod(remainingline[7])};

                if(stold(remainingline[8])< results[origsnpstr].minpvalue && (i==0 || !parameters.ancestry)){
                    results[origsnpstr].minpvalue = stold(remainingline[8]);
                }

                if(parameters.functional_scores && remainingline.size()>9){
                    results[origsnpstr].functional_score = stod(remainingline[9]);
                }
            }else{
                input >> skipline;
            }
        }
    };

    return results;
}


// Get iterator for positions using SNPdetails and a map

template<class T>
decltype(auto) GetPosIterator(const T& snp_positions, const SNPdetails& v)
{
    auto p = snp_positions.find(Posclass(ConvertChrStr(v.chromosome), v.position, v.baseline, v.effect));

    if(p==snp_positions.end()){
        p = snp_positions.find(Posclass(ConvertChrStr(v.chromosome), v.position, v.effect, v.baseline));
    }

    return p;
}


// Define a function that takes in an imputed class, residuals, and SNP details
VectorXd CreateAncestryNormalised(const Imputeclass& imputeclass, const vector<double>& residuals, const SNPdetails& snpdetails)
{
    // Declare a map that will store the indices of the input data points that belong to each ancestry class
    map<int, vector<int>> ancestry_indices;

    // Declare a vector to store the ancestry-normalized values for the input data
    VectorXd result(imputeclass.size());

    // Loop through each residual value and add its index to the corresponding ancestry class in the map
    // Also, copy the corresponding imputed class value to the result vector
    forc(i, residuals){
        // Round the current residual value to the nearest integer and use the result as the key to access the appropriate ancestry class in the map
        ancestry_indices[int(residuals[i]+0.5)].push_back(i);
        // Copy the corresponding imputed class value to the result vector
        result[i] = imputeclass[i];
    }
//{
//    ofstream output1("test1.txt");

//    output1 << result << '\n';
//}
//    cout << result << "\n\n";

    // Loop through each ancestry class (except the first one) and normalize the values in the result vector
    forv(i, 1, snpdetails.stats.size()){
        // Subtract the mean of the result values for the current ancestry class from all the values in the result vector for this ancestry class
        result(ancestry_indices[i]).array() -= result(ancestry_indices[i]).mean();
        // Calculate the norm of the values in the result vector for this ancestry class
        double norm = result(ancestry_indices[i]).norm();
        // Divide the result vector for this ancestry class by the appropriate scaling factor

        if(get<1>(snpdetails.stats[i]) <= 0.0){     // not analysed in dataset
           result(ancestry_indices[i]) *= 0.0;
        }
        else if(norm!=0){                           // if norm is 0 all values are 0 so can't be refactored
            result(ancestry_indices[i]) /= get<1>(snpdetails.stats[i]) * norm;
        }
    }

//    VectorXd test = result.normalized();
//{
//    ofstream output2("test2.txt");

//    output2 << result.normalized() << '\n';
//}
//    cout << result.normalized() << '\n';

//    exit(0);

    // Return the normalized result vector
    return result.normalized();
}


// This is the main routine to generate a list of SNPs to analyse
// It will print out the SNPs in order at which the the top SNP was put into the model

void SeveralCorrelatedAnalyses(const string& inputfile, const vector<string>& analysis_files, const Parameters& parameters)
{
    // Get SNP names and position in bgen files

    auto [snp_name_mappings, snp_positions] = GetPositions(inputfile, analysis_files, parameters);

    // Get SNP statistics for snps generate by GetPositions routine

    auto snp_statistics = GetSNPstatistics(snp_name_mappings, analysis_files, parameters);

    vector<SNPdetails> sorted_snps(snp_statistics.size());

    transform(begin(snp_statistics), end(snp_statistics), begin(sorted_snps), [](auto x){return x.second;});

    // stable sort is used to make top statistic generation the same no matter the p-value threshold
    // also use the minimum p-value when the maximum absolute statistics are the same

    stable_sort(begin(sorted_snps), end(sorted_snps), [&](auto& x, auto& y){
        double r2_x = 1/(x.eaf*(1-x.eaf)*Sqr(get<1>(x.stats[0]))*parameters.base_num);

        if(r2_x>1) r2_x = 1;

        double r2_y = 1/(y.eaf*(1-y.eaf)*Sqr(get<1>(y.stats[0]))*parameters.base_num);

        if(r2_y>1) r2_y = 1;

        if(parameters.functional_scores){
            r2_x *= x.functional_score;
            r2_y *= y.functional_score;
        }

        if(parameters.IsCojo()) r2_x = r2_y = 1.0;

        if(parameters.ancestry){
            return make_pair(-x.minpvalue/r2_x, First_abs(x.stats)) > make_pair(-y.minpvalue/r2_y, First_abs(y.stats));
        }
        else{
            return make_pair(-x.minpvalue/r2_x, Max_abs(x.stats)) > make_pair(-y.minpvalue/r2_y, Max_abs(y.stats));
        }
    });

    vector<SNPgrouping> snpgroups;

    BgenParser bgenParser(inputfile);       // this will read data from the bgen file

    Imputeclass imputeclass;                // this will contain data from the bgen file

    boost::math::normal gaussian;

    double zthreshold = (parameters.test_threshold<=0?-quantile(gaussian, parameters.threshold/2):-quantile(gaussian, parameters.test_threshold/2));       // zthreshold is only used for a cojo analysis looking for secondary hits that may be at a different significance level than the primary hits

    struct FlipCheck{
        string ref;
        string alt;
    };

    vector<double> residuals = (parameters.includefile==""?vector<double>():FileRead<vector<double>>(parameters.includefile));

    if(parameters.includefile!=""){
        if(!ifstream(parameters.includefile)){
            exit(1);
        }
    }

    vector<int> include;

    vector<double> residuals_filtered;

    bool residuals_used = false;

    for(auto& v:residuals){
        include.push_back(double(v!=0));

        if(v!=0){
            residuals_filtered.push_back(v);

            if(v!=1) residuals_used = true;
        }
    }

    VectorXd temp(parameters.includefile==""?bgenParser.number_of_samples():residuals_filtered.size());     // this contains the data from the present SNP that may be added to the snp list

    auto genmap = ReadGeneticMap(parameters.genmapfile);

    for(auto& v:sorted_snps){
        // set the chromosome for mapping
        v.chromosome = parameters.chromosome;

        // cycle through sorted_snps where the most significant value is first
        auto p = GetPosIterator(snp_positions, v);

        if(p==snp_positions.end()){
            continue;
        }

        bgenParser.JumpToPosition(p->second);

        if(parameters.logresults) clog << v.minpvalue << '\n';

        if(parameters.includefile==""){
            bgenParser.ReadAllImputeclass(imputeclass);
        }
        else{
            bgenParser.ReadImputeclass(imputeclass, include);
        }

        int flip = MatchAlleles(imputeclass, FlipCheck{v.baseline, v.effect});

        if(imputeclass.eaf == 0 || imputeclass.eaf==1) continue;     // don't use data if monomorphic in reference file

        if(parameters.ancestry){            
            temp = flip * CreateAncestryNormalised(imputeclass, residuals_filtered, v);
        }
        else{
            double mean = 2.0*imputeclass.eaf;

            if(mean==0 || mean==2) continue;        // don't use data if monomorphic in reference file

            if(residuals_used){
                mean = 0;
                double total_weight = 0;
                forc(i, imputeclass){
                    if(imputeclass[i]>-0.5){
                        mean += imputeclass[i]*residuals_filtered[i];
                        total_weight += residuals_filtered[i];
                    }
                }
                mean /= total_weight;

                forc(i, temp){
                    if(imputeclass[i]>-0.5){
                        temp[i] = flip*residuals_filtered[i]*(imputeclass[i] - mean);
                    }
                    else{
                        temp[i]=0;
                    }
                }
            }
            else{
                forc(i, temp){
                    if(imputeclass[i]>-0.5){
                        temp[i] = flip*(imputeclass[i] - mean);
                    }
                    else{
                        temp[i]=0;
                    }
                }
            }

    //        double sd_val = sqrt(temp.squaredNorm()/temp.size());

            temp /= temp.norm();
        }

        v.ref_r2 = imputeclass.r2;      // ref_r2 is the r2 in the reference file rather than the statistics file

        double ss_r2 = 1/(v.eaf*(1-v.eaf)*Sqr(get<1>(v.stats[0]))*parameters.base_num);

        bool is_bad = [&]{
//            if(parameters.novarqc) return false;

//            double sd_ss = 1.0/sqrt(0.5*Sqr(get<1>(v.stats[0]))*parameters.base_num);

//            return bool(sd_ss < (0.5 * sd_val) | sd_ss > (sd_val + 0.1) | sd_ss < 0.1 | sd_val < 0.05);

            bool result = false;

            try {
                result |= v.ref_r2 < parameters.qc.at(0);
                result |= ss_r2 > parameters.qc.at(1);
                result |= ss_r2 > parameters.qc.at(2)*v.ref_r2;
                double maf_ss = min(v.eaf,1-v.eaf);
                double maf_ref = min(imputeclass.eaf, 1-imputeclass.eaf);
                result |= maf_ss*v.ref_r2 < parameters.qc.at(3);
                result |= Sqr(maf_ss-maf_ref)/(maf_ss+maf_ref)/(2-maf_ss-maf_ref) > parameters.qc.at(4);

            } catch (...) {
            }

            return (result);
        }();


        if((ss_r2 < parameters.r2impute && parameters.correct_r2) || is_bad) continue;        // ignore if r2 looks too low in summary statistics file
                                                                                                    // this also include cases when studies are included as the conditional analysis still becomes difficult
        double bestr2 = 0;
        double best_secondary_r2 = 0;
        int index = -1;
        int secondary_index = 0;

        // Also compare the effect of the relevant SNP on the candidate SNP.
        // Very significant SNPs will have a stronger effect on the z-score.

        double besteffect = 0;
        int effect_index = -1;

        v.genpos = CalcGeneticPosition(genmap, v.position);

        forc(i, snpgroups){
            double distance = abs(snpgroups[i].snps[0].snpinfo.genpos - v.genpos);

            double r2 = (distance < parameters.maxdistance ? snpgroups[i].snps[0].snpvec.dot(temp) : 0);    // if distance is greater than max distance assume 0 otherwise calculate correlation

            if(parameters.use_effect_threshold && abs(r2) < parameters.minindr){
                r2 = 0.0;
            }

            double effect = abs(r2) * (parameters.ancestry ? First_abs(snpgroups[i].snps[0].snpinfo.stats) : Max_abs(snpgroups[i].snps[0].snpinfo.stats));  // effect of stat on SNP statistic - effect * r2

            if(abs(r2) > abs(bestr2)){
                index = i;
                bestr2 = r2;
            }

            if(effect > besteffect){
                effect_index = i;
                besteffect = effect;
            }

            forv(j, 1, snpgroups[i].snps.size()){
//                double distance2 = abs(snpgroups[i].snps[j].snpinfo.genpos - v.genpos);

                double r2 = (distance < parameters.maxdistance ? snpgroups[i].snps[j].snpvec.dot(temp) : 0);

                if(parameters.use_effect_threshold && abs(r2) < parameters.minindr){
                    r2 = 0.0;
                }

                double effect = abs(r2) * (parameters.ancestry ? First_abs(snpgroups[i].snps[0].snpinfo.stats) : Max_abs(snpgroups[i].snps[0].snpinfo.stats));  // effect of stat on SNP statistic - effect * r2

                if(abs(r2) > abs(best_secondary_r2)){
                    secondary_index = i;
                    best_secondary_r2 = r2;
                }

                if(effect > besteffect){
                    effect_index = i;
                    besteffect = effect;
                }
            }
        }

        if(abs(bestr2)< parameters.minr && abs(best_secondary_r2) < parameters.minr && (!parameters.use_effect_threshold || besteffect < parameters.minr)){   // test whether to start new group

            //sets up new group

            if(!parameters.IsCojo() || v.minpvalue < parameters.threshold){
                auto best_stat = parameters.ancestry ? First_abs(v.stats) : Max_abs(v.stats);

                auto bestindex = parameters.ancestry ? 1 : Max_index(v.stats);

                SNPvector snpvec{v, temp};

                double zfactor = 1.0 + 1.0/parameters.CalcMaxPsi(best_stat);

                MatrixXd corr = MatrixXd::Constant(1,1, zfactor);

                snpgroups.push_back(SNPgrouping{{snpvec},VectorXd::Constant(1,best_stat), VectorXd::Constant(1,bestindex), corr, corr});
            }
        }
        else{
            // test whether to add to existing group and choose index based on best r2

            if(max(abs(best_secondary_r2), abs(bestr2)) < parameters.max_individual_r){
                if(parameters.use_effect_threshold && besteffect > abs(best_secondary_r2) && besteffect > abs(bestr2)){
                    snpgroups[effect_index].UpdateIfSignificant(v, zthreshold, temp, parameters);
                }
                else if(abs(best_secondary_r2) <= abs(bestr2)){
                    snpgroups[index].UpdateIfSignificant(v, zthreshold, temp, parameters);
                }
                else{
                    snpgroups[secondary_index].UpdateIfSignificant(v, zthreshold, temp, parameters);
                }
            }
        }
    }

    // print out selected SNPs and their statistics

    for(auto& v:snpgroups){
        forc(i, v.snps){

            if(parameters.method=="pca"){
                Printdelim(cout, '\t', v.snps.size(), i+1, v.snps[i].snpinfo.rsnum, 0, v.snps[i].snpinfo.minpvalue, 0, parameters.chromosome, v.snps[i].snpinfo.position, v.snps[i].snpinfo.effect, v.snps[i].snpinfo.baseline, v.snps[i].snpinfo.eaf, v.snps[i].snpinfo.r2);

                for(auto& w:v.snps[i].snpinfo.stats){
                    cout << '\t' << get<0>(w) << '\t' << get<1>(w);
                }
            }
            else{                
               double pvalue;

               try{
                   pvalue = 2*cdf(gaussian, -v.best_statistics[i]);
               }
               catch(...){
                   pvalue = -1.0;
               }

                Printdelim(cout, '\t', v.snps.size(), i+1, v.snps[i].snpinfo.rsnum, v.best_statistics[i], pvalue, v.best_index[i], parameters.chromosome, v.snps[i].snpinfo.position, v.snps[i].snpinfo.effect, v.snps[i].snpinfo.baseline, v.snps[i].snpinfo.eaf, v.snps[i].snpinfo.r2);

                for(auto& w:v.snps[i].snpinfo.stats){
                    cout << '\t' << get<0>(w) << '\t' << get<1>(w);
                }

                if(parameters.functional_scores){
                    cout << '\t' << v.snps[i].snpinfo.functional_score;
                }
            }

            if(parameters.IsCojo()){
//                double best_ind_statistic = 0;

//                for(auto& w: v.snps[i].snpinfo.stats){
//                    if(get<1>(w)>0 && abs(get<0>(w)/get<1>(w))>best_ind_statistic){
//                        best_ind_statistic = abs(get<0>(w)/get<1>(w));
//                    }
//                }

                double best_ind_statistic = parameters.ancestry ? First_abs(v.snps[i].snpinfo.stats) : Max_abs(v.snps[i].snpinfo.stats);

                double pvalue;

                try{
                    pvalue = 2*cdf(gaussian, -best_ind_statistic);
                }
                catch(...){
                    pvalue = -1.0;
                }

                cout << '\t' << pvalue;

                auto conditional_stats = v.CalcConditionalStats(parameters.ancestry);

                for(auto& w:conditional_stats){
                    cout << '\t' << get<0>(w)[i] << '\t' << get<1>(w)[i];
                }
            }

            cout << '\n';
        }
    }
}


// This is the command line routine to select SNPs and get their statistics

int CL_CondAnalysis(int argc, char* argv[])
{
    vector<string> inputfiles;

    string datafile, chromosome, includefile="", method="s4", refsnpfile="", genmapfile="";

    // test_threshold is only used in a cojo analysis and is the threshold for secondary hits

    double maf = 0, threshold = 1.0e-5, test_threshold = -1.0, cojo_threshold = 1.0e-3, minr = 0.02, maxdistance = 1e7, maxr = sqrt(0.6), r2impute=0.5;

    bool correctvalue = false;

    int use_multiple_r2 = -1;

    double nvalue=0, zchange = 100, indiv_change = 100;

    double max_individual_r = 10;

    bool novarqc = false;

    double corr_damp = 1.0;

    bool logresults = false;

    bool use_effect_threshold = false;

    double subset_threshold = -1.0;

    string hybrid_file = "";

    double hybrid_threshold = -1;

    double minindr = 0.005;

    vector<double> qc;

    bool ancestry = false;

    vector<double> bounds;

    string excludesnps;

    double regval = 0.0;

    bool functional_scores = false;

    try
    {
        po::options_description desc;       // this is a boost command line variable that specifies how the variables are read from the command line

        desc.add_options()("help,h", "produce help message")
                (",a", po::value<vector<string>>(&inputfiles)->multitoken(), "analysis files")
                (",c", po::value<string>(&chromosome)->required(), "chromosome")
                (",m", po::value<double>(&maf), "minimum MAF")
                (",p", po::value<double>(&threshold), "p value threshold")
                (",t", po::value<double>(&test_threshold), "threshold for testing snp")
                (",r", po::value<double>(&minr), "threshold for putting into correlated group")
                ("r2impute,", po::value<double>(&r2impute), "imputation threshold for including snp")
                ("maxr,", po::value<double>(&maxr), "threshold for maximum absolute r to consider")
                ("dist,", po::value<double>(&maxdistance), "max distance to consider being in same correlated group")
                (",d", po::value<string>(&datafile)->required(), "data file")
                ("correct,", po::value<bool>(&correctvalue)->zero_tokens(), "correct r2 based on variance of statistic")
                (",n", po::value<double>(&nvalue), "average number value")
                ("multiple,", po::value<int>(&use_multiple_r2), "use multiple r2 if set to 1 (default 0 if method S4, 1 otherwise)")
                (",z", po::value<double>(&zchange), "maximum increase in absolute z score")
                (",i", po::value<string>(&includefile), "include samples from ld file")
                ("method,", po::value<string>(&method), "method for adding SNPs to model, default is for for S4 selection, alternatives are cojo and pca")
                ("refsnpfile,", po::value<string>(&refsnpfile), "file containing the SNP positions in the reference file")
                ("genmapfile,", po::value<string>(&genmapfile), "use specified genetic map file for distance rather than position")
                ("maxindr,", po::value<double>(&max_individual_r), "maximum absolute individual r")
                ("novarqc,", po::value<bool>(&novarqc)->zero_tokens(), "don't reject SNPs based on failing QC on expected variance")
                ("corrdamp,", po::value<double>(&corr_damp), "damp the correlation by specified amount")
                ("log,", po::value<bool>(&logresults)->zero_tokens(), "correct r2 based on variance of statistic")
                ("qc,", po::value<vector<double>>(&qc)->multitoken(), "QC parameters for accepting SNPs")
                ("zindiv,", po::value<double>(&indiv_change), "maximum increase in absolute z score for specific snp added")
                ("cojothres,", po::value<double>(&cojo_threshold), "threshold for considering snp in a COJO analysis")
                ("effect,", po::value<bool>(&use_effect_threshold)->zero_tokens(), "test for change in z-score for significant SNPs")
                ("subset,", po::value<double>(&subset_threshold), "subset threshold to test for multiple correlation and z-score increase")
                ("hybrid_file,", po::value<string>(&hybrid_file), "hybrid file containing variants to prefer")
                ("hybrid_threshold,", po::value<double>(&hybrid_threshold), "Choose variants in hybrid file below threshold")
                ("ancestry,", po::value<bool>(&ancestry)->zero_tokens(), "Data is based on multiple ancestries")
                ("psibounds,", po::value<vector<double>>(&bounds)->multitoken(), "Psi bounds for different zscores")
                ("minindr,", po::value<double>(&minindr), "minimum individual correlation to put into correlated group when using effect (default 0.005)")
                ("excludesnps,", po::value<string>(&excludesnps), "amb to exlude ambiguous, indel to exlude indels and ambindel or indelamb to exclude both")
                ("reg,", po::value<double>(&regval), "regularisation value")
                ("functional,", po::value<bool>(&functional_scores)->zero_tokens(), "use functional scores");

        po::variables_map vm;

        po::store(parse_command_line(argc, argv, desc), vm);

        if (vm.count("help"))
        {
            cout << desc;

            return 1;
        }

        notify(vm);
    }
    catch (const po::error& e)
    {
        cout << e.what() << '\n';

        return -1;
    }

    if(use_multiple_r2==-1){
        if(method=="s4"){
            use_multiple_r2 = 0;
        }
        else{
            use_multiple_r2 = 1;
        }
    }

    // select SNPs and get statistics using the parameters specified from the command line

    SeveralCorrelatedAnalyses(datafile, inputfiles, {maf, threshold, test_threshold, minr, maxr, maxdistance, r2impute, chromosome, correctvalue, nvalue, bool(use_multiple_r2), zchange, includefile, method, refsnpfile, genmapfile, 10, max_individual_r, novarqc, corr_damp, logresults, qc, indiv_change, cojo_threshold, use_effect_threshold, subset_threshold, hybrid_file, hybrid_threshold, ancestry, bounds, minindr, excludesnps, regval, functional_scores});

    return 0;
}


// This is the command line routine to select SNPs and get their statistics

int CL_SelectedLS(int argc, char* argv[])
{
    vector<string> inputfiles;

    string datafile, chromosome, includefile="", method="s4", refsnpfile="", genmapfile="";

    // test_threshold is only used in a cojo analysis and is the threshold for secondary hits

    double threshold = 1.0e-5, maxr = sqrt(0.6), r2impute=0.5;

    bool correctvalue = false;

    double nvalue=0;

    string hybrid_file = "";

    vector<double> qc;

    bool ancestry = false;

    string summarystats_file, pgm_file, bgen_reffile_prefix, bgenixreffile_prefix, include_file, geneticmap_file;

    double genetic_distance = 3.0;

    string outputprefix;

    bool uncorrected;

    try
    {
        po::options_description desc;       // this is a boost command line variable that specifies how the variables are read from the command line

        desc.add_options()("help,h", "produce help message")
                (",p", po::value<double>(&threshold), "p value threshold")
                (",r", po::value<double>(&maxr), "threshold for correlation between existing variants")
                ("r2impute,", po::value<double>(&r2impute), "imputation threshold for including snp")
                ("dist,", po::value<double>(&genetic_distance), "maximum genetic distance to consider as uncorrelated")
                (",d", po::value<string>(&bgen_reffile_prefix)->required(), "data file")
                ("correct,", po::value<bool>(&correctvalue)->zero_tokens(), "correct r2 based on variance of statistic")
                (",n", po::value<double>(&nvalue), "average number value")
                (",i", po::value<string>(&include_file), "include samples from ld file")
                ("refsnpfile,", po::value<string>(&bgenixreffile_prefix), "file containing the SNP positions in the reference file")
                ("genmapfile,", po::value<string>(&geneticmap_file), "use specified genetic map file for distance rather than position")
                ("hybridfile,", po::value<string>(&hybrid_file), "hybrid file containing variants to prefer")
                ("ancestry,", po::value<bool>(&ancestry)->zero_tokens(), "Data is based on multiple ancestries")
                (",o", po::value<string>(&outputprefix)->required(), "output prefix")
                (",s", po::value<string>(&summarystats_file)->required(), "summary stats files")
                ("pgm,", po::value<string>(&pgm_file)->required(), "pgm file")
                ("uncorrected,", po::value<bool>(&uncorrected)->zero_tokens(), "Get the uncorrected results, rather than the corrected")
                ;

        po::variables_map vm;

        po::store(parse_command_line(argc, argv, desc), vm);

        if (vm.count("help"))
        {
            cout << desc;

            return 1;
        }

        notify(vm);
    }
    catch (const po::error& e)
    {
        cout << e.what() << '\n';

        return -1;
    }

    // select SNPs and get statistics using the parameters specified from the command line

    GetLowSignificancePGS(summarystats_file, pgm_file, hybrid_file, bgen_reffile_prefix, bgenixreffile_prefix, include_file, geneticmap_file, genetic_distance, outputprefix, {threshold, maxr, r2impute, correctvalue, nvalue, uncorrected, ancestry});

    return 0;
}


// Get file positions using bgi file generated by bgenix
// The purpose of this is to store the positions in these files rather than the data and
// use the data when needed, as storing the data will use a lot of memory.

map<Posclass, tuple<int, int64_t>> GetFilePositions(const string& inputfile, const string& bgi_prefix){
    zifstream input(inputfile);

    SplitString splitstr;

    string rsnum;

    map<int, set<string>> rsids;

    while(input >> splitstr){

        try{
            int chr = (splitstr[6] =="X" ? 23 : stoi(splitstr[6]));

            rsids[chr].insert(splitstr[2]);
        }
        catch(...){
            cout << "Error in line processing\n";

            cout << splitstr.linestr << '\n';

            exit(0);
        }
    }

    map<Posclass, tuple<int, int64_t>> result;

    forv(chr,1,24){
        if(![&]{
            ifstream temp_input(Paste(bgi_prefix, chr, ".bgen"));

            return bool(temp_input);
        }()){
            continue;
        }

        genfile::bgen::SqliteIndexQuery query(Paste(bgi_prefix, chr, ".bgen.bgi"));

        BgenParser bgenParser(Paste(bgi_prefix, chr, ".bgen"));

        query.include_rsids(vector<string>(rsids[chr].begin(), rsids[chr].end()));

        query.initialise();

        for(size_t i=0;i<query.number_of_variants();i++){

            auto fileposition = query.locate_variant(i);

            bgenParser.JumpToPosition(fileposition.first);

            std::string chromosome ;
            uint32_t position ;
            std::string rsid ;
            std::vector< std::string > alleles ;

            bgenParser.read_variant( &chromosome, &position, &rsid, &alleles );

            result[Posclass(ConvertChrStr(chromosome), position, alleles[0], alleles[1])] = make_tuple(chr, fileposition.first);
        }
    }

    return(result);
}


// This estimates the PRS based on the weightings given by the data structure
// normally only one file is used but this routine was orignally written for several bgen files

MatrixXd GetMultiplePrincipalComponentValuesFromIds(const vector<string>& datafile, map<Posclass, vector<double>>& data, const vector<string>& rsids, int print=0)
{
    string linestr, rsnum;

    Posclass posclass;

    Imputeclass imputeclass;

    uint32_t position;

    vector<string> alleles;

    vector<vector<double>> probs;

    int nsnps = 0;

    MatrixXd prsmat;

    bool first_run = true;

    auto UpdatePRSdosages = [&](BgenParser& bgenparse){
        tie(imputeclass.ref, imputeclass.alt) = {alleles[0],alleles[1]};

        imputeclass.position = int(position);

        auto itor = data.find(Posclass{imputeclass});

        bool swapped = false;

    // insertion deletions are not mapped if the order is different so check whether this is the case and mark for flipping

        if(itor==data.end() && (imputeclass.ref.size()!=1 || imputeclass.alt.size()!=1)){
            swap(imputeclass.ref, imputeclass.alt);

            itor = data.find(Posclass{imputeclass});

            swap(imputeclass.ref, imputeclass.alt);

            swapped = true;
        }

        if(itor != data.end())
        {
            //                bgenParser.ReadAllProbs(imputeclass);

            //            bgenparse.read_probs(&probs);

            //            if(prs.size()==0) prs = VectorXd::Zero(probs.size());

            int factor = MatchAlleles(Posclass{imputeclass}, itor->first);

            if(swapped){
                factor = -1;        // insertion deletion that needs to be flipped
            }

//            bgenparse.ReadAllImputeclassDosages(imputeclass);

            if(first_run){
                prsmat.setZero(itor->second.size()-1, imputeclass.size());

                first_run = false;
            }

            Map<VectorXd> rate(itor->second.data(), itor->second.size()-1);

            bgenparse.UpdateMultiPGS(prsmat, factor, rate, 2*itor->second.back());

//            forc(i, imputeclass){
//                double meandiff = 2*(factor==-1) + factor*imputeclass[i] - 2*itor->second.back();

//                Map<VectorXd> eigenVector(itor->second.data(), itor->second.size()-1);

//                prsmat.col(i) += meandiff*eigenVector;
////                forl(j, itor->second.size()-1){
////                    prsmat(j,i) += meandiff * itor->second[j];
////                }
//            }

            if(print){
                ++nsnps;

                if(nsnps%print==0){
                    cout << nsnps << " " << imputeclass.rsnum << '\n';
                }
            }
        }
        else
        {
            bgenparse.ignore_probs();
        }
    };

    for(auto& v: datafile)
    {
        BgenParser bgenParser(v);

        if(!rsids.empty()){
            genfile::bgen::SqliteIndexQuery query(Paste(v,".bgi"));

            query.include_rsids(rsids);

            query.initialise();

            if(print) cout << query.number_of_variants() << '\n';

            for(size_t i=0;i<query.number_of_variants();i++){

                auto fileposition = query.locate_variant(i);

                bgenParser.JumpToPosition(fileposition.first);

                bgenParser.read_variant(&imputeclass.chromosome, &position, &imputeclass.rsnum, &alleles);

                UpdatePRSdosages(bgenParser);
            }
        }
        else{
            while(bgenParser.read_variant(&imputeclass.chromosome, &position, &imputeclass.rsnum, &alleles))
            {
                UpdatePRSdosages(bgenParser);
            }
        }
    }

    return prsmat;
}


// This estimates the PRS based on the weightings given by the data structure
// normally only one file is used but this routine was orignally written for several bgen files

MatrixXd GetMultiplePrincipalComponentValuesFromSeveralBgenFiles(const vector<string>& datafile, map<Posclass, vector<double>>& data, const string& bgenix_prefix="", const vector<string>& chrx_files={}, int print=0)
{
    string linestr, rsnum;

    Posclass posclass;

    Imputeclass imputeclass;

    vector<string> alleles;

    vector<vector<double>> probs;

    vector<string> rsids=[&]{
        if(bgenix_prefix!=""){
            auto unique_ids=ConvertPositionsToIds(bgenix_prefix, data);

            return vector<string>(unique_ids.begin(), unique_ids.end());
        }
        else{
            return vector<string>();
        }
    }();

    auto value = GetMultiplePrincipalComponentValuesFromIds(datafile, data, rsids, 100);

    if(chrx_files.size()==2){
        auto value2 = GetMultiplePrincipalComponentValuesFromIds({chrx_files[0]}, data, rsids, 100);

        zifstream input(chrx_files[1]);

        vector<pair<int,int>> conv_indices;

        for(pair<int,int> temppair;input >> temppair;){
            conv_indices.emplace_back(temppair);
        }

        if(value.size()==0){
            int vecsize = max_element(conv_indices.begin(), conv_indices.end())->first + 1;
            value = VectorXd::Zero(vecsize);
        }

        for(auto v: conv_indices){
            if(unsigned(v.first)>=unsigned(value.size()) || unsigned(v.second) >= unsigned(value2.size())){
                cout << "index out of bounds\n";
                exit(1);
            }

            value.col(v.first) += value2.col(v.second);
        }
    }

    return value;
}


// psi function used for continuous shrinkage estimate
// adapted from python code in PRScs

double psi(double x, double alpha, double lam){
    double f = -alpha*(cosh(x)-1)-lam*(exp(x)-x-1);
    return f;
}

// dpsi function used for continuous shrinkage estimate
// adapted from python code in PRScs

double dpsi(double x, double alpha, double lam){
    double f = -alpha*sinh(x)-lam*(exp(x)-1);
    return f;
}

// g function used for continuous shrinkage estimate
// adapted from python code in PRScs

double g(double x, double sd, double td, double f1, double f2){

    double f=0;

    if (x >= -sd && x <= td){
        f = 1;
    }
    else if(x > td){
        f = f1;
    }
    else if(x < -sd){
        f = f2;
    }

    return f;
}


// gig random number generator used for continuous shrinkage estimate
// adapted from python code in PRScs
// gen is the Mersenne twister random number generator

double gigrnd(double p, double a, double b, mt19937& gen){

    // setup -- sample from the two-parameter version gig(lam,omega)

    double lam = p;

    double omega = sqrt(a*b);

    bool swap;

    if(lam < 0){
        lam = -lam;
        swap = true;
    }
    else{
        swap = false;
    }

    double alpha = sqrt(omega*omega+lam*lam)-lam;

    // find t

    double x = -psi(1, alpha, lam);

    double t = 0;

    if (x >= 0.5 && x <= 2){
        t = 1.0;
    }
    else if(x > 2){
        if(alpha==0 && lam==0){
            t = 1.0;
        }
        else{
            t = sqrt(2.0/(alpha+lam));
        }
    }
    else if(x < 0.5){
        if(alpha==0 && lam==0){
            t = 1.0;
        }
        else{
            t = log(4.0/(alpha+2.0*lam));
        }
    }

    // find s

    x = -psi(-1.0, alpha, lam);

    double s = 0;

    if (x >= 0.5 && x <= 2){
        s = 1;
    }
    else if(x > 2){
        if(alpha==0 && lam==0){
            s = 1.0;
        }
        else{
            s = sqrt(4.0/(alpha*cosh(1.0)+lam));
        }
    }
    else if(x < 0.5){
        if(alpha==0 && lam==0){
            s = 1.0;
        }
        else if(alpha==0){
            s = 1.0/lam;
        }
        else if(lam==0){
            s = log(1.0+1.0/alpha+sqrt(1.0/pow(alpha,2)+2.0/alpha));
        }
        else{
            s = min(1.0/lam, log(1+1.0/alpha+sqrt(1.0/pow(alpha,2)+2.0/alpha)));
        }
    }

    // find auxiliary parameters

    double eta = -psi(t, alpha, lam);
    double zeta = -dpsi(t, alpha, lam);
    double theta = -psi(-s, alpha, lam);
    double xi = dpsi(-s, alpha, lam);

    p = 1/xi;
    double r = 1/zeta;

    double td = t-r*eta;
    double sd = s-p*theta;
    double q = td+sd;

    uniform_real_distribution<> randnum(0.0,1.0);

    double rnd = 0;

    // random variate generation
    while(1){
        double U = randnum(gen);
        double V = randnum(gen);
        double W = randnum(gen);
        if(U < q/(p+q+r)){
            rnd = -sd+q*V;
        }
        else if(U < (q+r)/(p+q+r)){
            rnd = td-r*log(V);
        }
        else{
            rnd = -sd+p*log(V);
        }

        double f1 = exp(-eta-zeta*(rnd-t));
        double f2 = exp(-theta+xi*(rnd+s));
        if(W*g(rnd, sd, td, f1, f2) <= exp(psi(rnd, alpha, lam))){
            break;
        }
    }

    // transform back to the three-parameter version gig(p,a,b)

    rnd = exp(rnd)*(lam/omega+sqrt(1+pow(lam,2)/pow(omega,2)));

    if(swap){
        rnd = 1.0/rnd;
    }

    rnd = rnd/sqrt(a/b);

    return rnd;
}


// stucture used to pass parameters to shrinkage model

struct MCMC_parameters{
    int n_iter;         // no. of iterations
    int n_burnin;       // n0. of burnin iterations
    int thin;           // no. of runs between each recorded value
    double a;           // a,b and phi are parameters passed to the continuous shrinkage model
    double b;
    double phi;
    double inflation_factor;        // inflation factor in case of population stratification (rarely used)
    double max_phi;                 // maximum phi value (rarely used)
    double region_significance;     // the region parameters could be used if other SNPs in the region are very significant,
    double region_factor;           // though the experiments so far didn't improve the model
    int region_distance;
    double sparse;
    int eaf_variance = 0;
    double zdiff;
    int threshold_to_increase_runs = -1;
    vector<double> approximation_threshold;
    vector<double> bounds;
    string method;
    vector<double> pvec;    // used for correlated phenotypes
    vector<double> wvec;    // used for correlated phenotypes
    double k = 1.0;         // extra option for generalised gamma (default will be 1)
    double min_psi = 0.0;   // check if bounding above zero reduces variance of MCMC sampling
    double regval = 0.0;    // regularisation factor for horseshoe model
    bool use_ss_r2;
    double eaf_power = 0.0;
    double r2_power = 1.0;
    bool functional_scores = false;

    double CalcMaxPsi(double zscore) const
    {
        if(bounds.empty()) return 1.0e50;

        double result = bounds[0];

        for(int i=2;i<int(bounds.size());i+=2){
            if(abs(zscore) < bounds[i-1]){
                result = bounds[i];
                break;
            }
        }

        return result;
    };
};


// This shrinks the estimates contained in the vector v based on the correlation matrix and the specified continuous shrinkage parameters
// Most of this code was adapted from the PRScs python code

VectorXd CalculateContinuousShrinkageEstimate(const vector<SNPdetails>& v, const MatrixXd& corr, const MCMC_parameters& mcmc_parameters, std::mt19937& gen){

    VectorXd b_marginal(v.size());

    VectorXd se_marginal(v.size());

    VectorXd phicorrection(v.size());

    boost::math::normal gaussian;

    // Put statistics in vectors b_marginal and se_marginal

    forc(i, b_marginal){
        b_marginal[i] = get<0>(v[i].stats[0]);

        se_marginal[i] = get<1>(v[i].stats[0])*mcmc_parameters.inflation_factor;
   }

    VectorXd zscores = b_marginal.cwiseQuotient(se_marginal);

    if(mcmc_parameters.method == "stepwise"){

        LLT<MatrixXd> zllt(corr);

        VectorXd beta_est = zllt.solve(zscores);

        return beta_est.cwiseProduct(se_marginal);
    }

    VectorXd maxpsi = zscores.unaryExpr([&](double x){return mcmc_parameters.CalcMaxPsi(x);});

    int n_pst = (mcmc_parameters.n_iter-mcmc_parameters.n_burnin)/mcmc_parameters.thin;


            // initialisation

    VectorXd beta(zscores.size());

    VectorXd beta_mean(zscores.size());

    VectorXd psi(zscores.size());

    VectorXd delta(zscores.size());

    VectorXd beta_est(zscores.size());


    beta = beta_mean = delta = beta_est = VectorXd::Zero(zscores.size());

    psi = VectorXd::Ones(zscores.size());

    VectorXd phivec = mcmc_parameters.phi*se_marginal.cwiseInverse().cwiseAbs2();

    forc(i, phivec){
        if(mcmc_parameters.use_ss_r2){
            phivec[i] *= pow(1.2*v[i].ss_r2, mcmc_parameters.r2_power) * v[i].functional_score;
        }
        else{
            phivec[i] *= pow(1.2*(v[i].ref_r2 > 0 ? v[i].ref_r2 : 0.1), mcmc_parameters.r2_power) *v[i].functional_score;     // for monomorphic reference set default of 0.1
        }
    }

    if(mcmc_parameters.eaf_variance==1){
        forc(i,phivec){
            if(v[i].ref_r2 > 0){                                    // check if monomorphic reference
                phivec[i] /= pow(v[i].eaf*(1-v[i].eaf)/0.0475, mcmc_parameters.eaf_power);

                phivec[i] *= v[i].functional_score;
            }
        }
    }

            // MCMC

    std::normal_distribution<> rnorm{0,1};

    int ntests = mcmc_parameters.n_iter - mcmc_parameters.n_burnin;

    if(zscores.size()<mcmc_parameters.threshold_to_increase_runs){
        ntests *= double(mcmc_parameters.threshold_to_increase_runs)/zscores.size();
    }

    vector<array<int,10>> nband(zscores.size());

    forv(itr, 1,mcmc_parameters.n_burnin+ntests+1){



        vector<int> indices, indices2, indices3;

        forc(i, psi){
            if(psi[i] > mcmc_parameters.approximation_threshold[0]){
                indices.push_back(i);
            }
            else if(psi[i]>mcmc_parameters.approximation_threshold[1]){
                indices3.push_back(i);
            }
            else{
                indices2.push_back(i);
            }
        }

//        cout << indices.size() << ' ' << indices3.size() << ' ' << indices2.size() << '\n';

//        if(indices3.size()==1){
//            indices.clear();
//            indices3.clear();

//            forl(i, zscores.size()-1){
//                indices.push_back(i);
//            }

//            indices3.push_back(zscores.size()-1);
//        }

        VectorXd psi2 = psi(indices);

        MatrixXd dinvt2 = MatrixXd(psi2.asDiagonal().inverse()) + corr(indices, indices);

        if(mcmc_parameters.regval > 1e-10){
            dinvt2.diagonal().array() += mcmc_parameters.regval;
        }
        else if(mcmc_parameters.method=="reg"){
            dinvt2 += maxpsi.asDiagonal().inverse();
        }

        MatrixXd transform2 = LLT<MatrixXd>(dinvt2).matrixL();

        VectorXd beta_tmp2 = transform2.triangularView<Lower>().solve(zscores(indices));

        VectorXd beta_tmp_mean2 = beta_tmp2;

        VectorXd normal_dist(zscores.size());

        forc(i, normal_dist){
            normal_dist[i] = rnorm(gen);
        }

        beta_tmp2 += normal_dist(indices);

        VectorXd beta2;

//        VectorXd beta2 = transform2.triangularView<Lower>().transpose().solve(beta_tmp2);

        MatrixXd ch2;

        VectorXd dch, beta_rem, beta_rem_tmp, beta_rem_result, beta_rem_mean, beta2_mean, beta3_mean;

        if(indices3.size()>0){
            ch2 = transform2.triangularView<Lower>().solve(corr(indices, indices3));

//            cout << ch2 << "\n\n";

//            cout << ch2.cwiseAbs2().colwise().sum() << "\n\n";

            dch = (psi(indices3).array().inverse() + 1.0 - ch2.cwiseAbs2().colwise().sum().transpose().array()).sqrt();

            beta_rem = (zscores(indices3) - ch2.transpose()*beta_tmp_mean2).cwiseQuotient(dch);

            beta_rem_tmp = beta_rem + normal_dist(indices3);

            beta_rem_result = beta_rem_tmp.cwiseQuotient(dch);

            beta2 = transform2.triangularView<Lower>().transpose().solve(beta_tmp2 - ch2*beta_rem_result);
        }
        else{
            beta2 = transform2.triangularView<Lower>().transpose().solve(beta_tmp2);
        }


        VectorXd psi3 = (psi(indices2).array().inverse() + 1).inverse();

        VectorXd beta_tmp3 = zscores(indices2).cwiseProduct(psi3.cwiseSqrt());

        VectorXd beta_tmp_mean3 = beta_tmp3;

        beta_tmp3 += normal_dist(indices2);

        VectorXd beta3 = beta_tmp3.cwiseProduct(psi3.cwiseSqrt());

        VectorXd beta_all(zscores.size());

        beta_all(indices) = beta2;

        beta_all(indices2) = beta3;

        beta_all(indices3) = beta_rem_result;

        if(!beta_all.allFinite()){

            cout << "Failed on run " << itr << '\n';

            forc(i, beta_all){
//                if(!isfinite(beta_all[i])){
                    cout << i+1 << '\t' << beta_all[i] << '\t' << v[i].rsnum << '\t' << get<0>(v[i].stats[0]) << '\t' << get<1>(v[i].stats[0]) << '\n';
//                }
            }
            exit(1);
        }



//        forc(i, beta_tmp2){
//            beta_tmp2[i] += rnorm(gen);
//        }


//        MatrixXd dinvt = MatrixXd(psi.asDiagonal().inverse()) + corr;

//        MatrixXd transform = LLT<MatrixXd>(dinvt).matrixL();

//        VectorXd beta_tmp = transform.triangularView<Lower>().solve(zscores);

//        VectorXd beta_tmp_mean = beta_tmp;

//        beta_tmp += normal_dist;

//////        forc(i, beta_tmp){
//////            beta_tmp[i] += rnorm(gen);
//////        }

//        beta = transform.triangularView<Lower>().transpose().solve(beta_tmp);

//        cout << beta_all - beta << "\n\n";


        if(abs(mcmc_parameters.k - 1.0) < 1.0e-5){
            forc(i, delta){
                delta[i] = gamma_distribution<>(mcmc_parameters.a+mcmc_parameters.b, 1.0/(psi[i]+phivec[i]))(gen);
            }
        }
        else{
            forc(i, delta){
                double random_gamma = gamma_distribution<>(mcmc_parameters.a+mcmc_parameters.b, 1.0)(gen);

                delta[i] = pow(random_gamma, 1.0/mcmc_parameters.k)/(psi[i]+phivec[i]);
            }
        }

        forc(i, psi){
            psi[i] = gigrnd(mcmc_parameters.a-0.5, 2*delta[i], Sqr(beta_all[i]), gen);

            if(psi[i] > maxpsi[i] && mcmc_parameters.method!="reg" && mcmc_parameters.regval<1e-10) psi[i] = maxpsi[i];

            if(psi[i] < mcmc_parameters.min_psi) psi[i] = mcmc_parameters.min_psi;
        }

        // posterior
        if(itr>mcmc_parameters.n_burnin && itr % mcmc_parameters.thin == 0){

//            beta_mean = transform.triangularView<Lower>().transpose().solve(beta_tmp_mean);

            beta_rem_mean = beta_rem.cwiseQuotient(dch);

            if(indices3.size()>0){
                beta2_mean = transform2.triangularView<Lower>().transpose().solve(beta_tmp_mean2 - ch2*beta_rem_mean);
            }
            else{
                beta2_mean = transform2.triangularView<Lower>().transpose().solve(beta_tmp_mean2);
            }

            beta3_mean = beta_tmp_mean3.cwiseProduct(psi3.cwiseSqrt());

            VectorXd beta_mean_all(zscores.size());

            beta_mean_all(indices) = beta2_mean;

            beta_mean_all(indices2) = beta3_mean;

            beta_mean_all(indices3) = beta_rem_mean;

            if(!beta_mean_all.allFinite()){
                forc(i, beta_mean_all){
                    cout << "Mean calculation\n";
                    if(!isfinite(beta_mean_all[i])){
                        cout << i+1 << '\t' << beta_mean_all[i] << '\t' << v[i].rsnum << '\n';
                    }
                }
                exit(1);
            }

//            cout << beta_mean - beta_mean_all << "\n\n";

            bool extreme_values = false;

            if(mcmc_parameters.zdiff > 0)
            forc(i, beta_mean_all){
                if(abs(beta_mean_all[i])-abs(zscores[i])>mcmc_parameters.zdiff){
                    extreme_values = true;
                }
            }

            if(extreme_values){
                itr--;
//                cout << "Extreme value found\n";
            }
            else{
                beta_est = beta_est + beta_mean_all/n_pst;
            }
        }
    }

//    forc(i, nband){
//        forc(j,nband[i]){
//            cout << i << '\t' << j << '\t' << nband[i][j] << '\n';
//        }
//    }

    return beta_est.cwiseProduct(se_marginal);

//    VectorXd result = beta_est.cwiseProduct(se_marginal);

//    if(mcmc_parameters.eaf_variance==1){
//        forc(i,result){
//            result[i] *= sqrt(v[i].eaf*(1-v[i].eaf)/0.09);
//        }
//    }

//    return result;
}


VectorXd SpikeAndSlab(const vector<SNPdetails>& v, const MatrixXd& corr, const MCMC_parameters& mcmc_parameters, std::mt19937& gen){

    VectorXd b_marginal(v.size());

    VectorXd se_marginal(v.size());

    boost::math::normal gaussian;

    // Put statistics in vectors b_marginal and se_marginal
    // perhaps remove region factors later to simplify code

    forc(i, b_marginal){
        b_marginal[i] = get<0>(v[i].stats[0]);

        se_marginal[i] = get<1>(v[i].stats[0]);
    }

    VectorXd zscores = b_marginal.cwiseQuotient(se_marginal);

            // initialization

    VectorXd beta(zscores.size());

    VectorXd beta_est(zscores.size());

    beta = beta_est = VectorXd::Zero(zscores.size());

            // MCMC

    std::normal_distribution<> rnorm{0,1};

    std::uniform_real_distribution<> runiform{0.0,1.0};

    int ntests = mcmc_parameters.n_iter - mcmc_parameters.n_burnin;

    forv(itr, 1,mcmc_parameters.n_burnin+ntests+1){

        forc(i, zscores){
            double res = zscores[i] + beta[i] - beta.dot(corr.col(i));

            double c1 = Sqr(mcmc_parameters.phi);

            double c2 = 1.0 / (1.0 + 1.0/c1);

            double c3 = c2*res;

            double c4 = sqrt(c2);

            double post_p_j = 1.0/(1.0 + sqrt(1+c1)*exp(-Sqr(c3/c4)/2)*(1-mcmc_parameters.b)/mcmc_parameters.b);

            if(runiform(gen) < post_p_j){
                beta[i] = c3 + rnorm(gen) * c4;
            }
            else{
                beta[i] = 0.0;
            }

            if(itr>mcmc_parameters.n_burnin){
                beta_est[i] += (c3*post_p_j) / ntests;
            }
        }
    }

    return beta_est.cwiseProduct(se_marginal);
}


MatrixXd MultiPhenoShrink(const vector<SNPdetails>& v, const MatrixXd& corr, const MCMC_parameters& mcmc_parameters, std::mt19937& gen){

    MatrixXd b_marginal(3, v.size());

    MatrixXd se_marginal(3, v.size());

    boost::math::normal gaussian;

    // Put statistics in vectors b_marginal and se_marginal
    // perhaps remove region factors later to simplify code

    forl(i, v.size()){
        forl(j, 3){
            b_marginal(j,i) = get<0>(v[i].stats[j]);

            se_marginal(j,i) = get<1>(v[i].stats[j]);
        }
    }

    MatrixXd zscores = b_marginal.cwiseQuotient(se_marginal);

            // initialization

    MatrixXd beta;

    MatrixXd beta_est;

    beta = beta_est = MatrixXd::Zero(3, v.size());

            // MCMC

    std::normal_distribution<> rnorm{0,1};

    std::uniform_real_distribution<> runiform{0.0,1.0};

    int ntests = mcmc_parameters.n_iter - mcmc_parameters.n_burnin;

    VectorXd res(3), c1(3);

    forc(i,c1){
        c1[i] = Sqr(mcmc_parameters.wvec[i]);
    }

    VectorXd c2 = 1.0 / (1.0 + 1.0/c1.array());

    VectorXd c4 = c2.cwiseSqrt();

    forv(itr, 1,mcmc_parameters.n_burnin+ntests+1){

        forl(i, v.size()){

            res[0] = zscores(0,i) + beta(0,i) - beta.row(0).dot(corr.col(i));

            res[1] = zscores(1,i) + beta(1,i) - beta.row(1).dot(corr.col(i));

            res[2] = zscores(2,i) + beta(2,i) - beta.row(2).dot(corr.col(i));

            VectorXd c3 = c2.array()*res.array();

            VectorXd ratio = c3.array()/c4.array();

            VectorXd exp_values = ratio.cwiseAbs2().array()/2.0;

            double max_value = max(exp_values[0], exp_values[1]+exp_values[2]);

            VectorXd lik(5);

// to avoid numerical overflow take exp of values with a maximum of 0

            lik[0] = exp(-max_value);
            lik[1] = exp(exp_values[0]-max_value)/sqrt(1+c1[0]);
            lik[2] = exp(exp_values[1]-max_value)/sqrt(1+c1[1]);
            lik[3] = exp(exp_values[2]-max_value)/sqrt(1+c1[2]);
            lik[4] = exp(exp_values[1]+exp_values[2]-max_value)/sqrt(1+c1[1])/sqrt(1+c1[2]);

//            exp_values.array()

//            VectorXd lik = exp(exp_values.array())/sqrt(1 + c1.array());

            VectorXd wt{{1.0/Sqr(se_marginal(1,i)), 1.0/Sqr(se_marginal(2,i))}};

            wt = wt/wt.sum();

            wt = wt.cwiseSqrt();

            VectorXd plik(5);

            plik[0] = (1.0 - mcmc_parameters.pvec[0] - mcmc_parameters.pvec[1] - mcmc_parameters.pvec[2] - mcmc_parameters.pvec[3])*lik[0];

            plik[1] = mcmc_parameters.pvec[0] * lik[1];

            plik[2] = mcmc_parameters.pvec[1] * lik[2];

            plik[3] = mcmc_parameters.pvec[2] * lik[3];

            plik[4] = mcmc_parameters.pvec[3] * lik[4];

            VectorXd prob = plik/plik.sum();

            if(isnan(prob[0])){
                cout << "Error in calculating probability\n";
            }


//            double post_p_j = 1.0/(1.0 + sqrt(1+c1)*exp(-Sqr(c3/c4)/2)*(1-mcmc_parameters.b)/mcmc_parameters.b);

            double rv = runiform(gen);

            if(rv < prob[0]){
                beta.col(i) = VectorXd::Zero(3);
            }
            else if(rv < prob[0] + prob[1]){
                beta(0,i) = c3[0] + rnorm(gen) *c4[0];
                beta(1,i) = wt[0] *beta(0,i);
                beta(2,i) = wt[1] * beta(0,i);
            }
            else if(rv < prob[0] + prob[1] + prob[2]){
                beta(1,i) = c3[1] + rnorm(gen) * c4[1];
                beta(0,i) = wt[0] * beta(1,i);
            }
            else if(rv < prob[0]+prob[1]+prob[2]+prob[3]){
                beta(2,i) = c3[2] * rnorm(gen) * c4[2];
                beta(0,i) = wt[1] * beta(2,i);
            }
            else{
                beta(1,i) = c3[1] + rnorm(gen) * c4[1];
                beta(2,i) = c3[2] * rnorm(gen) * c4[2];
                beta(0,i) = wt[0] * beta(1,i) + wt[1]*beta(2,i);
            }

            if(itr>mcmc_parameters.n_burnin){                
                beta_est(0,i) += (prob[1]*c3[0] + (prob[2]+prob[4])*c3[1]*wt[0]+(prob[3]+prob[4])*c3[2]*wt[1])/ntests;
                beta_est(1,i) += (prob[1]*c3[0]*wt[0] + (prob[2]+prob[4])*c3[1])/ntests;
                beta_est(2,i) += (prob[1]*c3[0]*wt[1] + (prob[3]+prob[4])*c3[2])/ntests;

                if(isnan(beta_est(0,i))){
                    cout << "Error in calculating probability\n";
                }
            }
        }
    }

    return beta_est.cwiseProduct(se_marginal);
}


// used to calculate stepwise estimate

tuple<VectorXd, VectorXd> CalculateConditionalStats(const VectorXd& b, const VectorXd& se, const MatrixXd& corr){
    MatrixXd semat = se.asDiagonal().inverse();

    MatrixXd bmat = semat * corr * semat;

    MatrixXd variance = bmat.inverse();

    VectorXd beta = variance*semat*semat*b;

    VectorXd se_result = variance.diagonal().cwiseSqrt();

    return make_tuple(beta, se_result);
}


// This calculates a stepwise estimate from summary statistics
// it adds values until the threshold is reached and then converts to conditional estimates based on the correlation matrix

VectorXd CalculateStepwiseEstimate(const vector<SNPdetails>& v, const MatrixXd& corr, double threshold){

    VectorXd b_marginal(v.size());

    VectorXd se_marginal(v.size());

    forc(i, b_marginal){
        b_marginal[i] = get<0>(v[i].stats[0]);

        se_marginal[i] = get<1>(v[i].stats[0]);
    }

    VectorXd b_temp;

    forc(i, v){

        auto [b, se] = CalculateConditionalStats(b_marginal.head(i+1), se_marginal.head(i+1), corr.topLeftCorner(i+1, i+1));

        if(abs(b[i]/se[i]) < threshold) break;

        b_temp = b;
    }

    VectorXd b_result = VectorXd::Zero(v.size());

    forc(i, b_temp){
        b_result[i] = b_temp[i];
    }

    return b_result;
}


// This structure adapts the correlation matrix based on passed parameters

struct Adjust_parameters{
    bool adjust;                // if true then correlations are adjusted
    double base_num;            // base num at which correlations are adjusted
    string weights_file;        // include file for LD reference individuals
    string prs_file;            // file to output prs weights
    string chrX_weights_file;   // include file for chromosome X reference individuals (can be different number of individuals if haplotypes are taken)
    string genmapprefix;
    double uncorrdist;
    double corr_damp;
    int seed = -1;
};


// This adjusts the correlation matrix based on some data being less well correlated

void AdjustMatrix(const Adjust_parameters& adjust_parameters, const vector<SNPdetails>& v, MatrixXd& corr){
    if(!adjust_parameters.adjust) return;

    VectorXd b_marginal(v.size());

    VectorXd se_marginal(v.size());

    forc(i, b_marginal){
        b_marginal[i] = get<0>(v[i].stats[0]);

        se_marginal[i] = get<1>(v[i].stats[0]);
    }

    forc(i, b_marginal){
        forl(j, i){
            if(v[i].ref_r2 > 0 && v[j].ref_r2>0){                       // only adjust correlation for entries when reference is not monomorphic
                double n1 = min(1.0, 1.0/(v[i].ref_r2*v[i].eaf*(1-v[i].eaf)*Sqr(get<1>(v[i].stats[0]))*adjust_parameters.base_num));

                double n2 = min(1.0, 1.0/(v[j].ref_r2*v[j].eaf*(1-v[j].eaf)*Sqr(get<1>(v[j].stats[0]))*adjust_parameters.base_num));

                double adjustment = sqrt(min(n1, n2)/max(n1, n2));

                corr(j,i) = corr(i,j) = adjust_parameters.corr_damp*adjustment*corr(j,i);
            }
        }
    }
}


// This shrinks the estimates based on the mcmc_parameters
// if mcmc_parameters is an MCMC_parameters structure then continuous shrinkage is used
// otherwise a stepwise estimate is used
// cv_num is the cross validation index, that in most cases won't be used if using an external test set

vector<double> GeneralPRSscores(const string& inputfile, const string& reference_file_prefix, const MCMC_parameters& mcmc_parameters, const Adjust_parameters& adjust_parameters, const string& analysis_group, const string& outputfile="", const vector<string>& testfile={}, const string& bgenix_prefix="", const vector<string>& chrx_files={}, const string& testphenotypefile=""){

// Read File and separate into correlated groups

    vector<vector<SNPdetails>> group_results;

    SplitString splitstr;

    zifstream input(inputfile);

    while(input>>splitstr){
        if(splitstr[1]=="1"){
            group_results.push_back(vector<SNPdetails>());
        }

        SNPdetails snpdetails{splitstr[2], splitstr[6], stoi(splitstr[7]), 0.0, splitstr[8], splitstr[9], stod(splitstr[10]), stod(splitstr[11]), -1, 1.0, 1.0, 1.0, {{stod(splitstr[12]), stod(splitstr[13])}}};

        if(splitstr.size()>14 + mcmc_parameters.functional_scores){
            for(int i=14;i<splitstr.size()-1-mcmc_parameters.functional_scores;i+=2){
                snpdetails.stats.push_back({stod(splitstr[i]), stod(splitstr[i+1])});
            }
        }

        if(mcmc_parameters.functional_scores){
            snpdetails.functional_score = stod(splitstr[splitstr.size()-1]);
        }

        if(snpdetails.functional_score>0 || !mcmc_parameters.functional_scores){
            group_results.back().push_back(snpdetails);
        }
    }

// Get list of chromosomes and positions to extract data

    auto snpposmap = GetFilePositions(inputfile, reference_file_prefix);

// Get estimates

    Imputeclass imputeclass, imputeclassfiltered;

    int num = 0;

    VectorXd prs;

    zofstream output;

    if(outputfile!=""){
        output.open(outputfile);
    }

    map<Posclass, vector<double>> prsmap;

    struct FlipCheck{
        string ref;
        string alt;
    };

    int nflips = 0;

    string weights_file = adjust_parameters.weights_file;

    vector<double> autosomal_residuals = (weights_file==""?vector<double>():FileRead<vector<double>>(weights_file));

    vector<int> autosomal_include;

    if(weights_file!=""){
        for(auto& v:autosomal_residuals){
            autosomal_include.push_back(v!=0);
        }
    }

    string chrX_weights_file = adjust_parameters.chrX_weights_file;

    // if chromosome X individuals not specified then assume autosomal individuals are to be used

    vector<double> chrX_residuals;

    vector<int> chrX_include;

    if(chrX_weights_file==""){
        chrX_residuals = autosomal_residuals;
        chrX_include = autosomal_include;
    }
    else{
        chrX_residuals = FileRead<vector<double>>(chrX_weights_file);

        for(auto& v:chrX_residuals){
            chrX_include.push_back(v!=0);
        }
    }

    map<int, map<double, double>> genmap;

    std::mt19937 gen;

    if(adjust_parameters.seed<0){
        std::random_device r;
        gen.seed(r());
    }
    else{
        gen.seed(adjust_parameters.seed);
    }

    // loop though the correlated groups

    for(auto& v:group_results){
        auto p = GetPosIterator(snpposmap, v[0]);

        if(p==snpposmap.end()){
            continue;
        }

        int chromosome = get<0>(p->second);

        cout << chromosome << ' ' << get<1>(p->second) << '\n';

        if(!genmap.count(chromosome)){
            genmap[chromosome] = ReadGeneticMap(adjust_parameters.genmapprefix==""?"":Paste(adjust_parameters.genmapprefix,chromosome,".txt"));
        }

//        auto genmap = ReadGeneticMap(adjust_parameters.genmapprefix==""?"":Paste(adjust_parameters.genmapprefix,chromosome,".txt"));

        vector<double> residuals = (chromosome==23?chrX_residuals:autosomal_residuals);

        vector<int> include = (chromosome==23?chrX_include:autosomal_include);

        BgenParser bgenParser(Paste(reference_file_prefix, chromosome, ".bgen"));

        int nsamples = bgenParser.number_of_samples();

        vector<double> residuals_filtered;

        bool is_weighted = false;

        if(!residuals.empty()){
            for(auto& v:residuals){

                if(v!=0){
                    residuals_filtered.push_back(v);

                    if(abs(v-1.0)> 1e-07){
                        is_weighted = true;
                    }
                }
            }
        }

//        MatrixXd mat(nsamples, v.size());

        MatrixXd matnorm(residuals.empty()?nsamples:residuals_filtered.size(), v.size());

        forc(i, v){
            auto p2 = GetPosIterator(snpposmap, v[i]);

            if(p2==snpposmap.end()){
                cerr << "Could not find snp " << v[i].rsnum << "\n";
                exit(1);
            }

//            bgenParser.JumpToPosition(get<1>(p2->second));

//            bgenParser.ReadAllImputeclass(imputeclass);

//            int flip = MatchAlleles(imputeclass, FlipCheck{v[i].baseline, v[i].effect});    // are the alleles in the same order as the test statistics

//            if(flip!=1){
//                cout << v[i].rsnum << ' ' << ++nflips << '\n';
//            }

            if(weights_file!=""){
                bgenParser.JumpToPosition(get<1>(p2->second));

                bgenParser.ReadImputeclass(imputeclassfiltered, include);

                int flip = MatchAlleles(imputeclassfiltered, FlipCheck{v[i].baseline, v[i].effect});

                if(v[0].stats.size()>1 && mcmc_parameters.method!="multipheno"){
                    matnorm.col(i) = flip * CreateAncestryNormalised(imputeclassfiltered, residuals_filtered, v[i]);

                    v[i].ref_r2 = imputeclassfiltered.r2;
                }
                else{
                    Eigen::Map<VectorXd> tempvec(imputeclassfiltered.genotypes.data(), imputeclassfiltered.size());  // doesn't matter if this is changed

                    if(is_weighted){
                        Eigen::Map<VectorXd> res_filtered(residuals_filtered.data(), residuals_filtered.size());

                        double weighted_mean = tempvec.dot(res_filtered)/tempvec.sum();

                        tempvec = (tempvec.array()-weighted_mean)*res_filtered.array();
                    }
                    else{
                        tempvec.array() -= tempvec.mean();
                    }

                    double norm = tempvec.norm();


//                    tempvec = flip * tempvec/(norm==0?1.0:norm);

                    matnorm.col(i) = flip*tempvec.transpose()/(norm==0?1.0:norm);

                    v[i].ref_r2 = imputeclassfiltered.r2;

//                    bgenParser.JumpToPosition(get<1>(p2->second));

//                    bgenParser.ReadAllImputeclass(imputeclass);

//                    double mean = 0;
//                    double total_weight = 0;
//                    forc(j, imputeclass){
//                        mean += imputeclass[j]*residuals[j];
//                        total_weight += residuals[j];
//                    }
//                    mean /= total_weight;

//                    if(isnan(mean)){
//                        cout << total_weight << '\n';

//                        exit(1);
//                    }

//                    bgenParser.JumpToPosition(get<1>(p2->second));

//                    bgenParser.ReadAllImputeclass(imputeclass);

//                    int flip = MatchAlleles(imputeclass, FlipCheck{v[i].baseline, v[i].effect});


//                    forl(j, nsamples){
//                        mat(j, i)  = flip*residuals[j]*(imputeclass[j] - mean);
//                    }

//                    if(mat.col(i).norm()>0.0){                                          // change this to account for 0 frequencies in multi-ethnic analyses
//                        mat.col(i) = mat.col(i)/mat.col(i).norm();
//                    }

//                    int count = 0;

//                    forc(j, include){
//                        if(include[j]){
//                            if(abs(tempvec[count]-matnorm(j,i))>1e-16){
//                                cout << j << ' ' << tempvec[count] - matnorm(j,i) << '\n';
//                            }
//                            count++;
//                        }
//                        else{
//                            if(matnorm(j,i)!=0){
//                                cout << j << ' ' << matnorm(j,i) << '\n';
//                            }
//                        }
//                    }


//                    VectorXd tempvec2(nsamples);

//                    forl(j, nsamples){

//                        tempvec2[j]  = flip*residuals[j]*(imputeclass[j] - mean);
//                    }

//                    tempvec2 /= tempvec2.norm();

//                    int k = 0;

//                    VectorXd tempvec3(tempvec.size());

//                    for(auto& v: tempvec2){
//                        if(v!= 0){
//                            tempvec3[k] = v;
//                            k++;
//                        }
//                    }

//                    VectorXd tempvec4 = tempvec.transpose()/norm;

////                    forc(j, tempvec3){
////                        if(abs(tempvec3[j] - tempvec4[j]) > 1e-20){
////                            cout << j << ' ' << tempvec3[j]-tempvec4[j] << '\n';
////                        }
////                    }

//                    matnorm.col(i) = tempvec3;


//                    bgenParser.JumpToPosition(get<1>(p2->second));

//                    bgenParser.ReadImputeclass(imputeclassfiltered, include);
                }
            }
            else{
                bgenParser.JumpToPosition(get<1>(p2->second));

                bgenParser.ReadAllImputeclass(imputeclass);

                int flip = MatchAlleles(imputeclass, FlipCheck{v[i].baseline, v[i].effect});    // are the alleles in the same order as the test statistics


                forl(j, nsamples){
                    matnorm(j,i) = flip*(imputeclass[j] - 2*imputeclass.eaf);
                }

                matnorm.col(i) = matnorm.col(i)/matnorm.col(i).norm();

                v[i].ref_r2 = imputeclass.r2;
            }

            v[i].ss_r2 = 1.0/(v[i].eaf*(1-v[i].eaf)*Sqr(get<1>(v[i].stats[0]))*adjust_parameters.base_num);

            if(v[i].ss_r2>1) v[i].ss_r2 = 1;
        }

        // as matnorm is normalised to have mean 0 and the sum of squares to be 1 can calculate the correlation matrix

//        MatrixXd corr = matnorm.transpose()*matnorm;

//        MatrixXd corr = mat.transpose()*mat;

        MatrixXd corr = matnorm.transpose()*matnorm;

//        if((corr.reshaped()-corr2.reshaped()).norm()>1e-14){
//            cout << (corr.reshaped()-corr2.reshaped()).norm() << "\n\n";

//            cout << corr << "\n\n" << corr2 << "\n\n";
//        }

        forl(i, corr.rows()){                                                         // change this as 0 frequencies not normalised
            corr(i,i) = 1.0;
        }

        AdjustMatrix(adjust_parameters, v, corr);       // Adjust matrix to take account of different correlation accuracies

        cout << v.size() << '\n';

        VectorXd estimate;

        MatrixXd estimate_multi;

        // if constexpr determines the T and chooses the corresponding function for shrinkage


        if(mcmc_parameters.method == "spikeandslab"){
            estimate = SpikeAndSlab(v, corr, mcmc_parameters, gen);
        }
        else if(mcmc_parameters.method == "multipheno"){
            estimate_multi = MultiPhenoShrink(v, corr,mcmc_parameters, gen);
        }
        else{
            if(adjust_parameters.seed>=0){
                gen.seed(adjust_parameters.seed);
            }

            estimate = CalculateContinuousShrinkageEstimate(v, corr, mcmc_parameters, gen);
        }


        if(outputfile!=""){
            if(mcmc_parameters.method == "multipheno"){
                forc(i, v){

                    string chromosome = v[i].chromosome;

                    if(chromosome=="X") chromosome = "23";

                    Printtabline(output, v[i].rsnum, chromosome, v[i].position, v[i].baseline, v[i].effect, estimate_multi(0,i), estimate_multi(1, i), estimate_multi(2, i), v[i].eaf);
                }
            }
            else{
                forc(i, v){

                    string chromosome = v[i].chromosome;

                    if(chromosome=="X") chromosome = "23";

                    Printtabline(output, v[i].rsnum, chromosome, v[i].position, v[i].baseline, v[i].effect, estimate[i], v[i].eaf);
                }
            }
        }

        if((!testfile.empty() || chrx_files.size()==2) && mcmc_parameters.method != "multipheno"){
            forc(i, v){

                string chromosome = v[i].chromosome;

                if(chromosome=="X") chromosome = "23";

                prsmap[Posclass(stoi(chromosome), v[i].position, v[i].baseline, v[i].effect)] = {estimate[i], v[i].eaf};
            }
        }

        cout << ++num << '\n';
    }

    output.pop();   // make sure output file is written to before generating the PGS scores which can take time

    // need to add an option for regression or possibly not use cv_num


    if((!testfile.empty() || chrx_files.size()==2) && mcmc_parameters.method != "multipheno"){

        return GetPRSstats(testphenotypefile, testfile, bgenix_prefix, chrx_files, analysis_group, prsmap, adjust_parameters.prs_file)[0];
    }
    else{
        return vector<double>(4, 0);
    }
}


// CL_mcmc_runs shrinks the estimates based on parameters passed by the command line
// used by S4_shrink program

int CL_mcmc_runs(int argc, char* argv[])
{
    string inputfile;

    string outputfile;

    double a, b, phi;

    double k = 1.0;

    double inflation = 1;

    int n_tests = 2000;

    string analysis_group;

    vector<string> testfile;

    string bgenix_prefix;

    string testphenotypefile = "";

    string weights_file = "";

    string chrX_weights_file = "";

    string prs_file = "";

    string genmapprefix = "";

    double uncorrdist = -1.0;

    string reference_file_prefix;

    bool correctvalue;

    double nvalue = 0;

    double max_phi = 1.0e50;

    double region_significance = 1e-200;

    double region_factor = 1.0;

    double sparse = 10.0;

    int region_distance = 0;

    bool eaf_variance = false;

    double corr_damp = 1.0;

    vector<string> chrx_files;

    double zdiff = -1.0;

    int seed = -1;

    int threshold_to_increase_runs = -1;

    vector<double> approximation_threshold{-1.0,-1.0};

    vector<double> bounds;

    string method = "horseshoe";

    double min_psi = 0.0;

    vector<double> pvec, wvec;

    double regvalue = 0.0;

    int n_burnin = 1000;

    bool use_ss_r2 = false;

    double eaf_power = 1.0;

    double r2_power = 1.0;

    bool functional_scores = false;

    try
    {
        po::options_description desc;                           // this is a boost class that handles command line parameters

        desc.add_options()("help,h", "produce help message")
                (",a", po::value<double>(&a), " a value")
                (",b", po::value<double>(&b), "b value")
                (",g", po::value<string>(&analysis_group), "analysis group")
                (",p", po::value<double>(&phi), "phi value")
                ("burnin,",po::value<int>(&n_burnin), "Number of burnin runs")
                (",n", po::value<int>(&n_tests), "number of recorded runs")
                (",i", po::value<string>(&inputfile), "input file prefix")
                (",o", po::value<string>(&outputfile)->required(), "output file")
                (",f", po::value<vector<string>>(&testfile)->multitoken(), "test file for checking PRS score")
                ("ref,", po::value<string>(&bgenix_prefix), "bgenix list files containing snp ids")
                (",r", po::value<string>(&reference_file_prefix)->required(), "reference files prefix")
                ("pheno,", po::value<string>(&testphenotypefile), "phenotype file for testing")
                ("correct,", po::value<bool>(&correctvalue)->zero_tokens(), "correct r2 based on variance of statistic")
                ("nvalue,", po::value<double>(&nvalue), "average number value")
                (",w", po::value<string>(&weights_file), "weights file for reference samples")
                ("wchrx,", po::value<string>(&chrX_weights_file), "chrX weights file for reference samples")
                ("o2,", po::value<string>(&prs_file), "output file for prs scores")
                ("inflation,", po::value<double>(&inflation), "inflation factor for statistics (default 1)")
                ("maxphi,", po::value<double>(&max_phi), "maximum phi to protect against numerical instability (default large value)")
                ("region_significance,", po::value<double>(&region_significance), "significance threshold of top SNP for region to be significant")
                ("region_factor,", po::value<double>(&region_factor), "phi multiplicative change if region is significant and SNP is in boundary of top SNP")
                ("region_distance,", po::value<int>(&region_distance), "The distance that the region extends from the top SNP")
                ("sparse,", po::value<double>(&sparse), "values for which correlation is set to be 0")
                ("mapprefix,", po::value<string>(&genmapprefix), "genetic map prefix")
                ("uncorrdist,", po::value<double>(&uncorrdist), "distance at which correlation is assumed to be 0")
                ("eafvariance,", po::value<bool>(&eaf_variance)->zero_tokens(), "assume variance is larger for rarer SNPs")
                ("corrdamp,", po::value<double>(&corr_damp), "damp the correlation by specified amount")
                (",z", po::value<double>(&zdiff), "maximum increase in absolute z to be allowed in MCMC run")
                (",x", po::value<vector<string>>(&chrx_files)->multitoken(), "files used to calibrate individuals if chromosome X files used")
                ("seed,", po::value<int>(&seed), "seed to initialise random number generate, default is a computer generated seed")
                ("runthreshold,", po::value<int>(&threshold_to_increase_runs), "if number of snps in a group is below this threshold then increase runs as less computationally expensive")
                ("approx,", po::value<vector<double>>(&approximation_threshold)->multitoken(), "values for approximation in MCMC calculation (default no approximation)")
                ("psibounds,", po::value<vector<double>>(&bounds)->multitoken(), "Psi bounds for different zscores")
                ("method,", po::value<string>(&method), "Method for shrinkage (default horseshoe)")
                (",k", po::value<double>(&k), "Second shape parameter if using generalised gamma (default = 1)")
                ("minpsi,", po::value<double>(&k), "minimum psi (check if this reduces variance of MCMC sampling")
                ("pvec,", po::value<vector<double>>(&pvec)->multitoken(), "p values for correlated phenotypes if used")
                ("wvec,", po::value<vector<double>>(&wvec)->multitoken(), "w values for correlated phenotypes if used")
                ("reg,", po::value<double>(&regvalue), "regularisation factor for horseshoe model")
                ("ssr2,", po::value<bool>(&use_ss_r2)->zero_tokens(), "Use the estimated summary statistics r2")
                ("eafpower,", po::value<double>(&eaf_power), "the exponent for how the variance is larger for rarer SNPs")
                ("r2power,", po::value<double>(&r2_power), "the exponent for how a more accurate r2 increase the variance in the horseshoe model")
                ("functional,", po::value<bool>(&functional_scores)->zero_tokens(), "use functional scores");

        po::variables_map vm;

        po::store(parse_command_line(argc, argv, desc), vm);

        if (vm.count("help"))
        {
            cout << desc;

            return 1;
        }

        notify(vm);
    }
    catch (const po::error& e)
    {
        cout << e.what() << '\n';

        return -1;
    }

    vector<double> prsstats_sum(8, 0.0);

    while(approximation_threshold.size()<2){
        approximation_threshold.push_back(-1.0);
    }

    if(bounds.empty() && max_phi < 1.0e10) bounds.push_back(max_phi);

    MCMC_parameters mcmc_parameters{n_tests+n_burnin, n_burnin, 1, a, b, phi, inflation, max_phi, region_significance, region_factor, region_distance, sparse, int(eaf_variance), zdiff, threshold_to_increase_runs, approximation_threshold, bounds, method, pvec, wvec, k, min_psi, regvalue, use_ss_r2, eaf_power, r2_power, functional_scores};

    Adjust_parameters adjust_parameters{correctvalue, nvalue, weights_file, prs_file, chrX_weights_file, genmapprefix, uncorrdist, corr_damp, seed};

    auto prsstats = GeneralPRSscores(inputfile, reference_file_prefix, mcmc_parameters, adjust_parameters, analysis_group, outputfile, testfile, bgenix_prefix, chrx_files, testphenotypefile);

    if(!testfile.empty() || chrx_files.size()==2){
        Printtabline(cout, a, b, phi, prsstats[0], prsstats[1], prsstats[2], prsstats[3], prsstats[4], prsstats[0]*prsstats[4], (prsstats[0]-1.96*prsstats[1])*prsstats[4], (prsstats[0]+1.96*prsstats[1])*prsstats[4]);
    }

    return 0;
}


// CL_stepwise_runs generates estimates from a threshold specified by the command line parameters
// used by Stepwise program

int CL_stepwise_runs(int argc, char* argv[])
{
    string inputfile;

    string outputfile;

    double threshold;

    int c = 5;

    string analysis_group;

    vector<string> testfile;

    string bgenix_prefix;

    vector<string> chrx_files;

    string testphenotypefile = "";

    string weights_file = "";

    string reference_file_prefix;

    bool correctvalue;

    double nvalue = 0;


    try
    {
        po::options_description desc;

        desc.add_options()("help,h", "produce help message")
                (",t", po::value<double>(&threshold)->required(), " stepwise z-score threshold")
                (",c", po::value<int>(&c), "number of cv groups")
                (",g", po::value<string>(&analysis_group), "analysis group")
                (",i", po::value<string>(&inputfile), "input file prefix")
                (",o", po::value<string>(&outputfile)->required(), "output file")
                (",f", po::value<vector<string>>(&testfile)->multitoken(), "test file for checking PRS score")
                ("ref,", po::value<string>(&bgenix_prefix), "file containing ids in test file")
                (",r", po::value<string>(&reference_file_prefix)->required(), "reference files prefix")
                ("pheno,", po::value<string>(&testphenotypefile), "phenotype file for testing")
                ("correct,", po::value<bool>(&correctvalue)->zero_tokens(), "correct r2 based on variance of statistic")
                ("nvalue,", po::value<double>(&nvalue), "average number value")
                (",w", po::value<string>(&weights_file), "weights file for reference samples")
                (",x", po::value<vector<string>>(&chrx_files)->multitoken());



        po::variables_map vm;

        po::store(parse_command_line(argc, argv, desc), vm);

        if (vm.count("help"))
        {
            cout << desc;

            return 1;
        }

        notify(vm);
    }
    catch (const po::error& e)
    {
        cout << e.what() << '\n';

        return -1;
    }

    double auc_sum = 0;

    double prschi2sum = 0;

    double effect_sum = 0;

    double se_sum = 0;

    Adjust_parameters adjust_parameters{correctvalue, nvalue, weights_file};

//    auto prsstats = GeneralPRSscores(inputfile, reference_file_prefix, threshold, adjust_parameters, analysis_group, outputfile, testfile, bgenix_prefix, chrx_files, testphenotypefile);

//    if(!testfile.empty() || chrx_files.size()==2){
//        Printtabline(cout, threshold, prsstats[0], prsstats[1], prsstats[2], prsstats[3]);
//    }

    return 0;
}


// CL_PRSstats gives estimates for how well the PRS scores match the phenotype data
// It takes either a bgen and analysis weights file or just a dosage file, together with the phenotypes file
// used by PRSstats program

int CL_PRSstats(int argc, char* argv[])
{
    string analysis_group;

    vector<string> testfile;

    string bgenix_prefix;

    vector<string> chrx_files;

    string testphenotypefile = "";

    vector<string> prs_file;

    try
    {
        po::options_description desc;

        desc.add_options()("help,h", "produce help message")
                (",g", po::value<string>(&analysis_group), "analysis group")
                (",f", po::value<vector<string>>(&testfile)->multitoken(), "test file for checking PRS score")
                (",i", po::value<string>(&bgenix_prefix), "reference files")
                (",w", po::value<vector<string>>(&prs_file)->multitoken(), "weights for prs")
                ("pheno,", po::value<string>(&testphenotypefile), "phenotype file for testing")
                (",x", po::value<vector<string>>(&chrx_files)->multitoken());


        po::variables_map vm;

        po::store(parse_command_line(argc, argv, desc), vm);

        if (vm.count("help"))
        {
            cout << desc;

            return 1;
        }

        notify(vm);
    }
    catch (const po::error& e)
    {
        cout << e.what() << '\n';

        return -1;
    }

    PrintPRSstats(testphenotypefile, testfile, bgenix_prefix, chrx_files, analysis_group, prs_file, -1);

    return 0;
}


// CL_PRSstats gives estimates for how well the PRS scores match the phenotype data
// It takes either a bgen and analysis weights file or just a dosage file, together with the phenotypes file
// used by PRSstats program

int CL_PRSgenerate(int argc, char* argv[])
{
    vector<string> data_files;

    string prs_file;

    string bgenix_prefix;

    string outputfile;

    vector<string> chrx_files;

    string group_file;

    bool dragen = false;

    int chromosome = 0;

    bool use_positions = false;

    bool leading_zero = false;

    bool leading_chr = false;

    bool numeric_x = false;

    try
    {
        po::options_description desc;

        desc.add_options()("help,h", "produce help message")
                (",d", po::value<vector<string>>(&data_files)->multitoken(), "files containing imputed data")
                (",w", po::value<string>(&prs_file), "weights for prs")
                (",r", po::value<string>(&bgenix_prefix), "reference files created by bgenix list")
                (",o", po::value<string>(&outputfile), "output file containing prs weights")
                (",x", po::value<vector<string>>(&chrx_files)->multitoken(), "files to calibrate individuals if chromosome X files used")
                (",g", po::value<string>(&group_file), "group file to average missing values (only needed for possible missingness)")
                ("dragen,", po::value<bool>(&dragen)->zero_tokens(), "get id names from dragen naming convention")
                (",c", po::value<int>(&chromosome), "assume chromosome of data file")
                ("pos,p", po::value<bool>(&use_positions)->zero_tokens(), "use positions rather than ids")
                ("leadzero,", po::value<bool>(&leading_zero)->zero_tokens(), "when using positions assume 1 is 01 etc (used in biobank)")
                ("leadchr,", po::value<bool>(&leading_chr)->zero_tokens(), "when using positions assume 1 is chr1 etc")
                ("numx,", po::value<bool>(&numeric_x)->zero_tokens(), "when using positions assume X is coded as 23");

        po::variables_map vm;

        po::store(parse_command_line(argc, argv, desc), vm);

        if (vm.count("help"))
        {
            cout << desc;

            return 1;
        }

        notify(vm);
    }
    catch (const po::error& e)
    {
        cout << e.what() << '\n';

        return -1;
    }

    auto prsmap = GeneratePRSmap(prs_file);

    bool multipgs = !prsmap.empty() && prsmap.begin()->second.size() > 2;

    if(multipgs){
        auto value = GetMultiplePrincipalComponentValuesFromSeveralBgenFiles(data_files, prsmap, bgenix_prefix, chrx_files, 100);

        cout << "Created prs values\n";

        zofstream output(outputfile);

        forl(i, value.rows()){
            if(i) output << '\t';

            output << "pgs" << i+1;
        }

        output << '\n';

        forl(i, value.cols()){
            forl(j, value.rows()){
                if(j) output << '\t';

                output << value(j,i);
            }

            output << '\n';
        }
    }
    else{
        auto value = GetPrincipalComponentValuesFromSeveralBgenFiles(data_files, prsmap, {100, group_file, dragen, use_positions, leading_zero, leading_chr, numeric_x, chromosome}, bgenix_prefix, chrx_files);

        cout << "Created prs values\n";

        zofstream output(outputfile);

        output << "prs\n";

        forc(i, value){
            output << value[i] << '\n';
        }
    }

    return 0;
}


// CL_GenerateIds generates ids for the S4 selected SNPs by using bgenix files
// It takes a list of bgenix files together with the S4 snps file and the output file

int CL_GenerateIds(int argc, char* argv[])
{
    string bgenix_prefix;

    string snp_file;

    string outputfile;

    try
    {
        po::options_description desc;

        desc.add_options()("help,h", "produce help message")
                (",b", po::value<string>(&bgenix_prefix), "bgenix list files")
                (",s", po::value<string>(&snp_file), "file generated by S4_select")
                (",o", po::value<string>(&outputfile), "output file containing prs weights");

        po::variables_map vm;

        po::store(parse_command_line(argc, argv, desc), vm);

        if (vm.count("help"))
        {
            cout << desc;

            return 1;
        }

        notify(vm);
    }
    catch (const po::error& e)
    {
        cout << e.what() << '\n';

        return -1;
    }

    set<Posclass> snpset;

    zifstream input(snp_file);

    Posclass posclass;

    while(skipdelim(input, '\t', 6)){
        input >> posclass >> skipline;

        snpset.insert(posclass);

        if(posclass.ref.size()!=1 || posclass.alt.size()!=1){
            swap(posclass.ref, posclass.alt);

            snpset.insert(posclass);
        }
    }

    auto snplist = ConvertPositionsToIds(bgenix_prefix, snpset);

    zofstream output(outputfile);

    for(auto& v:snplist){
        output << v << '\n';
    }

    return 0;
}

void CompareReads()
{
    BgenParser parser1("/home/jpt34/USrefPanel/us_chrall23.bgen");

    BgenParser parser2("/home/jpt34/USrefPanel/us_chrall23.bgen");

//    auto include = FileRead<vector<int>>("/home/jpt34/Downloads/ocac_include.txt");

    Imputeclass imputeclass1, imputeclass2;

    while(parser1.ReadAllImputeclass(imputeclass1)){
        parser2.ReadAllImputeclassCompare(imputeclass2);

        if(abs(imputeclass1.r2-imputeclass2.r2)>1.0e-12){
            cout << imputeclass1.rsnum << ' '  << imputeclass2.rsnum << ' ' << imputeclass1.eaf - imputeclass2.eaf << ' ' << imputeclass1.r2 - imputeclass2.r2 << '\n';
        }

        forl(i, imputeclass1.size()){
            if(imputeclass1[i] != imputeclass2[i]){
                cout << imputeclass1.rsnum << ' '  << imputeclass2.rsnum << ' ' << i << ' ' << imputeclass1[i] << ' ' << imputeclass2[i] << '\n';
            }
        }
    }
}


bool DampCorrelationMatrix(SpMat& corr, const VectorXd& bstandard, const VectorXd& beta_standard, const SumStatsParameters& params)
{
    vector<int> indices_to_change;

    forc(i, beta_standard){

        if(abs(beta_standard[i])>abs(bstandard[i]) && abs(beta_standard[i]-bstandard[i])>10){
            indices_to_change.push_back(i);
        }
    }

    for(auto & v: indices_to_change){
        corr.row(v) *= params.corr_damp;
        corr.col(v) *= params.corr_damp;
        corr.coeffRef(v,v) = 1.0;
    }

    return !indices_to_change.empty();
}


// Calculate PRS from matrix from summary statisitics

vector<double> MatrixPRSstats(SpMat& corrmat, const VectorXd& w, const VectorXd& b, const VectorXd& s)
{
    forc(i, b){
        if(s[i]<1e-5 || !isfinite(s[i]) || !isfinite(b[i]) || !isfinite(w[i])){
            cout << i << ' ' << s[i] << ' ' << b[i] << ' ' << w[i] << '\n';
        }
    }

    DiagonalMatrix<double, -1> smat = s.asDiagonal();

    SpMat c = smat.inverse() * corrmat * smat.inverse();

    double numerator = (w.cwiseQuotient(s.cwiseAbs2())).dot(b);

    double quadprod = (w.transpose() * c * w).value();

    double uncorrelated_denom = w.cwiseQuotient(s).squaredNorm();

    return {numerator, quadprod, uncorrelated_denom};
}


void TestSummaryPRS()
{
    SpMat corrmat(2,2);

    using Tlet = Eigen::Triplet<double>;

    vector<Tlet> triplet;

    triplet.push_back(Tlet(0,0,1));
    triplet.push_back(Tlet(1,1,1));
    triplet.push_back(Tlet(0,1,0.7));
    triplet.push_back(Tlet(1,0,0.7));

    corrmat.setFromTriplets(triplet.begin(), triplet.end());

    VectorXd b{{0.1,0.05}};

    VectorXd s{{0.01,0.01}};

    VectorXd w{{1,1}};

//    auto prs_stats = MatrixPRSstats(corrmat, w, b, s);

//    cout << prs_stats.first << '\t' << prs_stats.second << '\n';
}


// Use summary statistics to predict the accuracy of PRS weights
// prs_file is the same as the file used to generate PRS weights


void EstimatePRSaccuracyFromSummaryStatistics(const string &summarystats_file, const string& bgenix_reffile, const string& ref_file, const string prs_file, const string& include_file, string geneticmap_file, const SumStatsParameters& params)
{
    auto prsmap = GeneratePRSmap(prs_file);

    vector<string> rsids=[&]{
        if(bgenix_reffile!=""){
            auto unique_ids=ConvertPositionsToIds(bgenix_reffile, prsmap, false);

            return vector<string>(unique_ids.begin(), unique_ids.end());
        }
        else{
            return vector<string>();
        }
    }();

    auto position_map = GetPositionsFromIds(rsids, ref_file);

    auto summary_stats = GetSelectedSNPstatistics(prsmap, position_map, summarystats_file, ref_file, include_file, geneticmap_file);

    int ndim = summary_stats.size();

    using Tlet = Eigen::Triplet<double>;

    vector<Tlet> triplet;

    forl(i, ndim){
        triplet.push_back(Tlet(i,i,1));
        forl(j, i){
            if(abs(summary_stats[i].genetic_position - summary_stats[j].genetic_position) < params.genetic_distance){
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

    VectorXd b(ndim), s(ndim), w(ndim);

    forc(i, b){
        b[i] = summary_stats[i].b;
        s[i] = summary_stats[i].s;
        w[i] = summary_stats[i].w[0];

        if(!isfinite(s[i])){
            cout << summary_stats[i].rsnum << ' ' << summary_stats[i].posclass << '\n';
        }
    }

    auto prs_stats = MatrixPRSstats(corrmat, w, b , s);

    double denom = params.corr_damp*(prs_stats[1]-prs_stats[2])+prs_stats[2];

    cout << prs_stats[0]/denom << '\t' << 1/sqrt(denom) << '\t' << Sqr(prs_stats[0])/denom << '\t' << prs_stats[0] << '\t' << prs_stats[1] << '\t' << prs_stats[2] << '\n';
}


// Command line to use summary statistics to predict the accuracy of PRS weights

int CL_SummaryStatisticsPRSstats(int argc, char* argv[])
{
    string summarystats_file;

    string bgenix_reffile;

    string bgen_ref_file;

    string prs_file;

    string include_file;

    string geneticmap_file;

    double genetic_dist;

    double corr_damp = 1.0;

    try
    {
        po::options_description desc;

        desc.add_options()("help,h", "produce help message")
                (",s", po::value<string>(&summarystats_file), "summary statistics file for phenotype of interest")
                (",b", po::value<string>(&bgen_ref_file), "bgen reference data file")
                (",r", po::value<string>(&bgenix_reffile), "reference files created by bgenix list")
                (",w", po::value<string>(&prs_file), "weights for prs")
                (",i", po::value<string>(&include_file), "include file")
                ("genmapfile,", po::value<string>(&geneticmap_file), "genetic map file")
                ("dist,", po::value<double>(&genetic_dist), "max distance to consider being in same correlated group")
                ("corrdamp,", po::value<double>(&corr_damp), "damp the correlation by specified amount");


        po::variables_map vm;

        po::store(parse_command_line(argc, argv, desc), vm);

        if (vm.count("help"))
        {
            cout << desc;

            return 1;
        }

        notify(vm);
    }
    catch (const po::error& e)
    {
        cout << e.what() << '\n';

        return -1;
    }

    EstimatePRSaccuracyFromSummaryStatistics(summarystats_file, bgenix_reffile, bgen_ref_file, prs_file, include_file, geneticmap_file, {genetic_dist, corr_damp});

    return 0;
}


int CL_CalcSNPcontributions(int argc, char* argv[])
{
    string bgenix_reffile_prefix;

    string bgen_ref_file_prefix;

    string pgm_file;

    string include_file;

    string geneticmap_file;

    double genetic_dist;

    try
    {
        po::options_description desc;

        desc.add_options()("help,h", "produce help message")
                (",b", po::value<string>(&bgen_ref_file_prefix), "bgen reference data file")
                (",r", po::value<string>(&bgenix_reffile_prefix), "reference files created by bgenix list")
                (",w", po::value<string>(&pgm_file), "weights for pgm")
                (",i", po::value<string>(&include_file), "include file")
                ("genmapfile,", po::value<string>(&geneticmap_file), "genetic map file")
                ("dist,", po::value<double>(&genetic_dist), "max distance to consider being in same correlated group");


        po::variables_map vm;

        po::store(parse_command_line(argc, argv, desc), vm);

        if (vm.count("help"))
        {
            cout << desc;

            return 1;
        }

        notify(vm);
    }
    catch (const po::error& e)
    {
        cout << e.what() << '\n';

        return -1;
    }

    CalcSNPcontributions(pgm_file, bgenix_reffile_prefix, bgen_ref_file_prefix, include_file, geneticmap_file, genetic_dist);

    return 0;
}


int CL_MultiPGSsumstats(int argc, char* argv[])
{
    string summarystats_file;

    string bgenix_reffile;

    string bgen_ref_file;

    vector<string> pgm_files;

    string include_file;

    string chrX_include_file;

    string geneticmap_file;

    string corr_file_prefix = "";

    double genetic_dist = 3.0;

    try
    {
        po::options_description desc;

        desc.add_options()("help,h", "produce help message")
                (",s", po::value<string>(&summarystats_file)->required(), "summary statistics file for phenotype of interest")
                (",b", po::value<string>(&bgen_ref_file), "bgen reference data file")
                (",r", po::value<string>(&bgenix_reffile), "reference files created by bgenix list")
                (",w", po::value<vector<string>>(&pgm_files)->multitoken()->required(), "weights for polygenic models")
                (",i", po::value<string>(&include_file), "include file")
                ("i2,", po::value<string>(&chrX_include_file), "include file for chromosome X (only specify if different from main include file)")
                ("genmapfile,", po::value<string>(&geneticmap_file), "genetic map file")
                ("dist,", po::value<double>(&genetic_dist), "max distance to consider being in same correlated group")
                ("corr,c", po::value<string>(&corr_file_prefix), "correlation file");


        po::variables_map vm;

        po::store(parse_command_line(argc, argv, desc), vm);

        if (vm.count("help"))
        {
            cout << desc;

            return 1;
        }

        notify(vm);
    }
    catch (const po::error& e)
    {
        cout << e.what() << '\n';

        return -1;
    }

    if(corr_file_prefix!=""){
        MultiPGSregressionFromCorrelationFiles(summarystats_file, pgm_files, corr_file_prefix);
    }
    else{
        MultiPGSregression(summarystats_file, pgm_files, bgenix_reffile, bgen_ref_file, include_file, chrX_include_file, geneticmap_file, genetic_dist);
    }

    return 0;
}


int CL_CorrelationMatrix(int argc, char* argv[])
{
    string summarystats_file;

    string bgenix_reffile;

    string bgen_ref_file;

    vector<string> pgm_files;

    string include_file;

    string chrX_include_file;

    string geneticmap_file;

    string output_file;

    double genetic_dist = 3.0;

    try
    {
        po::options_description desc;

        desc.add_options()("help,h", "produce help message")
                (",b", po::value<string>(&bgen_ref_file), "bgen reference data file")
                (",r", po::value<string>(&bgenix_reffile), "reference files created by bgenix list")
                (",w", po::value<vector<string>>(&pgm_files)->multitoken(), "weights for polygenic models")
                (",i", po::value<string>(&include_file), "include file")
                ("genmapfile,", po::value<string>(&geneticmap_file), "genetic map file")
                ("dist,", po::value<double>(&genetic_dist), "max distance to consider being in same correlated group")
                (",o", po::value<string>(&output_file), "output file containing correlations");


        po::variables_map vm;

        po::store(parse_command_line(argc, argv, desc), vm);

        if (vm.count("help"))
        {
            cout << desc;

            return 1;
        }

        notify(vm);
    }
    catch (const po::error& e)
    {
        cout << e.what() << '\n';

        return -1;
    }

    PrintCorrelationMatrix(bgenix_reffile, bgen_ref_file, pgm_files, include_file, geneticmap_file, genetic_dist, output_file);

    return 0;
}


void TestMahalobisDistance()
{
    int ndim = 5000;

    MatrixXd a = MatrixXd::Random(ndim, ndim);

    VectorXd w = VectorXd::Random(ndim);

    MatrixXd r = a*a.transpose();

    MachineTimer time4;

    MatrixXd rinv = r.llt().solve(MatrixXd::Identity(ndim, ndim));

//    MatrixXd rinv = r.inverse();

    time4.ResetPrintElapsed();

    MachineTimer time1;

    double result = double(w.transpose() * rinv * w);

    time1.ResetPrintElapsed();

    MachineTimer time5;

    LLT<MatrixXd> rllt(r);

    time5.ResetPrintElapsed();

    MachineTimer time2;

    VectorXd t = rllt.solve(w);

    double result2 = w.dot(t);

    time2.ResetPrintElapsed();

    MatrixXd transform = rllt.matrixL();

    MachineTimer time3;

    double result3 = transform.triangularView<Lower>().solve(w).norm();

//    double result3 = t2.squaredNorm();

    time3.ResetPrintElapsed();

    cout << result << '\n';

    cout << result2 << '\n';

    cout << Sqr(result3) << '\n';

}

int main(int argc, char* argv[])
{
//    TestMahalobisDistance();

//    TestSummaryPRS();

//    CompareReads();
//    MatrixXd b = MatrixXd::Random(5,5);

//    MatrixXd a = b*b.transpose();

//    MatrixXd ainv = a.inverse();

//    VectorXd c = VectorXd::Random(5);

//    double value1 = c.transpose()*ainv*c;

//    VectorXd transform = a.llt().solve(c);

//    double value2 = c.dot(transform);


//    cout << value1 << '\n' << value2 << '\n';



//    zifstream input("/home/jpt34/BCAC/chr22_selected_snps.txt");

//    SplitStringView splitview;

//    while(input >> splitview){
//        cout << splitview.m_linesplit[0] << '\n';
//    }

//    return CL_GenerateIds(argc, argv);

//    CL_stepwise_runs(argc, argv);

    return CL_CondAnalysis(argc, argv);

//    return CL_bfdp_runs(argc, argv);

//    return CL_mcmc_runs(argc, argv);

//    return CL_SelectedLS(argc, argv);

//    return CL_PRSstats(argc, argv);

//    return CL_CalcSNPcontributions(argc, argv);

//    return CL_MultiPGSsumstats(argc, argv);

//    return CL_CorrelationMatrix(argc, argv);

//    return CL_PRSgenerate(argc, argv);

//    return CL_SummaryStatisticsPRSstats(argc, argv);

//    genfile::bgen::SqliteIndexQuery query("height_r02_all.bgen.bgi");

//    vector<string> rsids = {"rs6690515"};

//    query.include_rsids(rsids);

//    query.initialise();
}

