#ifndef ROUTINES_H
#define ROUTINES_H

#include <experimental/string_view>
#include <Eigen/Dense>
#include <jlibrary.h>
#include <VectorclassNew4/vectorclass.h>
#include <Posclass.h>

using namespace std;
using namespace Eigen;


struct PhenotypeData
{
    vector<int> include;
    VectorXi status;
    VectorXd regval;
    MatrixXd cov;
    VectorXi study;
    vector<string> studyconv;
};


struct SumStatsPRS{
    string rsnum;
    Posclass posclass;
    double genetic_position;
    vector<double> w;   // weights from PRS weights file
    double b;   // summary statistics coefficient estimate
    double s;   // summary statistics standard error estimate
    VectorXd ref_dosages; // reference dosages for SNP
    int corr_index;
};


struct GenerateParameters
{
    int print;
    string group_file;
    bool dragen;
    bool use_positions;
    bool leading_zero = false;
    bool leading_chr = false;
    bool numeric_x = false;
    int chromosome=-1;
};


struct SplitStringView
{
    std::string m_linestr;
    std::vector<std::experimental::string_view> m_linesplit;
    std::string m_delim;

    SplitStringView(const std::string& vdelim="\t"): m_delim(vdelim){}

//    char& operator[](int index){return linesplit[index];}

    const std::experimental::string_view& operator[](int index) const {return m_linesplit[index];}

    std::string str() const {return m_linestr;}

    size_t size() const {return m_linesplit.size();}

    std::vector<std::experimental::string_view>::iterator begin(){return m_linesplit.begin();}

    std::vector<std::experimental::string_view>::const_iterator cbegin() const {return m_linesplit.cbegin();}

    std::vector<std::experimental::string_view>::iterator end(){return m_linesplit.end();}

    std::vector<std::experimental::string_view>::const_iterator cend() const {return m_linesplit.cend();}
};


// adapted from https://www.bfilipek.com/2018/07/string-view-perf-followup.html



void GetCovariateDataFromFileWithStudyInfo(const string& inputfile, ArrayXi& status, VectorXd& regval, MatrixXd& cov, ArrayXi& study, vector<string>& studyconv, const string& group);

void CheckValidStatus(const vector<int> &statustemp, int ndata);

int GetNumIndividuals(const string& inputfile);

vector<int> GetCvGroup(const string& inputfile);

PhenotypeData GetPhenotypeData(const string& inputfile, const string &group);

void ConvertCategoricalVariables(const MatrixXi& catcov, MatrixXi& convcov);

void ConvertToSseMatrix(const MatrixXd& cov, MatrixXd& cov_sse);

map<double, double> ReadGeneticMap(const string& inputfile);

double CalcGeneticPosition(const map<double, double> genmap, double physical_position);


template<typename Derived>
void ConvertToStudyNumbers(ArrayBase<Derived>& status, const vector<string>& studystr, ArrayXi& study, vector<string>& studyconv);

void ConvertToStudyNumbersForRegressionValues(VectorXd& regval, const vector<string>& studystr, ArrayXi& study, vector<string>& studyconv);

std::istream& operator>>(std::istream& is, SplitStringView& v);

std::ostream& operator<<(std::ostream& os, const SplitStringView& v);

std::istream& skipdelim(std::istream& is, char c, int numtimes);

map<Posclass, vector<double>> GeneratePRSmap(const string& inputfile);

map<Posclass, long> GetPositionsFromIds(vector<string> rsids, const string& ref_file);

vector<SumStatsPRS> GetSelectedSNPstatistics(const map<Posclass, vector<double>>& prsmap, const map<Posclass, long>& position_map, const string& summarystats_file, const string& ref_file, const string& include_file, const string& geneticmap_file);

VectorXd GetPrincipalComponentValuesFromSeveralBgenFiles(const vector<string>& datafile, map<Posclass, vector<double>>& data, const GenerateParameters& generate_parameters, const string& bgenix_prefix="", const vector<string>& chrx_files={});

VectorXd GetPrincipalComponentValuesFromPositions(const vector<string>& datafile, map<Posclass, vector<double>>& data, const GenerateParameters &generate_parameters);

std::vector<std::string> getVariantStrings(const std::map<Posclass, std::vector<double>>& variantMap, int targetChromosome);

VectorXd GetPrincipalComponentValuesFromIds(const vector<string>& datafile, map<Posclass, vector<double>>& data, const vector<string>& rsids, const GenerateParameters& generate_parameters);

vector<vector<double> > GetPRSchi2(const PhenotypeData& phenotypedata, const vector<VectorXd> &p);

vector<vector<double> > GetRegressionPRS(const PhenotypeData& phenotypedata, const vector<VectorXd> &p);

double GetStratRankSumStat(const VectorXi& status, const VectorXi& study, const VectorXd& p);

double GetPearsonCorrelation(const VectorXd& regval, const VectorXd& prs);


template<class T>
set<string> ConvertPositionsToIds(const string& bgenix_prefix, T& prsmap, bool multi_file = true)
{
    set<string> result;

    map<Posclass, string> partial_match;

    set<Posclass> exact_match;

    Posclass posclass;

    SplitStringView splitstr;

    forv(i, 1, (multi_file?24:2)){
        ifstream bgenix_input(multi_file?Paste(bgenix_prefix,i,".txt"):bgenix_prefix);

        if(!bgenix_input) continue;

        while(bgenix_input >> splitstr){
            if(splitstr.size()>=7){
                try {
                    posclass = Posclass(ConvertChrStr(string(splitstr[2])), stoi(string(splitstr[3])), splitstr[5], splitstr[6]);

                    //                    posclass = Posclass(ConvertChrStr(splitstr[2]), stoi(splitstr[3]), splitstr[5], splitstr[6]);


                    if(prsmap.count(posclass)){
                        result.emplace(splitstr[1]);
                        exact_match.emplace(posclass);
                    }
                    else{
                        if(splitstr[5].size()!=1 || splitstr[6].size()!=1){
                            posclass = Posclass(ConvertChrStr(string(splitstr[2])), stoi(string(splitstr[3])), splitstr[6], splitstr[5]);

                            if(prsmap.count(posclass)){
                                partial_match.emplace(posclass, splitstr[1]);
                            }
                        }
                    }
                } catch (...) {
                    //                    cerr << "Error reading line containing\n" << splitstr.str() << "\n";
                }
            }
        }
    }

    for(auto& [v,w]:partial_match){
        if(!exact_match.count(v)){
            result.emplace(w);
        }
    }

    return result;
}

//typename map<Posclass, T>::const_iterator

template<class T>
auto GetPosItor(const T& posmap, Posclass posclass){
    auto itor = posmap.find(posclass);

    if(itor==posmap.end() && (posclass.ref.size()!=1 || posclass.alt.size()!=1)){
        swap(posclass.ref, posclass.alt);
        itor = posmap.find(posclass);
    }

    return itor;
}

#endif // ROUTINES_H
