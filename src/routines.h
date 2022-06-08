#ifndef ROUTINES_H
#define ROUTINES_H

#include <experimental/string_view>
#include <Eigen/Dense>
#include <jlibrary.h>
#include <VectorclassNew3/vectorclass.h>


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

PhenotypeData GetPhenotypeData(const string& inputfile, const string &group, int cv_group);

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


#endif // ROUTINES_H
