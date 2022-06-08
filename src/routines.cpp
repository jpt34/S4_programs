#include "routines.h"

#include <Eigen/Dense>
#include <array>
#include <boost/algorithm/string.hpp>
#include <boost/math/distributions/chi_squared.hpp>
#include "Mlogit.h"
#include <jlibrary.h>
#include <zstream.h>
#include "bgen/parser.hpp"

using namespace std;
namespace bo = boost;


// This functions gets covariate data from the phenotype file
// The second row contains how the field in the phenotype file is used


void GetCovariateDataFromFileWithStudyInfo(const string& inputfile, ArrayXi& status, VectorXd& regval, MatrixXd& cov, ArrayXi& study, vector<string>& studyconv, const string& group)
{
    ifstream input(inputfile);

    CheckOpen(input, inputfile);

    SplitString splitstr;

    input >> splitstr;

    string tempstr;

    getline(input, tempstr);

    if(tempstr.find(".") != string::npos)
    {
        cout << "There is a string containing \".\" in the second line\nThis may indicate a decimal number placed there by mistake\n";

        exit(1);
    }

    stringstream strm(tempstr);

    vector<string> parameters((istream_iterator<string>(strm)), istream_iterator<string>());

    if(splitstr.size() != parameters.size())
    {
        cout << "The number of fields in the first row is " << splitstr.size() << " and the number of fields in the second row is " << parameters.size() << '\n';

        cout << "The data may be space-separated rather than tab-separated\n";

        exit(1);
    }

    int ndata = 0;

    while(input.good())
    {
        input.ignore(1e9,'\n');

        input >> ws;

        ndata++;
    }

    input.clear();

    input.seekg(0,ios::beg);

    for(int i=0;i<2;i++) input.ignore(1e9,'\n');

    int ncontcov = count(parameters.begin(), parameters.end(), "2") + count(parameters.begin(), parameters.end(), "C");

    int ncatcov = count(parameters.begin(), parameters.end(), "discrete");

    int nstudies = count(parameters.begin(), parameters.end(), "3") + count(parameters.begin(), parameters.end(), "S");

    int nassoc = count(parameters.begin(), parameters.end(), "1") + count(parameters.begin(), parameters.end(), "B") + count(parameters.begin(), parameters.end(), "M");

    int nregval = count(parameters.begin(), parameters.end(), "4") + count(parameters.begin(), parameters.end(), "P");


    if(nassoc > 1)
    {
        forc(i, parameters)
        {
            if(splitstr[i]!=group && (parameters[i]=="1" || parameters[i]=="B" || parameters[i]=="M"))
            {
                parameters[i] = "0";

                nassoc--;
            }
        }
    }

    if(nregval > 1)
    {
        forc(i, parameters)
        {
            if(splitstr[i]!=group && (parameters[i]=="4" || parameters[i]=="P"))
            {
                parameters[i] = "0";

                nregval--;
            }
        }
    }

    if(nassoc + nregval != 1)
    {
        cout << "Either no associations to test or too many\n";

        cout << "You may need to specify a group if more than one possible association\n";

        exit(1);
    }

    if(nstudies > 1)
    {
        cout << "Can only analyse one study - the rest should be coded as dummy variables\n";

        exit(1);
    }

    MatrixXd contcov;

    MatrixXi catcov;

    double *pcontcov = 0;

    int *pcatcov = 0;

    int *pstatus = 0;

    double *pregval = 0;

    if(nassoc > 0)
    {
        status.resize(ndata);

        pstatus = &status[0];
    }

    if(nregval > 0){
        regval.resize(ndata);

        pregval = &regval[0];
    }

    if(ncontcov > 0)
    {
        contcov.resize(ncontcov, ndata);

        pcontcov = &contcov(0,0);
    }

    if(ncatcov > 0)
    {
        catcov.resize(ncatcov, ndata);

        pcatcov = &catcov(0,0);
    }

    vector<string> studystr = (nstudies?vector<string>(ndata,""):vector<string>(ndata,"Default"));

    vector<string> linesplit;

    for(int i=0;i<ndata;i++)
    {
        getline(input, tempstr);

        bo::algorithm::split(linesplit, tempstr, bo::algorithm::is_any_of("\t"));

        if(tempstr.size()==0)
        {
            cout << "Empty line in middle of phenotype file\n";

            exit(1);
        }

        for(auto& v: linesplit) bo::algorithm::trim(v);

        if(linesplit.size() < parameters.size())		// if case of errors in creating tab separated file
        {
            for(int j=linesplit.size();j<parameters.size();j++) linesplit.push_back("");
        }

        for(int j=0;j<parameters.size();j++)
        {
//			input >> tempstr >> ws;

            if(parameters[j]=="1" || parameters[j]=="B")
            {
                try
                {
                    *pstatus = stoi(linesplit[j]);

//                    if(*pstatus>1) *pstatus=1;

                    if(*pstatus<0) *pstatus=-1;

                    pstatus++;
                }
                catch(...)
                {
                    *pstatus++ = -1;
                }
            }
            else if(parameters[j]=="2" || parameters[j]=="C")
            {
                try
                {
                    *pcontcov = stod(linesplit[j]);

                    pcontcov++;
                }
                catch(...)
                {
                    *pcontcov++ = -99;
                }
            }
            else if(parameters[j]=="discrete")
            {
                try
                {
                    *pcatcov = stoi(linesplit[j]);

                    pcatcov++;
                }
                catch(...)
                {
                    *pcatcov++ = -99;
                }
            }
            else if(parameters[j]=="4" || parameters[j] == "P")
            {
                try
                {
                    *pregval = stod(linesplit[j]);

                    pregval++;
                }
                catch(...)
                {
                    *pregval++ = -99;
                }
            }
            else if(parameters[j]=="3" || parameters[j] == "S")
            {
                studystr[i] = linesplit[j];
            }
            else
            {
                // don't try to convert string
            }
        }
    }

    MatrixXi convcov;

    MatrixXd tempcov;

    if(ncatcov > 0) ConvertCategoricalVariables(catcov, convcov);

    if(ncontcov + ncatcov > 0) tempcov.resize(contcov.rows() + convcov.rows(), ndata);

    if(ncontcov > 0) tempcov.topRows(contcov.rows()) = contcov;

    if(ncatcov > 0) tempcov.bottomRows(convcov.rows()) = convcov.cast<double>();

// Remove variables that don't vary on non-missing data

    vector<int> includeindices;

    for(int i=0;i<tempcov.rows();i++)
    {
        double diff;

        double max = -1e5;

        double min = 1e5;

        for(int j=0;j<tempcov.cols();j++)
        {
            if(status.size()!=0){
                if(status[j] > -0.5)
                {
                    if(tempcov(i,j) > max) max = tempcov(i,j);

                    if(tempcov(i,j) < min && tempcov(i,j) > -98) min = tempcov(i,j);
                }
            }
            else{
                if(regval[j] != -99)
                {
                    if(tempcov(i,j) > max) max = tempcov(i,j);

                    if(tempcov(i,j) < min && tempcov(i,j) > -98) min = tempcov(i,j);
                }
            }

            diff = max-min;
        }

        if(diff > 1.0e-5)
        {
            includeindices.push_back(i);
        }
    }

    if(!includeindices.empty())
    {
        cov.resize(includeindices.size(), ndata);

        for(int i=0;i<includeindices.size();i++)
        {
            cov.row(i) = tempcov.row(includeindices[i]);
        }
    }

    if(cov.rows() > 0)
    {
        for(int i=0;i<status.size();i++)
        {
            if((abs(cov.col(i).array()+99)<1.0e-5).any())
            {
                status[i] = -1;
            }
        }

        for(int i=0;i<regval.size();i++)
        {
            if((abs(cov.col(i).array()+99)<1.0e-5).any())
            {
                regval[i] = -99;
            }
        }
    }

    if(status.size()!=0){
        ConvertToStudyNumbers(status, studystr, study, studyconv);
    }
    else{
        ConvertToStudyNumbersForRegressionValues(regval, studystr, study, studyconv);
    }
}


// this checks whether there has been a problem in reading either the genotype file or the phenotype file

void CheckValidStatus(const vector<int>& statustemp, int ndata)
{
    if(statustemp.size() != ndata)
    {
        cout << "Mismatch in number of individuals in phenotype and genotype files\n";

        cout << statustemp.size() << " in phenotype file, " << ndata << " in genotype file\n";

        exit(1);
    }
    else
    {
        int nvalid = isum(statustemp, a);

        if(nvalid==0)
        {
            cout << "No valid individuals\n";

            exit(1);
        }
    }
}


// works out how many individuals are in the data file

int GetNumIndividuals(const string& inputfile)
{   
    int result;

    if(inputfile.size()>=4 && inputfile.substr(inputfile.size()-4)=="bgen")
    {
        BgenParser bgenParser(inputfile);

        result = bgenParser.number_of_samples();
    }
    else
    {
        zifstream input(inputfile);

        string tempstr;

        getline(input, tempstr);

        stringstream strm(tempstr);

        int ndata = 0;

        while(strm >> tempstr)
        {
            ndata++;
        }

        CheckOpen(input, inputfile);

        if(ndata==6)
        {
            result = tempstr.size();
        }
        else if((ndata-5)%3!=0 && (ndata-6)%3!=0)
        {
            cout << "Number of columns is " << ndata << " which does not correspond to a valid number of individuals in an IMPUTE file\n";

            exit(1);
        }
        else
        {
            result = (ndata-5)/3;   // if extra field result is still the same by rounding
        }
    }

    return result;
}


// Gets a vector of cv groups from the phenotype file

vector<int> GetCvGroup(const string& inputfile){
    zifstream input(inputfile);

    SplitString splitstr;

    input >> splitstr;

    int indexnum = -1;

    forc(i, splitstr){
        if(splitstr[i]=="Group" || splitstr[i]=="group"){
            indexnum = i;
        }
    }

    if(indexnum ==-1){
        return vector<int>();
    }
    else{
        input >> splitstr;

        vector<int> result;

        while(input >> splitstr){
            result.push_back(stoi(splitstr[indexnum]));
        }

        return result;
    }
}


// Gets phenotype data from file and put into structure

PhenotypeData GetPhenotypeData(const string& inputfile, const string& group, int cv_group)
{
    ArrayXi status;

    VectorXd regval;

    MatrixXd cov;

    ArrayXi study;

    vector<string> studyconv;

    GetCovariateDataFromFileWithStudyInfo(inputfile, status, regval, cov, study, studyconv, group);

    PhenotypeData phenotypedata;

    if(status.size()!=0){
        phenotypedata.include.resize(status.size());
    }
    else{
        phenotypedata.include.resize(regval.size());
    }

    auto cvvec = GetCvGroup(inputfile);

    int nvalid = 0;

    if(status.size()>0){
        forc(i, status){
            if(status[i]>=0 && (cvvec.size()!=status.size() || cvvec[i]==cv_group || cv_group==-1)){
                ++nvalid;
            }
        }
    }
    else{
        forc(i, regval){
            if(regval[i]!=-99 && (cvvec.size()!=regval.size() || cvvec[i]==cv_group || cv_group==-1)){
                ++nvalid;
            }
        }
    }

    if(nvalid==0)
    {
        cout << "No studies have both cases and controls\n";

        exit(1);
    }

    int num = 0;

    if(status.size()>0){
        phenotypedata.status.resize(nvalid);
    }
    else{
        phenotypedata.regval.resize(nvalid);
    }

    phenotypedata.cov.resize(cov.rows(), nvalid);

    phenotypedata.study.resize(nvalid);

    phenotypedata.studyconv = studyconv;


    if(status.size()>0){
        forc(i, status)
        {
            if(status[i] >= 0 && (cvvec.size()!=status.size() || cvvec[i]==cv_group || cv_group==-1))
            {
                phenotypedata.include[i] = 1;

                phenotypedata.status[num] = status[i];

                if(cov.rows()>0) phenotypedata.cov.col(num) = cov.col(i);

                phenotypedata.study[num] = study[i];

                num++;
            }
            else
            {
                phenotypedata.include[i] = 0;
            }
        }
    }
    else{
        forc(i, regval)
        {
            if(regval[i] != -99 && (cvvec.size()!=regval.size() || cvvec[i]==cv_group || cv_group==-1))
            {
                phenotypedata.include[i] = 1;

                phenotypedata.regval[num] = regval[i];

                if(cov.rows()>0) phenotypedata.cov.col(num) = cov.col(i);

                phenotypedata.study[num] = study[i];

                num++;
            }
            else
            {
                phenotypedata.include[i] = 0;
            }
        }
    }


    set<int> invalid_indices;

    forl(i, cov.rows())
    {
        if(abs(phenotypedata.cov.row(i).maxCoeff()-phenotypedata.cov.row(i).minCoeff()) < 1e-9) // if no change of covariate on filtered data
        {
             invalid_indices.insert(i);
        }
    }

    if(!invalid_indices.empty())
    {
        MatrixXd newcov(phenotypedata.cov.rows()-invalid_indices.size(), phenotypedata.cov.cols());

        for(int i=0,j=0;i<phenotypedata.cov.rows();i++)
        {
            if(!invalid_indices.count(i))
            {
                newcov.row(j) = phenotypedata.cov.row(i);

                j++;
            }
        }

        phenotypedata.cov = newcov;
    }

    return phenotypedata;
}


// This converts a categorical variable into an integer vector starting from 0 to number_used-1.

void ConvertCategoricalVariables(const MatrixXi& catcov, MatrixXi& convcov)
{
    VectorXi maxvalues = catcov.rowwise().maxCoeff();

    VectorXi minvalues = catcov.rowwise().minCoeff();

    int nconv = maxvalues.sum() - minvalues.sum();

    int offset = 0;

    convcov.resize(nconv, catcov.cols());

    for(int i=0;i<catcov.rows();i++)
    {
        for(int j=minvalues[i]+1;j<=maxvalues[i];j++)
        {
            convcov.row(offset+j-minvalues[i]-1) = (catcov.row(i).array()==j).cast<int>();
        }

        offset += maxvalues[i]-minvalues[i];
    }
}


// This converts study strings to values starting from 0, provided that at least one case and control is in the study set

template<typename Derived>
void ConvertToStudyNumbers(ArrayBase<Derived>& status, const vector<string>& studystr, ArrayXi& study, vector<string>& studyconv)
{
    map<string, std::array<int,2>> counts;

    for(int i=0;i<status.size();i++)
    {
        if(status[i]!=-1) ++counts[studystr[i]][(status[i]!=0)];
    }

    int num = 0;

    map<string, int> convert;

    for(auto& v: counts)
    {
        if(v.second[0] > 0 && v.second[1] > 0 && v.first!="" && v.first!="NA" && v.first!="-99" && v.first!="\"\"" && v.first!="\"NA\"" && v.first!="\"-99\"")
        {
            studyconv.push_back(v.first);
            convert[v.first] = num++;
        }
        else
        {
            convert[v.first] = -1;
        }
    }

    study.resize(status.size());

    for(int i=0;i<status.size();i++)
    {
        if(convert[studystr[i]]==-1)
        {
            study[i] = 0;
            status[i] = -1;	// make sure not analysed for small numbers in study
        }
        else
        {
            study[i] = convert[studystr[i]];
        }
    }
}


// This converts study strings to values starting from 0 based on regression values rather than binary outcomes

void ConvertToStudyNumbersForRegressionValues(VectorXd& regval, const vector<string>& studystr, ArrayXi& study, vector<string>& studyconv)
{
    map<string, int> counts;

    for(int i=0;i<regval.size();i++)
    {
        if(regval[i]!=-99) ++counts[studystr[i]];
    }

    int num = 0;

    map<string, int> convert;

    for(auto& v: counts)
    {
        if(v.second>1 && v.first!="" && v.first!="NA" && v.first!="-99" && v.first!="\"\"" && v.first!="\"NA\"" && v.first!="\"-99\"")
        {
            studyconv.push_back(v.first);
            convert[v.first] = num++;
        }
        else
        {
            convert[v.first] = -1;
        }
    }

    study.resize(regval.size());

    for(int i=0;i<regval.size();i++)
    {
        if(convert[studystr[i]]==-1)
        {
            study[i] = 0;
            regval[i] = -99;	// make sure not analysed for small numbers in study
        }
        else
        {
            study[i] = convert[studystr[i]];
        }
    }
}


// converts matrix to a form suitable for sse calculations

void ConvertToSseMatrix(const MatrixXd& cov, MatrixXd& cov_sse)
{
    cov_sse.resize(cov.rows()*2, (cov.cols()+1)/2);

    for(int indiv=0;indiv<2*(cov.cols()/2);indiv+=2)
    {
        forl(i, cov.rows())
        {
            cov_sse(i*2, indiv/2) = cov(i, indiv);

            cov_sse(i*2+1, indiv/2) = cov(i, indiv+1);
        }
    }

    if(cov.cols()%2)
    {
        forl(i, cov.rows())
        {
            cov_sse(i, cov_sse.cols()-1) = cov(i, cov.cols()-1);
        }
    }
}


map<double, double> ReadGeneticMap(const string& inputfile)
{
    if(inputfile=="") return map<double, double>();

    zifstream input(inputfile);

    input >> skipline;

    double f1, f2, f3;

    map<double, double> result;

    while(input >> f1 >> f2 >> f3){
        result[f1]=f3;
    }

    return result;
}


double CalcGeneticPosition(const map<double, double> genmap, double physical_position)
{
    if(genmap.empty()) return physical_position;

    auto itor = genmap.lower_bound(physical_position);

    double result = 0;

    if(itor!=genmap.end()){

        if(itor==genmap.begin()){
            result = itor->second;
        }
        else{
            auto itorprev = prev(itor);

            result = (physical_position-itorprev->first)*(itor->second-itorprev->second)/(itor->first-itorprev->first) + itorprev->second;
        }
    }
    else{
        result = genmap.rbegin()->second;
    }

    return result;
}


std::istream& operator>>(std::istream& is, SplitStringView& v)
{
    int n = 0;

    if(std::getline(is, v.m_linestr))
    {
//        if(v.merge_delim==boost::algorithm::token_compress_on)
//        {
//            boost::algorithm::trim(v.linestr);
//        }

        size_t first = 0;

        v.m_linesplit.resize(0);

        experimental::string_view vstr(v.m_linestr);

        while (first < vstr.size())
        {
            const auto second = vstr.find_first_of(v.m_delim, first);

            if (first != second)
                v.m_linesplit.emplace_back(vstr.substr(first, second-first));

            if (second == std::string_view::npos)
                break;

            first = second + 1;
        }
    }

    return is;
}


std::ostream& operator<<(std::ostream& os, const SplitStringView& v)
{
    os << v.m_linestr;

    return os;
}


std::istream& skipdelim(std::istream& is, char c, int numtimes)
{
    for(int i=0;i<numtimes;i++){
        is.ignore(std::numeric_limits<int>::max(), c);
    }

    return is;
}
