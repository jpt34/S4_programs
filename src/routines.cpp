#include "routines.h"

#include <Eigen/Dense>
#include <array>
#include <boost/algorithm/string.hpp>
#include <boost/math/distributions/chi_squared.hpp>
#include "Mlogit.h"
#include <jlibrary.h>
#include <zstream.h>
#include "bgen/parser.hpp"
#include "bgen/IndexQuery.hpp"

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

PhenotypeData GetPhenotypeData(const string& inputfile, const string& group)
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

    int nvalid = 0;

    if(status.size()>0){
        forc(i, status){
            if(status[i]>=0){
                ++nvalid;
            }
        }
    }
    else{
        forc(i, regval){
            if(regval[i]!=-99){
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
            if(status[i] >= 0)
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
            if(regval[i] != -99)
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

    string f2;

    double f1, f3;

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


map<Posclass, vector<double>> GeneratePRSmap(const string& inputfile)
{
    SplitString splitstr;

    zifstream input(inputfile);

    map<Posclass, vector<double>> result;

    vector<double> vec;

    while(input >> splitstr){
        if(splitstr.size()==7){
            result[Posclass(splitstr[1]=="X"?23:stoi(splitstr[1]), stoi(splitstr[2]), splitstr[3], splitstr[4])] = {stod(splitstr[5]), stod(splitstr[6])};
        }
        else{
            vec.resize(splitstr.size()-5);

            forc(i,vec){
                vec[i] = stod(splitstr[i+5]);

                result[Posclass(splitstr[1]=="X"?23:stoi(splitstr[1]), stoi(splitstr[2]), splitstr[3], splitstr[4])] = vec;
            }
        }
    }

    return result;
}


map<Posclass, long> GetPositionsFromIds(vector<string> rsids, const string& ref_file)
{
    BgenParser bgenParser(ref_file);

    genfile::bgen::SqliteIndexQuery query(ref_file +".bgi");

    query.include_rsids(rsids);

    query.initialise();

    map<Posclass, long> result;

    for(size_t i=0;i<query.number_of_variants();i++){

        auto fileposition = query.locate_variant(i);

        bgenParser.JumpToPosition(fileposition.first);

        std::string chromosome ;
        uint32_t position ;
        std::string rsid ;
        std::vector< std::string > alleles ;
        std::vector< std::vector< double > > probs ;

        bgenParser.read_variant( &chromosome, &position, &rsid, &alleles );

        result[Posclass(ConvertChrStr(chromosome), position, alleles[0], alleles[1])] = fileposition.first;
    }

    return result;
}


// Get SNP statisitics of snps that are in the reference file
// Column 2 chromosome
// Column 3 position
// Column 4 baseline
// Column 5 effect
// Column 8 odds ratio
// Column 9 standard error

vector<SumStatsPRS> GetSelectedSNPstatistics(const map<Posclass, vector<double>>& prsmap, const map<Posclass, long>& position_map, const string& summarystats_file, const string& ref_file, const string& include_file, const string& geneticmap_file)
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

    zifstream input(summarystats_file);

    while(input >> splitstr){
        try {
            Posclass posclass(ConvertChrStr(string(splitstr[1])), stoi(string(splitstr[2])), splitstr[4], splitstr[3]);

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
                            double bvalue = flip * stod(string(splitstr[7]));

                            double svalue = stod(string(splitstr[8]));

                            if(isfinite(bvalue) && isfinite(svalue) && svalue>0){
                                results.push_back(SumStatsPRS());

                                results.back().w = itor_prsmap->second;

                                results.back().posclass = posclass;

                                results.back().b = bvalue;

                                results.back().s = svalue;

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

                                int flip2 = MatchAlleles(imputeclass, FlipCheck{itor_prsmap->first.ref, itor_prsmap->first.alt});

                                if(flip2 == -1){
                                    for(auto& value:results.back().w)
                                    value *= -1;    // won't matter if frequency as not used here
                                }
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

// This estimates the PRS based on the weightings given by the data structure
// normally only one file is used but this routine was orignally written for several bgen files

VectorXd GetPrincipalComponentValuesFromSeveralBgenFiles(const vector<string>& datafile, map<Posclass, vector<double>>& data, const GenerateParameters& generate_parameters, const string& bgenix_prefix, const vector<string>& chrx_files)
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

    auto value = rsids.empty() && generate_parameters.use_positions ?
                    GetPrincipalComponentValuesFromPositions(datafile, data, generate_parameters):
                    GetPrincipalComponentValuesFromIds(datafile, data, rsids, generate_parameters);

    if(chrx_files.size()==2){
        auto value2 = rsids.empty() && generate_parameters.use_positions ?
                        GetPrincipalComponentValuesFromPositions({chrx_files[0]}, data, generate_parameters):
                        GetPrincipalComponentValuesFromIds({chrx_files[0]}, data, rsids, generate_parameters);

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

            value[v.first] += value2[v.second];
        }
    }

    return value;
}


std::vector<std::string> getVariantStrings(
    const std::map<Posclass, std::vector<double>>& variantMap,
    int targetChromosome) // e.g., 1 (for chr1)
{
    std::vector<std::string> result;

    for (const auto& [posClass, scores] : variantMap) {
        if (posClass.chr == targetChromosome || targetChromosome <=0) {

             result.push_back(
                 "chr" + std::to_string(posClass.chr) + ":" +
                 std::to_string(posClass.position) + ":" +
                 posClass.ref + ":" +
                 posClass.alt
             );
        }
    }

    return result;
}



// This estimates the PRS based on the weightings given by the data structure
// normally only one file is used but this routine was orignally written for several bgen files

VectorXd GetPrincipalComponentValuesFromIds(const vector<string>& datafile, map<Posclass, vector<double>>& data, const vector<string>& rsids, const GenerateParameters &generate_parameters)
{
    string linestr, rsnum;

    Posclass posclass;

    Imputeclass imputeclass;

    uint32_t position;

    vector<string> alleles;

    vector<vector<double>> probs;

    int nsnps = 0;

    vector<double> prsvec;

    auto group_vec = (generate_parameters.group_file==""?vector<int>():FileRead<vector<int>>(generate_parameters.group_file));

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

            if(group_vec.empty()){
                bgenparse.UpdatePRS(prsvec, factor, itor->second[0], 2*itor->second[1]);
            }
            else{
                bgenparse.UpdatePGSMissingAdjusted(prsvec, group_vec, factor, itor->second[0], 2*itor->second[1]);
            }



            //            double value;

            //            forc(i, probs)
            //            {
            //                if(probs[i].size()==3)
            //                {
            //                    value = 2*probs[i][2]+probs[i][1];
            //                }
            //                else if(probs[i].size()==2)
            //                {
            //                    value = 2*probs[i][1];
            //                }
            //                else
            //                {
            //                    value = -1;
            //                }

            //                if(value >= -1e-9)
            //                {
            //                    double meandiff = 2*(factor==-1) + factor*value - 2*itor->second[1];

            //                    prs[i] += meandiff * itor->second[0];
            //                }
            //            }

            if(generate_parameters.print){
                ++nsnps;

                if(nsnps%generate_parameters.print==0){
                    cout << nsnps << " " << imputeclass.rsnum << '\n';
                }
            }
        }
        else
        {
            bgenparse.ignore_probs();
        }
    };

    int chromosome = 0;

    for(auto& v: datafile)
    {
        BgenParser bgenParser(v);

        if(!bgenParser.is_bgenformat){
            cerr << v << " is not in BGEN format\n";

            exit(1);
        }

        if(!rsids.empty()){
            genfile::bgen::SqliteIndexQuery query(Paste(v,".bgi"));

            query.include_rsids(rsids);

            query.initialise();

            if(generate_parameters.print) cout << query.number_of_variants() << '\n';

            for(size_t i=0;i<query.number_of_variants();i++){

                auto fileposition = query.locate_variant(i);

                bgenParser.JumpToPosition(fileposition.first);

                bgenParser.read_variant(&imputeclass.chromosome, &position, &imputeclass.rsnum, &alleles);

                UpdatePRSdosages(bgenParser);
            }
        }
        else if(generate_parameters.dragen){

            // Determine chromosome
            chromosome = (datafile.size() == 22) ? chromosome + 1 : generate_parameters.chromosome;

            set<string> found_snps;
            vector<string> dragen_ids = getVariantStrings(data, chromosome);

            // Track matches
            size_t original_matches = 0;
            size_t swapped_matches = 0;
            size_t flipped_matches = 0;

            // First pass - original alleles
            auto first_pass = [&]() {
                genfile::bgen::SqliteIndexQuery query(Paste(v,".bgi"));
                query.include_rsids(dragen_ids);
                query.initialise();
                original_matches = query.number_of_variants();

                for(size_t i = 0; i < original_matches; i++) {
                    auto fileposition = query.locate_variant(i);
                    bgenParser.JumpToPosition(fileposition.first);
                    bgenParser.read_variant(&imputeclass.chromosome, &position, &imputeclass.rsnum, &alleles);
                    found_snps.insert(imputeclass.rsnum);
                    UpdatePRSdosages(bgenParser);
                }
            };

            // Lambda to process missing variants with optional flipping
            auto process_missing_variants = [&](const vector<string>& missing_ids, bool flip_alleles = false) {
                vector<string> processed_ids;
                map<string, string> processed_to_original;

                for (const auto& id : missing_ids) {
                    // Parse the variant ID
                    vector<string> parts;
                    size_t start = 0, end = id.find(':');
                    while (end != string::npos) {
                        parts.push_back(id.substr(start, end - start));
                        start = end + 1;
                        end = id.find(':', start);
                    }
                    parts.push_back(id.substr(start));

                    if (parts.size() == 4) {
                        string ref = parts[2];
                        string alt = parts[3];

                        string processed_id;

                        if (flip_alleles) {
                            // Skip A/T and C/G combinations
                            if (!((ref == "A" && alt == "T") || (ref == "T" && alt == "A") ||
                                 (ref == "C" && alt == "G") || (ref == "G" && alt == "C"))) {
                                // Flip the alleles
                                auto flip = [](const string& a) {
                                    if (a == "A") return string("T");
                                    if (a == "T") return string("A");
                                    if (a == "C") return string("G");
                                    if (a == "G") return string("C");
                                    return a;
                                };

                                processed_id = parts[0] + ":" + parts[1] + ":" + flip(ref) + ":" + flip(alt);
                                processed_ids.push_back(processed_id);
                                processed_to_original[processed_id] = id;

                                processed_id = parts[0] + ":" + parts[1] + ":" + flip(alt) + ":" + flip(ref);
                                processed_ids.push_back(processed_id);
                                processed_to_original[processed_id] = id;
                            }
                        } else {
                            // Simple swap of ref/alt
                            processed_id = parts[0] + ":" + parts[1] + ":" + alt + ":" + ref;
                            processed_ids.push_back(processed_id);
                            processed_to_original[processed_id] = id;
                        }
                    }
                }

                // Query the processed variants
                genfile::bgen::SqliteIndexQuery query(Paste(v, ".bgi"));
                query.include_rsids(processed_ids);
                query.initialise();

                vector<string> found_variants;
                for (size_t i = 0; i < query.number_of_variants(); i++) {
                    auto fileposition = query.locate_variant(i);
                    bgenParser.JumpToPosition(fileposition.first);
                    bgenParser.read_variant(&imputeclass.chromosome, &position, &imputeclass.rsnum, &alleles);

                    UpdatePRSdosages(bgenParser);
                    found_variants.push_back(processed_to_original[imputeclass.rsnum]);
                }

                return found_variants;
            };

            // First pass - original alleles
            first_pass();

            original_matches = found_snps.size();

            vector<string> final_missing;

            // Second pass - swapped alleles for missing variants
            if (found_snps.size() < dragen_ids.size()) {
                vector<string> missing_ids;
                for (const auto& id : dragen_ids) {
                    if (!found_snps.count(id)) {
                        missing_ids.push_back(id);
                    }
                }

                auto new_found_list = process_missing_variants(missing_ids);

                swapped_matches = new_found_list.size();

                set<string> new_found(new_found_list.begin(), new_found_list.end());

                // Third pass - flipped alleles for still missing variants
                if (new_found.size() < missing_ids.size()) {
                    vector<string> still_missing;
                    for (const auto& id : missing_ids) {
                        if (!new_found.count(id)) {
                            still_missing.push_back(id);
                        }
                    }

                    auto flipped_found = process_missing_variants(still_missing, true);

                    set<string> new_flipped_found(new_found_list.begin(), new_found_list.end());

                    flipped_matches = flipped_found.size();

                    // Third pass - flipped alleles for still missing variants
                    if (flipped_found.size() < still_missing.size()) {
                        for (const auto& id : still_missing) {
                            if (!new_flipped_found.count(id)) {
                                final_missing.push_back(id);
                            }
                        }
                    }
                }
            }

            if(generate_parameters.print) {
                cout << "Variant matching results:\n"
                     << "  Original matches: " << original_matches << '\n'
                     << "  Swapped matches: " << swapped_matches << '\n'
                     << "  Flipped matches: " << flipped_matches << '\n'
                     << "  Total snps: " << original_matches+swapped_matches+flipped_matches << "/" << dragen_ids.size() << '\n';

                if(!final_missing.empty())
                {
                    cout << "Missing variants\n";

                    for(auto v:final_missing){
                        cout << v << '\n';
                    }
                }
            }
        }
        else{
            while(bgenParser.read_variant(&imputeclass.chromosome, &position, &imputeclass.rsnum, &alleles))
            {
                UpdatePRSdosages(bgenParser);
            }
        }
    }

    VectorXd prs(prsvec.size());

    forc(i, prsvec){
        prs[i] = prsvec[i];
    }

    return prs;
}


VectorXd GetPrincipalComponentValuesFromPositions(const vector<string>& datafile, map<Posclass, vector<double>>& data, const GenerateParameters& generate_parameters)
{
    string linestr, rsnum;

    Posclass posclass;

    Imputeclass imputeclass;

    uint32_t position;

    vector<string> alleles;

    vector<vector<double>> probs;

    int nsnps = 0;

    vector<double> prsvec;

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

            bgenparse.UpdatePRS(prsvec, factor, itor->second[0], 2*itor->second[1]);

            //            forc(i, imputeclass){
            //                double meandiff = 2*(factor==-1) + factor*imputeclass[i] - 2*itor->second.back();

            //                Map<VectorXd> eigenVector(itor->second.data(), itor->second.size()-1);

            //                prsmat.col(i) += meandiff*eigenVector;
            ////                forl(j, itor->second.size()-1){
            ////                    prsmat(j,i) += meandiff * itor->second[j];
            ////                }
            //            }

            if(generate_parameters.print){
                ++nsnps;

                if(nsnps%generate_parameters.print==0){
                    cout << nsnps << " " << imputeclass.rsnum << '\n';
                }
            }
        }
        else
        {
            bgenparse.ignore_probs();
        }
    };

    auto Format_chromosome = [&](int chr) -> std::string {
        std::string result;

        // Handle special chromosomes (X, Y, M/MT)
        if (chr == 23 || chr == 24 || chr == 25) {
            if (generate_parameters.numeric_x) {
                result = std::to_string(chr);
            }
            else{
                result = "X";
            }
        }
        // Handle regular chromosomes (1-22)
        else {
            if (generate_parameters.leading_zero && chr < 10) {
                result = "0" + std::to_string(chr);
            } else {
                result = std::to_string(chr);
            }
        }

        // Add 'chr' prefix if needed
        if (generate_parameters.leading_chr) {
            result = "chr" + result;
        }

        return result;
    };

    forc(i, datafile) {
        const int chromosome = generate_parameters.chromosome > 0 ? generate_parameters.chromosome : i + 1;

        // First identify all positions we need to query
        std::vector<std::tuple<std::string, uint32_t>> ranges;
        for (auto& v : data) {
            if (v.first.chr == chromosome || (datafile.size() != 22 && generate_parameters.chromosome==0)) {
                ranges.emplace_back(Format_chromosome(v.first.chr), v.first.position);
            }
        }

        const size_t BATCH_SIZE = 900;
        const size_t total_variants = ranges.size();
        size_t processed = 0;

        if (generate_parameters.print) {
            cout << "Total variants to process: " << total_variants << '\n';
        }

        while (processed < total_variants) {
            // Create new parser and query for each batch
            BgenParser bgenParser(datafile[i]);
            genfile::bgen::SqliteIndexQuery query(Paste(datafile[i], ".bgi"));

            // Add current batch of ranges
            const size_t batch_end = std::min(processed + BATCH_SIZE, total_variants);
            for (size_t j = processed; j < batch_end; j++) {
                query.include_range({get<0>(ranges[j]), get<1>(ranges[j]), get<1>(ranges[j])});
            }

            query.initialise();
            const size_t variants_in_batch = query.number_of_variants();

//            if (generate_parameters.print) {
//                cout << "Processing batch: " << processed << " to " << (processed + variants_in_batch - 1)
//                     << " (" << variants_in_batch << " variants)" << '\n';
//            }

            for (size_t j = 0; j < variants_in_batch; j++) {
                auto fileposition = query.locate_variant(j);
                bgenParser.JumpToPosition(fileposition.first);
                bgenParser.read_variant(&imputeclass.chromosome, &position, &imputeclass.rsnum, &alleles);
                UpdatePRSdosages(bgenParser);
            }

            processed = batch_end;
        }
    }


    VectorXd prs(prsvec.size());

    forc(i, prsvec){
        prs[i] = prsvec[i];
    }

    return prs;
}

// Get logistic statistics for calculated PRS values

vector<vector<double>> GetPRSchi2(const PhenotypeData& phenotypedata, const vector<VectorXd>& p)
{
    Mlogit logit(phenotypedata, p.size());

    MatrixXd cov(phenotypedata.cov.rows()+p.size(), phenotypedata.cov.cols());

    cov << MatrixXd::Zero(p.size(), phenotypedata.cov.cols()), phenotypedata.cov;

    forc(j, p){

// This checks if the two vectors have the same size. If they don’t, it means there is a mismatch between the weightings file and the phenotype file and the program terminates.
        if(p[j].size()!=phenotypedata.include.size()){
            cout << "Size of weightings file is " << p[j].size() << " and of phenotype file is " << phenotypedata.include.size() << "\n";

            exit(1);
        }

        int count = 0;

        forc(i, phenotypedata.include){
            if(phenotypedata.include[i]){
                cov(j, count) = p[j][i];
                ++count;
            }
        }
    }

    double lik;

    VectorXd x;

    MatrixXd hessian;

    logit.LogitStats(cov, lik, x, hessian);

    vector<vector<double>> result;

    forc(i, p){
        result.push_back({x[i], sqrt(hessian(i,i)), 2*(lik-logit.liknull)});
    }

    return result;
}


// Get regression statistics for calculated PRS values

vector<vector<double>> GetRegressionPRS(const PhenotypeData& phenotypedata, const vector<VectorXd>& p)
{
    int nstudies = phenotypedata.study.maxCoeff()+1;

    int nparam = phenotypedata.cov.rows() + p.size() + nstudies;

    MatrixXd cov(nparam, phenotypedata.cov.cols());

    cov << MatrixXd::Zero(p.size(),phenotypedata.cov.cols()), MatrixXd::Zero(nstudies, phenotypedata.cov.cols()), phenotypedata.cov;

    forc(j, p){

// This checks if the two vectors have the same size. If they don’t, it means there is a mismatch between the weightings file and the phenotype file and the program terminates.
        if(p[j].size()!=phenotypedata.include.size()){
            cout << "Size of weightings file is " << p[j].size() << " and of phenotype file is " << phenotypedata.include.size() << "\n";

            exit(1);
        }

        int count = 0;

        forc(i, phenotypedata.include){
            if(phenotypedata.include[i]){
                cov(j, count) = p[j][i];
                cov(p.size()+phenotypedata.study[count], count)=1;
                ++count;
            }
        }
    }

    MatrixXd var = (cov*cov.transpose());

    vector<vector<double>> result;

    if(var.determinant()!=0){
        var = var.inverse();

        VectorXd c = cov*phenotypedata.regval;

        VectorXd b = var*c;

        double regsum = phenotypedata.regval.squaredNorm();

        double rss = regsum - b.dot(c);

        forc(i, p){
            double se = sqrt(rss*var(i,i)/(phenotypedata.study.size()-cov.rows()));

            result.push_back({b[i], se, Sqr(b[i]/se)});
        }

//        VectorXd pred = cov.transpose()*b;

//        double rss = (phenotypedata.regval - pred).squaredNorm();

        int nrowcalc = cov.rows()-p.size();

        MatrixXd covnull = cov.bottomRows(nrowcalc);

        MatrixXd varnull = (covnull*covnull.transpose()).inverse();

        VectorXd cnull = covnull*phenotypedata.regval;

        VectorXd bnull = varnull*cnull;

        double mss = regsum - bnull.dot(cnull);

//        VectorXd pred2 = cov.bottomRows(nrowcalc).transpose()*b.tail(nrowcalc);

//        double mss = (phenotypedata.regval - pred2).squaredNorm();

        double r2 = 1-rss/mss;

        for(auto& v: result){
            v.push_back(r2);
        }

        return result;
    }
    else{
        return vector<vector<double>>(p.size(), vector<double>{0.0,0.0,0.0});
    }
}


// Get stratified area under the curve that takes into account the study groups

double GetStratRankSumStat(const VectorXi& status, const VectorXi& study, const VectorXd& p)
{
    map<int, map<double, vector<int>>> rankp;

    forc(i, p)
    {
        if(status[i]==0 || status[i]==1) rankp[study[i]][p[i]].push_back(status[i]);
    }

    double weightsum = 0.0;

    double ustrata_strat = 0.0;

    for(auto& x: rankp)
    {
        double ranksum[2] = {0,0};

        int studycounts[2] = {0,0};

        double rankvalue = 0.5;

        for(auto& v: x.second)
        {
            for(auto& w: v.second)
            {
                ranksum[w] += rankvalue + 0.5*v.second.size();

                ++studycounts[w];
            }

            rankvalue += v.second.size();
        }

        double weight = (double(studycounts[0])*studycounts[1])/(studycounts[0]+studycounts[1]+1);

        weightsum += weight;

        if(weight > 0)
        {
            ustrata_strat += weight*(ranksum[1] - studycounts[1]*double(studycounts[1]+1)/2.0)/(double(studycounts[0])*studycounts[1]);
        }
    }

    ustrata_strat /= weightsum;

    return ustrata_strat;
}


// Get pearson correlation between prs values and outcome

double GetPearsonCorrelation(const VectorXd& regval, const VectorXd& prs)
{
    VectorXd regval_std = (regval.array()-regval.mean());

    regval_std /= regval_std.norm();

    VectorXd prs_std = (prs.array()-prs.mean());

    regval_std /= prs_std.norm();

    return regval_std.dot(prs_std);
}
