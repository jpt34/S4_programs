#pragma once

#include <iostream>
#include <set>
#include <map>
#include <vector>
#include <algorithm>
#include <numeric>
#include <string>
#include <fstream>
#include <sstream>
#include <cmath>
#include <cstdarg>
#include <boost/algorithm/string.hpp>
#include <boost/foreach.hpp>
#include <boost/spirit/include/qi.hpp>

#include <emmintrin.h>

#ifdef __SSSE3__
#include <tmmintrin.h>
#endif

#define INPUT(input, inputfile) std::ifstream input(inputfile);CheckOpen(input, inputfile)

#define OUTPUT(output, outputfile) std::ofstream output(outputfile);CheckOpen(output, outputfile)


#define forc(i, vec) for(int i=0, f=int(vec.size());i<f;i++)

#define forl(i, num) for(int i=0;i<num;i++)

#define forv(i, start, finish) for(int i=start;i<finish;i++)

#define msort(vec, ...) std::sort(vec.begin(), vec.end(), [&](decltype(*vec.cbegin())& a, decltype(*vec.cbegin())& b){return __VA_ARGS__;})

#define mtransform(vec, vec2, command) std::transform(vec.cbegin(), vec.cend(), vec2.begin(), [&](decltype(*vec.cbegin())& a){return command;})

#define mtransform2(vec, vec2, vec3, command) std::transform(vec.cbegin(), vec.cend(), vec2.cbegin(), vec3.begin(), [&](decltype(*vec.cbegin())& a, decltype(*vec2.cbegin())& b){return command;})

#define mcount_if(vec, command) std::count_if(vec.begin(), vec.end(), [&](decltype(*vec.begin())& a){return command;})

#define mremove_if(vec, command) std::remove_if(vec.begin(), vec.end(), [&](decltype(*vec.begin())& a){return command;})

#define mremovemap_if(vec, command) for(auto i = vec.begin();(i = std::find_if(i, vec.end(), [&](decltype(*vec.begin())& a){return command;})) != vec.end();vec.erase(i++))

#define isum(vec, command) std::accumulate(vec.begin(), vec.end(), 0, [&](int& b, decltype(*vec.begin())& a){return b + command;})

#define dsum(vec, command) std::accumulate(vec.begin(), vec.end(), 0.0, [&](double& b, decltype(*vec.begin())& a){return b + command;})

#define foreach BOOST_FOREACH

template<class InputIterator1, class InputIterator2, class Fn2>
inline Fn2 foreach2_impl(InputIterator1 first1, InputIterator1 last1, InputIterator2 first2, Fn2 fn)
{
	for (; first1 != last1; ++first1, ++first2)
	{
		fn(*first1, *first2);
	}
	return fn;
}

template<class InputIterator1, class InputIterator2, class InputIterator3, class Fn3>
inline Fn3 foreach3_impl(InputIterator1 first1, InputIterator1 last1, InputIterator2 first2, InputIterator3 first3, Fn3 fn)
{
    for (; first1 != last1; ++first1, ++first2, ++first3)
    {
        fn(*first1, *first2, *first3);
    }
    return fn;
}

#define foreach2(vec1, vec2, command) foreach2_impl(vec1.begin(), vec1.end(), vec2.begin(), [&](decltype(*vec1.begin())& a, decltype(*vec2.begin())& b){command;})

#define for2(vec1, vec2) foreach2_impl(vec1.begin(), vec1.end(), vec2.begin(), [&](decltype(*vec1.begin())& a, decltype(*vec2.begin())& b)

#define for3(vec1, vec2, vec3) foreach3_impl(vec1.begin(), vec1.end(), vec2.begin(), vec3.begin(), [&](decltype(*vec1.begin())& a, decltype(*vec2.begin())& b, decltype(*vec3.begin())& c)

template<class T> inline T Sqr(T x) {return x*x;}

template<class T>
double CalcShannonInfo(const T& freq)
{
    auto Xlogx = [](double x)
    {
        return (x<=0 || x>=1 ? 0 : x*log(x));
    };

    auto N = freq.size();

    auto x = dsum(freq, Xlogx(a/N) + Xlogx((1-a)/N));

    auto z = -log(double(N));

    auto probsum = dsum(freq, a)/N;

    auto z2 = Xlogx(probsum) + Xlogx(1-probsum);

    auto shannon_info = x - z - z2;

    return shannon_info;
}

// might not need namespace std but need to check
namespace std
{
template<class T, class U>
std::istream& operator>>(std::istream& is, std::pair<T,U>& v)
{
	is >> v.first >> v.second;

	return is;
}

template<class T, class U>
std::ostream& operator<<(std::ostream& os, const std::pair<T,U>& v)
{
    os << v.first << '\t' << v.second;

    return os;
}

template<std::size_t> struct int_{};

template <class Tuple, size_t Pos>
std::istream& read_tuple(std::istream& in, Tuple& t, int_<Pos> ) {
in >> std::get< std::tuple_size<Tuple>::value-Pos >(t);
return read_tuple(in, t, int_<Pos-1>());
}

template <class Tuple>
std::istream& read_tuple(std::istream& in, Tuple& t, int_<1> ) {
return in >> std::get<std::tuple_size<Tuple>::value-1>(t);
}

template <class... Args>
std::istream& operator>>(std::istream& in, std::tuple<Args...>& t) {
return read_tuple(in, t, int_<sizeof...(Args)>());
}

template <class Tuple, size_t Pos>
std::ostream& print_tuple(std::ostream& out, const Tuple& t, int_<Pos> ) {
out << std::get< std::tuple_size<Tuple>::value-Pos >(t) << '\t';
return print_tuple(out, t, int_<Pos-1>());
}

template <class Tuple>
std::ostream& print_tuple(std::ostream& out, const Tuple& t, int_<1> ) {
return out << std::get<std::tuple_size<Tuple>::value-1>(t);
}

template <class... Args>
std::ostream& operator<<(std::ostream& out, const std::tuple<Args...>& t) {
return print_tuple(out, t, int_<sizeof...(Args)>());
}
}

template< class T, std::size_t N>
std::istream& operator>>(std::istream& in, std::array<T,N>& value)
{
    for(int i=0;i<int(N);i++)
    {
        in >> value[i];
    }

    return in;
}

template< class T, std::size_t N>
std::ostream& operator<<(std::ostream& out, const std::array<T,N>& value)
{
    for(int i=0;i<int(N);i++)
    {
        if(i) out << '\t';

        out << value[i];
    }

    return out;
}

template<class T>
T FileRead(const std::string& inputfile)
{
	std::ifstream input(inputfile);

	if(!input)
	{
		std::cout << "Could not open file " << inputfile << '\n';
	}

	return T((std::istream_iterator<typename T::value_type>(input)), std::istream_iterator<typename T::value_type>());
}


template<class T>
T FileReadH(const std::string& inputfile)
{
	std::ifstream input(inputfile);

	if(!input)
	{
		std::cout << "Could not open file " << inputfile << '\n';
	}

	input.ignore(1e9,'\n');

	return T((std::istream_iterator<typename T::value_type>(input)), std::istream_iterator<typename T::value_type>());
}


template<class T>
bool Getlinesplit(std::istream& input, std::string& linestr, T& linesplit, const std::string& delim = "\t")
{
	if(getline(input, linestr))
	{
		boost::algorithm::split(linesplit, linestr, boost::algorithm::is_any_of(delim));

		return true;
	}
	else
	{
		return false;
	}
}


template <class T, class U>
auto Find(T& data, const U& finddata)->decltype(&data.begin()->second)
{
	auto itor = data.find(finddata);
	
	if(itor != data.end())
	{
		return &itor->second;
	}
	else
	{
		return 0;
	}
}


template<class T, class U>
bool Contains(const T& data, const U& finddata)
{
	return data.find(finddata) != data.end();
}


template<class T, class U>
bool InRegion(const T& data, const U& finddata)
{
	auto p = data.lower_bound(finddata);

	return (p!=data.end() && p->second <= finddata);
}


template<class T>
void CheckOpen(const T& filestream, const std::string& filename)
{
	if(!filestream)
	{
		std::cout << "Could not open " << filename << '\n';

		exit(1);
	}
}


template<class T>
void CheckFileError(const T& filestream, const std::string& filename)
{	
	if(filestream.bad())
	{
		std::cout << "I/0 error in " << filename << '\n';
		exit(1);
	}
	else if(!filestream.eof() && filestream.fail())
	{
		std::cout << "Error in processing " << filename << '\n';
		exit(1);
	}
}

template<typename T>
void Printdelim(std::ostream& output, char u, const T& arg2)
{
    output << arg2;
}

template<typename T, typename... Args> void Printdelim(std::ostream& output, char u, const T& arg2, const Args&...args)
{
    output << arg2 << u;

    Printdelim(output, u, args...);
}

template<typename T, typename... Args> void Printdelimline(std::ostream& output, char u, const T& arg2, const Args&...args)
{
    Printdelim(output, u, arg2, args...);

    output << '\n';
}

template<typename T>
void Printtabline(std::ostream& output, const T& arg2)
{
    output << arg2 << '\n';
}

template<typename T, typename... Args> void Printtabline(std::ostream& output, const T& arg2, const Args&...args)
{
    output << arg2 << '\t';

    Printtabline(output, args...);
}


template<typename T>
void PasteHelper_(std::ostream& output, const T& arg2)
{
    output << arg2;
}

template<typename T, typename... Args>
void PasteHelper_(std::ostream& output, const T& arg2, const Args&...args)
{
    output << arg2;

    PasteHelper_(output, args...);
}

template<typename... Args>
std::string Paste(const Args&...args)
{
    std::ostringstream strm;

    PasteHelper_(strm, args...);

    return strm.str();
}

template<class T>
struct PushStream_
{
    PushStream_(std::vector<T>& t):data(t){};

    std::vector<T>& data;
};

template<class T>
PushStream_<T> PushStream(std::vector<T>& data)
{
    return PushStream_<T>(data);
}

template<class T>
std::istream &operator>>(std::istream &input, PushStream_<T>&& vec)
{
    T temp;

    if(input >> temp)
    {
        vec.data.emplace_back(temp);
    }

    return input;
}

template<class T>
struct InsertStream_
{
    InsertStream_(std::set<T>& t):data(t){};

    std::set<T>& data;
};

template<class T>
InsertStream_<T> InsertStream(std::set<T>& data)
{
    return InsertStream_<T>(data);
}

template<class T>
std::istream &operator>>(std::istream &input, InsertStream_<T>&& vec)
{
    T temp;

    if(input >> temp)
    {
        vec.data.emplace(temp);
    }

    return input;
}

namespace
{
void CheckTrue(bool condition, const std::string& message)
{
    if(!condition)
    {
        std::cout << message << '\n';

        exit(1);
    }
}


void zprintf(std::stringstream& os, const char* s)
{
    while (*s)
    {
        if (*s == '%') throw std::runtime_error("invalid format string: missing arguments");

        os << *s++;
    }
}
}

template<typename T, typename... Args>
void zprintf(std::stringstream& os, const char* s, const T& value, const Args&... args)
{
    while (*s)
    {
        if (*s == '%')
        {
            os << value;

            return zprintf(os, ++s, args...);
        }

        os << *s++;
    }

    throw std::runtime_error("extra arguments provided to zprintf");
}


template<typename... Args>
std::string zformat(const char* s, const Args&... args)
{
    std::stringstream os;

    zprintf(os, s, args...);

    return os.str();
}


class ifstream2 : public std::ifstream
{
public:
	explicit ifstream2( const char* filename, ios_base::openmode mode = ios_base::in ): basic_ifstream(filename, mode)
	{
		if(!*this)
		{
			std::cout << "Could not open " << filename << '\n';
		}
	}

	explicit ifstream2( const std::string& filename, ios_base::openmode mode = ios_base::in ): basic_ifstream(filename, mode)
	{
		if(!*this)
		{
			std::cout << "Could not open " << filename << '\n';
		}
	}
};


class ofstream2 : public std::ofstream
{
public:
	explicit ofstream2( const char* filename, ios_base::openmode mode = ios_base::out ): basic_ofstream(filename, mode)
	{
		if(!*this)
		{
			std::cout << "Could not open " << filename << '\n';
		}
	}

	explicit ofstream2( const std::string& filename, ios_base::openmode mode = ios_base::out ): basic_ofstream(filename, mode)
	{
		if(!*this)
		{
			std::cout << "Could not open " << filename << '\n';
		}
	}
};


static int ConvertChrStr(const std::string& chr)
{
	int result;

	try
	{
		if(chr.length()>3){
			result = std::stoi(chr.substr(3));
		}
		else{	
			result = std::stoi(chr);
		}
	}
	catch(...)
	{
		if(chr == "X")
		{
			result = 23;
		}
		else if(chr == "Y")
		{
			result = 24;
		}
		else if(chr == "XY")
		{
			result = 25;
		}
		else if(chr == "MT")
		{
			result = 26;
		}
		else
		{
			result = 27;
		}
	}

	return result;
}


static std::string ConvertStrChr(const int chr)
{
    std::string result;

    if(chr == 23)
    {
        result = "X";
    }
    else if(chr == 24)
    {
        result = "Y";
    }
    else if(chr == 25)
    {
        result = "XY";
    }
    else if(chr == 26)
    {
        result = "MT";
    }
    else
    {
        result = std::to_string(chr);
    }

    return result;
}


struct SNPclass
{
	std::string rsnum;
	int chr;
	int position;
	std::string ref;
	std::string alt;
	std::string genotypes;

	SNPclass(): rsnum(), chr(), position(), ref(), alt(), genotypes(){}

//	SNPclass(const SNPclass& other): rsnum(other.rsnum), chr(other.chr), position(other.position), ref(other.ref), alt(other.alt), genotypes(other.genotypes){}

//	SNPclass& operator=(const SNPclass& other)
//	{
//		rsnum = other.rsnum;
//		chr = other.chr;
//		position = other.position;
//		ref = other.ref;
//		alt = other.alt;
//		genotypes = other.genotypes;

//		return *this;
//	}

//	SNPclass& operator=(SNPclass&& other)
//	{
//		rsnum = std::move(other.rsnum);
//		chr = std::move(other.chr);
//		position = other.position;
//		ref = std::move(other.ref);
//		alt = std::move(other.alt);
//		genotypes = std::move(other.genotypes);

//		return *this;
//	}

//	SNPclass(SNPclass&& other): rsnum(std::move(other.rsnum)), chr(std::move(other.chr)), position(other.position), ref(std::move(other.ref)), alt(std::move(other.alt)), genotypes(std::move(other.genotypes)){}

//	~SNPclass(){}

	std::pair<int,int> GetPosition() const{return std::make_pair(chr,position);}

	size_t size() const {return genotypes.size();}

    void resize(int num_elements) {genotypes.resize(num_elements);}

	std::string::iterator begin(){return genotypes.begin();}

	std::string::const_iterator cbegin() const {return genotypes.cbegin();}

	std::string::iterator end(){return genotypes.end();}

        std::string::const_iterator cend() const {return genotypes.cend();}

    std::string GetSNPinfo() const
    {
        return rsnum + "\t" + std::to_string(chr) + "\t" + std::to_string(position) + "\t" + ref + "\t" + alt;
    }

    std::string GetSNPinfoeffectfirst() const
    {
        return rsnum + "\t" + std::to_string(chr) + "\t" + std::to_string(position) + "\t" + alt + "\t" + ref;
    }

	double GetEaf() const
	{
		int sum = 0;
		int ntested = 0;

		forc(i, genotypes)
		{
			if(genotypes[i]!='0')
			{
				sum += genotypes[i]-'1';
				++ntested;
			}
		}

		if(ntested > 0)
		{
			return 0.5*double(sum)/ntested;
		}
		else
		{
			return -1;
		}
	}

    double GetCallRate() const
    {
        return double(size() - count(genotypes.begin(), genotypes.end(), '0'))/size();
    }

    double GetHWE() const
    {
        int counts[4] = {0,0,0,0};

        for(auto& v: genotypes)
        {
            ++counts[v-48];
        }

        int countA;
        int countB;
        double statNum = 0.0;
        double statDenom = 0.0;
        double stat;
        double h;
        double v;


        int nSubj = counts[1] + counts[2] + counts[3];

        if(nSubj > 1)
        {
            h = (double(counts[1])*counts[3] - 0.25*counts[2]*(counts[2] - 1))/(nSubj-1);

            countA = 2*counts[1] + counts[2];
            countB = 2*counts[3] + counts[2];

            v = double(countA)*(countA - 1) * double(countB)*(countB -1)/(8.0*Sqr(nSubj-1)*(2.0*nSubj-3));

            statNum   += h;
            statDenom += v;
        }

        if(statDenom > 1.0e-7)
        {
            stat = (statNum * statNum) / statDenom;
        }
        else
        {
            stat = -1;
        }

        return stat;
    }

	bool operator<(const SNPclass& other) const
	{
		return std::make_pair(chr, position) < std::make_pair(chr, other.position);
	}

	bool operator==(const std::pair<int,int>& pos) const
	{
		return GetPosition() == pos;
	}

	char& operator[](int index){return genotypes[index];}

	const char& operator[](int index) const {return genotypes[index];} 
};


struct Imputeclass
{
    double eaf;
    double r2;
    int position;
    bool extra_field = false;
    std::string rsnum;
    std::string chromosome;
    std::string ref;
    std::string alt;
    std::vector<double> genotypes;

    Imputeclass(){}

//    Imputeclass(const Imputeclass& other): rsnum(other.rsnum), chromosome(other.chromosome), position(other.position), ref(other.ref), alt(other.alt), genotypes(other.genotypes), eaf(other.eaf), r2(other.r2){}

//    Imputeclass& operator=(const Imputeclass& other)
//    {
//        rsnum = other.rsnum;
//        chromosome = other.chromosome;
//        position = other.position;
//        ref = other.ref;
//        alt = other.alt;
//        genotypes = other.genotypes;
//        eaf = other.eaf;
//        r2 = other.r2;

//        return *this;
//    }

//    Imputeclass& operator=(Imputeclass&& other)
//    {
//        rsnum = std::move(other.rsnum);
//        chromosome = std::move(other.chromosome);
//        position = other.position;
//        ref = std::move(other.ref);
//        alt = std::move(other.alt);
//        genotypes = std::move(other.genotypes);
//        eaf = other.eaf;
//        r2 = other.r2;

//        return *this;
//    }

//    Imputeclass(Imputeclass&& other): rsnum(std::move(other.rsnum)), chromosome(std::move(other.chromosome)), position(other.position), ref(std::move(other.ref)), alt(std::move(other.alt)), genotypes(std::move(other.genotypes)), eaf(other.eaf), r2(other.r2){}

//    ~Imputeclass(){}

    size_t size() const {return genotypes.size();}

    double GetCallRate() const
    {
        return double(size() - mcount_if(genotypes, a<-0.5))/size();
    }

    double& operator[](int index){return genotypes[index];}

    const double& operator[](int index) const {return genotypes[index];}
};


template<class T>
auto Eaf(const T& t, int)->decltype(t.eaf, double())
{
    return t.eaf;
}

template<class T>
auto Eaf(const T& t, char)->decltype(t.GetEaf(), double())
{
    return t.GetEaf();
}

template<class T>
auto Eaf(const T& t, long)->decltype(double())
{
    return 0.5;
}

template<class T>
auto Eaf(const T& t)->decltype(Eaf(t, '0'), double())
{
    return Eaf(t, '0');
}


template<class T, class U>
int MatchAlleles(const T& lhs, const U& rhs, double ambinclusion = -1.0, double ambfactor = 1.0e-5)
{
    const static std::map<std::string, int> conv{{"A", 1}, {"T",2}, {"C",4}, {"G",8}};

    try
    {
        int lhs1 = conv.at(lhs.ref);

        int lhs2 = conv.at(lhs.alt);

        int rhs1 = conv.at(rhs.ref);

        int rhs2 = conv.at(rhs.alt);

        int value1 = lhs1 | lhs2;

        int value2 = rhs1 | rhs2;

        int check = value1 ^ value2;

        if((check == 0 || check == 15) && abs(value1 - value2) < 8)
        {
            if(value1 < 4 || value1 > 11)
            {
                if(abs(Eaf(lhs) - 0.5) < ambinclusion || abs(Eaf(rhs) - 0.5) < ambinclusion)
                {
                    return 0;
                }
                else
                {
                    return 2*((Eaf(lhs) - 0.5)*(Eaf(rhs) - 0.5) > (1 - 2*(lhs.ref == rhs.ref))*ambfactor) - 1;
                }
            }
            else
            {
                return 2*((lhs1<3)==(rhs1<3)) - 1;
            }
        }
        else
        {
            return 2;
        }
    }
    catch(...)
    {
        if(lhs.ref == "N" || rhs.ref == "N") return 2;

        if(lhs.ref == rhs.ref && lhs.alt == rhs.alt)
        {
            return 1;
        }
        else if(lhs.ref == rhs.alt && lhs.alt == rhs.ref)
        {
            return -1;
        }
        else
        {
            return 2;
        }
    }
}


struct Conv
{
	std::map<std::string, int> conv;

	double ambfactor;

	Conv(double factor=0): ambfactor(factor)
	{
		conv["A"] = 1;
		conv["T"] = 2;
		conv["C"] = 4;
		conv["G"] = 8;
	}

	template<class T> 
	int operator()(const T& lhs, const T& rhs, double ambinclusion = -1) const
	{
		try
		{
			int lhs1 = conv.at(lhs.ref);

			int lhs2 = conv.at(lhs.alt);

			int rhs1 = conv.at(rhs.ref);

			int rhs2 = conv.at(rhs.alt);

			int value1 = lhs1 | lhs2;

			int value2 = rhs1 | rhs2;

			int check = value1 ^ value2;

			if((check == 0 || check == 15) && abs(value1 - value2) < 8)
			{
				if(value1 < 4 || value1 > 11)
				{
					if(abs(lhs.eaf - 0.5) < ambinclusion || abs(rhs.eaf - 0.5) < ambinclusion)
					{
						return 0;
					}
					else
					{
						return 2*((lhs.eaf - 0.5)*(rhs.eaf - 0.5) > (1 - 2*(lhs.ref == rhs.ref))*ambfactor) - 1;
					}
				}
				else
				{
					return 2*((lhs1<3)==(rhs1<3)) - 1;
				}
			}
			else
			{
				return 0;
			}
		}
		catch(...)
		{
			if(lhs.ref == "N" || rhs.ref == "N") return 0;
			
			if(lhs.ref == rhs.ref && lhs.alt == rhs.alt)
			{
				return 1;
			}
			else if(lhs.ref == rhs.alt && lhs.alt == rhs.ref)
			{
				return -1;
			}
			else
			{
				return 0;
			}
		}
	}
};


struct SplitString
{
    std::string linestr;
    std::vector<std::string> linesplit;
    std::string delim;
    boost::algorithm::token_compress_mode_type merge_delim;

    SplitString(const std::string& vdelim="\t", boost::algorithm::token_compress_mode_type vmerge_delim=boost::algorithm::token_compress_off): delim(vdelim), merge_delim(vmerge_delim){}

//    char& operator[](int index){return linesplit[index];}

    const std::string& operator[](int index) const {return linesplit[index];}

    std::string str() const {return linestr;}

    size_t size() const {return linesplit.size();}

    std::vector<std::string>::iterator begin(){return linesplit.begin();}

    std::vector<std::string>::const_iterator cbegin() const {return linesplit.cbegin();}

    std::vector<std::string>::iterator end(){return linesplit.end();}

    std::vector<std::string>::const_iterator cend() const {return linesplit.cend();}
};


namespace
{

double Parsedouble(const char *&p) {
    double r = 0.0;
    bool neg = false;

    double convvec[4] = {1,0.1,0.01,0.001};

    while(*p==' ' || *p=='\t') p++;

    if (*p == '-') {
        neg = true;
        ++p;
    }
    while (*p >= '0' && *p <= '9') {
        r = (r*10.0) + (*p - '0');
        ++p;
    }
    if (*p == '.') {
//        double f = 0.0;
        int n = 0;
        ++p;
        while (*p >= '0' && *p <= '9') {
            r = (r*10.0) + (*p - '0');
            ++p;
            ++n;
        }
        if(n<4)
        {
            r *= convvec[n];
        }
        else
        {
            r /= std::pow(10.0, n);
        }
    }
    if (neg) {
        r = -r;
    }
    return r;
}


std::istream& operator>>(std::istream& is, SplitString& v)
{
    if(std::getline(is, v.linestr))
    {
        if(v.merge_delim==boost::algorithm::token_compress_on)
        {
            boost::algorithm::trim(v.linestr);
        }

        boost::algorithm::split(v.linesplit, v.linestr, boost::algorithm::is_any_of(v.delim), v.merge_delim);
    }

    return is;
}


std::ostream& operator<<(std::ostream& os, const SplitString& v)
{
    os << v.linestr;

    return os;
}


std::istream& operator>>(std::istream& is, SNPclass& v)
{
	is >> v.rsnum >> v.chr >> v.position >> v.ref >> v.alt >> v.genotypes;

	return is;
}


std::ostream& operator<<(std::ostream& os, const SNPclass& v)
{
	os << v.rsnum << '\t' << v.chr << '\t' << v.position << '\t' << v.ref << '\t' << v.alt << '\t' << v.genotypes;

	return os;
}


std::istream& operator>>(std::istream& is, Imputeclass& v)
{
    is >> v.chromosome;

    if(v.extra_field)
    {
        std::string tempstr;

        is >> tempstr;
    }

    is >> v.rsnum >> v.position >> v.ref >> v.alt >> std::ws;

    std::string linestr;

    getline(is, linestr);

    double p0, p1, p2, p1sum=0, p2sum=0, sqrsum=0, eafsum = 0;

    int ncount = 0;

    auto iter = linestr.begin();

    auto end_iter = linestr.end();

    v.genotypes.resize(0);

    while(boost::spirit::qi::phrase_parse(iter, end_iter, boost::spirit::qi::double_ >> boost::spirit::qi::double_ >> boost::spirit::qi::double_, boost::spirit::qi::space, p0, p1, p2))
    {
        double w;

        if(isnan(p0)) p0 = 0;

        if(isnan(p1)) p1 = 0;

        if(isnan(p2))   // use for chromosome X SNP
        {
            p2 = p1;
            p1 = 0;
        }

        double den = p0+p1+p2;

        if(den > 0.9)
        {
            p1 /= den;

            p2 /= den;

            w = (2*p2 + p1);

            sqrsum += Sqr(w);

            p1sum += p1;

            p2sum += p2;

            ++ncount;

            eafsum += w;

            v.genotypes.push_back(w);
        }
        else
        {
//            w = -1;

            v.genotypes.push_back(-1);
        }
    }
/*
    v.genotypes.resize(v.nindiv);

    for(auto& w: v.genotypes)
    {
        p0 = Parsedouble(str);

        p1 = Parsedouble(str);

        p2 = Parsedouble(str);

        double den = p0+p1+p2;

        if(den > 0.9)
        {
            p1 /= den;

            p2 /= den;

            w = (2*p2 + p1);

            sqrsum += Sqr(w);

            p1sum += p1;

            p2sum += p2;

            ++ncount;

            eafsum += w;
        }
        else
        {
            w = -1;
        }
    }
*/
    eafsum /= 2*ncount;

    v.eaf = eafsum;

    v.r2 = (sqrsum-4.0*ncount*eafsum*eafsum)/(p1sum + 4*p2sum - 4.0*ncount*eafsum*eafsum);

    if(eafsum == 0 || eafsum==1) v.r2 = -1;

    return is;
}


template<class V, class T>
void FilterInPlace(V& v, const T& include)
{
    int index = 0;

    forc(i, include)
    {
        if(include[i]==1)
        {
            v[index++] = v[i];
        }
    }

    v.resize(index);
}


template<class V, class T>
V FilterVec(const V& v, const T& include)
{
    int num = 0;

    forc(i, include)
    {
        num += (include[i]==1);
    }

    V result;

    result.resize(num);

    int index = 0;

    forc(i, include)
    {
        if(include[i]==1)
        {
            result[index++] = v[i];
        }
    }

    return result;
}


template<class V, class T>
V FilterMat(const V& v, const T& include)
{
    int num = 0;

    forc(i, include)
    {
        num += (include[i]==1);
    }

    V result;

    result.resize(v.rows(), num);

    int index = 0;

    if(v.rows()>0)
    {
        forc(i, include)
        {
            if(include[i]==1)
            {
                result.col(index++) = v.col(i);
            }
        }
    }

    return result;
}


template<class T>
bool ReadSNPclass(std::istream& is, SNPclass& v, const T& include)
{
    std::string tempstr;

    if(is >> v.rsnum >> v.chr >> v.position >> v.ref >> v.alt >> tempstr)
    {
        v.genotypes.resize(0);

        forc(i, include)
        {
            if(include[i]) v.genotypes.push_back(tempstr[i]);
        }
    }

    return static_cast<bool>(is);
}

template<class T>
bool ReadImputeclass(std::istream& is, Imputeclass& imputeclass, const T& include, bool read_position=true)
{
    if(read_position)
    {
        is >> imputeclass.chromosome;

        if(imputeclass.extra_field)
        {
            std::string tempstr;

            is >> tempstr;
        }

        is >> imputeclass.rsnum >> imputeclass.position;
    }

    if(is >> imputeclass.ref >> imputeclass.alt >> std::ws)
    {
        std::string linestr;

        getline(is, linestr);

        double p0, p1, p2, p1sum=0, p2sum=0, sqrsum=0, eafsum = 0;

        int ncount = 0;

        imputeclass.genotypes.resize(0);

        auto iter = linestr.begin();

        auto end_iter = linestr.end();

        for(auto v: include)
        {
            double w;

//            p0 = Parsedouble(str);

//            p1 = Parsedouble(str);

//            p2 = Parsedouble(str);

            boost::spirit::qi::phrase_parse(iter, end_iter, boost::spirit::qi::double_ >> boost::spirit::qi::double_ >> boost::spirit::qi::double_, boost::spirit::qi::space, p0, p1, p2);

            if(isnan(p0)) p0 = 0;

            if(isnan(p1)) p1 = 0;

            if(isnan(p2))   // use for chromosome X SNP
            {
                p2 = p1;
                p1 = 0;
            }

            double den = p0+p1+p2;

            if(v)
            {
                if(den > 0.9)
                {
                    p1 /= den;

                    p2 /= den;

                    w = (2*p2 + p1);

                    sqrsum += Sqr(w);

                    p1sum += p1;

                    p2sum += p2;

                    ++ncount;

                    eafsum += w;

                    imputeclass.genotypes.push_back(w);
                }
                else
                {
                    //            w = -1;

                    imputeclass.genotypes.push_back(-1);
                }
            }
        }

        eafsum /= 2*ncount;

        imputeclass.eaf = eafsum;

        imputeclass.r2 = (sqrsum-4.0*ncount*eafsum*eafsum)/(p1sum + 4*p2sum - 4.0*ncount*eafsum*eafsum);

        if(eafsum == 0 || eafsum==1) imputeclass.r2 = -1;
    }

    return static_cast<bool>(is);
}


Imputeclass ConvertSNPclass(const SNPclass& snpclass)
{
    Imputeclass result;

    result.rsnum = snpclass.rsnum;

    result.chromosome = std::to_string(snpclass.chr);

    result.position = snpclass.position;

    result.ref = snpclass.ref;

    result.alt = snpclass.alt;

    result.eaf = snpclass.GetEaf();

    result.r2 = snpclass.GetCallRate();

    result.genotypes.resize(snpclass.size());

    mtransform(snpclass, result.genotypes, a-49);

    return result;
}


template<class T>
bool ReadImputeclassFromJdata(std::istream& is, Imputeclass& v, const T& include)
{
    SNPclass snpclass;

    ReadSNPclass(is, snpclass, include);

    v = ConvertSNPclass(snpclass);

    return static_cast<bool>(is);
}


std::ostream& PrintSNP(std::ostream& os, const SNPclass& v)
{
	os << v.rsnum << '\t' << v.chr << '\t' << v.position << '\t' << v.ref << '\t' << v.alt;

	return os;
}

template< class SinglePassRange, class Value>
inline auto count(const SinglePassRange& rng, const Value& val)->decltype(std::count(rng.cbegin(), rng.cend(), val))
{
    return std::count(rng.cbegin(), rng.cend(), val);
}

template< class SinglePassRange>
inline void sort(SinglePassRange& rng)
{
    std::sort(rng.begin(), rng.end());
}


std::istream& skipline(std::istream& is)
{
	is.ignore(std::numeric_limits<int>::max(), '\n');

	return is;
}


std::map<std::pair<int,int>,std::pair<int,int>> ReadRegionFile(const std::string& inputfile)
{
	std::map<std::pair<int,int>, std::pair<int,int>> regionmap;

	int chromosome, position1, position2;

	std::ifstream input(inputfile);

	while(input >> chromosome >> position1 >> position2)
	{
		if(regionmap.insert(std::make_pair(std::make_pair(chromosome, position2), std::make_pair(chromosome, position1))).second==false)
		{
			if(std::make_pair(chromosome, position1) < regionmap[std::make_pair(chromosome, position2)])
			{
				regionmap[std::make_pair(chromosome, position2)] = std::make_pair(chromosome, position1);
			}
		}
	}

	forl(i, 2)
	{
		for(auto itor=regionmap.begin();itor!=regionmap.end();)
		{
			auto p = regionmap.upper_bound(itor->first);

			if(p!=regionmap.end() && p->second <= itor->second)
			{
				regionmap.erase(itor++);
			}
			else
			{
				++itor;
			}
		}
	}

	return regionmap;
}


#pragma warning(disable: 4996)

std::string sformat(const std::string &fmt, ...) {
    int size = 100;
    std::string str;
    va_list ap;
    while (1) {
        str.resize(size);
        va_start(ap, fmt);
        int n = vsnprintf((char *)str.c_str(), size, fmt.c_str(), ap);
        va_end(ap);
        if (n > -1 && n < size) {
            str.resize(n);
            return str;
        }
        if (n > -1)
            size = n + 1;
        else
            size *= 2;
    }
    return str;
}

#pragma warning(default: 4996)


int GetNumIndivsFromFile(const std::string& inputfile)
{
    std::ifstream input(inputfile);

    std::string tempstr;

    forl(i,6) input >> tempstr;

    return int(tempstr.size());
}

inline double SumSimd(const __m128d& value)
{
    double temp[2];

    _mm_storeu_pd(temp, value);

    return temp[0] + temp[1];
}

inline long long SumSimdInt(const __m128i& value)
{
    long long temp[2];

    _mm_storeu_si128((__m128i*) temp, value);

    return temp[0] + temp[1];
}

double CalcCorrFromVec(const std::vector<double>& data1, const std::vector<double>& data2)
{
    __m128d x2_mm = _mm_setzero_pd();

    __m128d y2_mm = _mm_setzero_pd();

    __m128d xy_mm = _mm_setzero_pd();

    __m128d xsum_mm = _mm_setzero_pd();

    __m128d ysum_mm = _mm_setzero_pd();

    __m128d missing = _mm_set1_pd(-0.5);

    __m128i numvalid_mm = _mm_setzero_si128();

    for(int i=0;i<2*(data1.size()/2);i+=2)
    {
        auto x = _mm_loadu_pd(&data1[i]);

        auto y = _mm_loadu_pd(&data2[i]);

        auto include = _mm_and_pd(_mm_cmpgt_pd(x,missing), _mm_cmpgt_pd(y,missing));

        x = _mm_and_pd(x, include);

        y = _mm_and_pd(y, include);

        x2_mm = _mm_add_pd(x2_mm, _mm_mul_pd(x, x));

        xsum_mm = _mm_add_pd(xsum_mm, x);

        y2_mm = _mm_add_pd(y2_mm, _mm_mul_pd(y, y));

        ysum_mm = _mm_add_pd(ysum_mm, y);

        xy_mm = _mm_add_pd(xy_mm, _mm_mul_pd(x, y));

        numvalid_mm = _mm_add_epi64(numvalid_mm, (__m128i&) include);

//        if(status[i] > -0.5)
//        {
//            *p1++ = data1[i];
//            *p2++ = data2[i];
//        }
    }

    double x2 = SumSimd(x2_mm);

    double y2= SumSimd(y2_mm);

    double xy= SumSimd(xy_mm);

    double xsum = SumSimd(xsum_mm);

    double ysum = SumSimd(ysum_mm);

    long long numvalid = -SumSimdInt(numvalid_mm);


    for(int i=2*(data1.size()/2);i<data1.size();i++)
    {
        if(data1[i] > -0.5 && data2[i] > -0.5)
        {
            x2 += Sqr(data1[i]);

            xsum += data1[i];

            y2 += Sqr(data2[i]);

            ysum += data2[i];

            xy += data1[i]*data2[i];

            ++numvalid;
        }
    }

    double nrecip = 1.0/numvalid;

    double r2 = copysign(Sqr(xy - nrecip*xsum*ysum)/((x2-nrecip*Sqr(xsum))*(y2-nrecip*Sqr(ysum))), xy - nrecip*xsum*ysum);

    return r2;
}


double CalcCorr(const Imputeclass& data1, const Imputeclass& data2)
{
    return CalcCorrFromVec(data1.genotypes, data2.genotypes);
}

inline int SumSimdChar(const __m128i& value)
{
    char temp[16] __attribute__((aligned(16)));

    _mm_store_si128((__m128i*) temp, value);

    int sum = 0;

    forl(i,16) sum += int(temp[i]);

    return sum;
}

#ifdef __SSSE3__
double CalcCorr(const SNPclass& data1, const SNPclass& data2)
{
    __m128i corr = _mm_set1_epi8('2');

    __m128i missing = _mm_set1_epi8(-2);

    double x2 = 0;

    double y2 = 0;

    double xy = 0;

    double xsum = 0;

    double ysum = 0;

    int nummissing = 0;

    const int maxnumber = 2032;

    for(int j=0;j<=maxnumber*(data1.size()/maxnumber);j+=maxnumber)
    {
        __m128i x2_mm = _mm_setzero_si128();

        __m128i y2_mm = _mm_setzero_si128();

        __m128i xy_mm = _mm_setzero_si128();

        __m128i xsum_mm = _mm_setzero_si128();

        __m128i ysum_mm = _mm_setzero_si128();

        __m128i nummissing_mm = _mm_setzero_si128();

        for(int i=j;i<std::min(int(16*(data1.size()/16)), j+maxnumber);i+=16)
        {
            auto x = _mm_loadu_si128((const __m128i*)&data1[i]);

            auto y = _mm_loadu_si128((const __m128i*)&data2[i]);

            x = _mm_sub_epi8(x, corr);

            y = _mm_sub_epi8(y, corr);

            auto exclude = _mm_or_si128(_mm_cmpeq_epi8(x,missing), _mm_cmpeq_epi8(y,missing));

            x = _mm_andnot_si128(exclude, x);

            x2_mm = _mm_add_epi8(x2_mm, _mm_abs_epi8(x));

            xsum_mm = _mm_add_epi8(xsum_mm, x);

            y = _mm_andnot_si128(exclude, y);

            y2_mm = _mm_add_epi8(y2_mm, _mm_abs_epi8(y));

            ysum_mm = _mm_add_epi8(ysum_mm, y);

            nummissing_mm = _mm_add_epi8(nummissing_mm, exclude);

            xy_mm = _mm_add_epi8(xy_mm, _mm_sign_epi8(x, y));
        }

        x2 += SumSimdChar(x2_mm);

        y2 += SumSimdChar(y2_mm);

        xy += SumSimdChar(xy_mm);

        xsum += SumSimdChar(xsum_mm);

        ysum += SumSimdChar(ysum_mm);

        nummissing -= SumSimdChar(nummissing_mm);
    }

    for(int i=16*(data1.size()/16);i<data1.size();i++)
    {
        if(data1[i]!='0' && data2[i]!='0')
        {
            double v1 = data1[i]-'2';

            double v2 = data2[i]-'2';

            x2 += Sqr(v1);

            xsum += v1;

            y2 += Sqr(v2);

            ysum += v2;

            xy += v1*v2;
        }
        else
        {
            ++nummissing;
        }
    }

    double nrecip = 1.0/(data1.size() - nummissing);

    double r2 = copysign(Sqr(xy - nrecip*xsum*ysum)/((x2-nrecip*Sqr(xsum))*(y2-nrecip*Sqr(ysum))), xy - nrecip*xsum*ysum);

    if(abs(x2-nrecip*Sqr(xsum)) < 1.0e-10 || abs(y2-nrecip*Sqr(ysum)) < 1.0e-10) r2 = 0;

    if(std::isinf(r2))
    {
        std::cout << "infinite\n";

        std::cout << abs(y2-nrecip*Sqr(ysum)) << '\n';
    }

    return r2;
}
#endif
}
