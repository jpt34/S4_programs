#ifndef IMPUTECLASSDETAILED_H
#define IMPUTECLASSDETAILED_H

#include <jlibrary.h>

using namespace std;

class ImputeclassDetailed : public Imputeclass
{
public:
    vector<array<double, 3>> probs;

    ImputeclassDetailed();

    tuple<double, double, int, int, bool> GetR2statistics(const vector<int> &include, int start, int finish);

    map<string, tuple<double, double, int, int, bool>> GetR2statisticsGroup(const vector<string> &group);

    void Filter(const vector<int>& include, Imputeclass &imputeclass);
};


std::istream& operator>>(std::istream& is, ImputeclassDetailed& v);

#endif // IMPUTECLASSDETAILED_H
