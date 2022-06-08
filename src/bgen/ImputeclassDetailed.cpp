#include "ImputeclassDetailed.h"

ImputeclassDetailed::ImputeclassDetailed()
{
}

tuple<double, double, int, int, bool> ImputeclassDetailed::GetR2statistics(const vector<int> &include, int start, int finish)
{
    double p1sum=0, p2sum=0, sqrsum=0, eafsum = 0;

    int ncount = 0, ngenotyped = 0;

    forv(i, start, finish)
    {
        if(include[i] && genotypes[i] > -0.5)
        {
            sqrsum += Sqr(genotypes[i]);

            p1sum += probs[i][1];

            p2sum += probs[i][2];

            ++ncount;

            eafsum += genotypes[i];

            ngenotyped += (probs[i][0]==1 || probs[i][1]==1 || probs[i][2]==1);
        }
    }

    eafsum /= 2*ncount;

    double r2value = (sqrsum-4.0*ncount*eafsum*eafsum)/(p1sum + 4*p2sum - 4.0*ncount*eafsum*eafsum);

    if(eafsum == 0 || eafsum==1) r2value = -1;

    bool is_genotyped = r2value > 0.99 && double(ngenotyped)/ncount > 0.94;

    return make_tuple(r2value, eafsum, ngenotyped, ncount, is_genotyped);
}


map<string, tuple<double, double, int, int, bool>> ImputeclassDetailed::GetR2statisticsGroup(const vector<string> &group)
{
    map<string, double> p1sum, p2sum, sqrsum, eafsum;

    map<string, int> ncount, ngenotyped;

    map<string, tuple<double, double, int, int, bool>> result;

    forc(i, group)
    {
        string v = group[i];

        if(group[i]!="-1" && group[i]!="-99" && group[i]!="NA" && genotypes[i] > -0.5)
        {
            sqrsum[v] += Sqr(genotypes[i]);

            p1sum[v] += probs[i][1];

            p2sum[v] += probs[i][2];

            ++ncount[v];

            eafsum[v] += genotypes[i];

            ngenotyped[v] += (probs[i][0]==1 || probs[i][1]==1 || probs[i][2]==1);

            if(sqrsum[v] > p1sum[v] + 4*p2sum[v])
            {
                cout << genotypes[i] << ' ' << sqrsum[v] << ' ' << p1sum[v] + 4*p2sum[v] << endl;
            }
        }
    }

    for(auto& w:p1sum)
    {
        auto v = w.first;

        if(v=="All") continue;

        sqrsum["All"] += sqrsum[v];

        p1sum["All"] += p1sum[v];

        p2sum["All"] += p2sum[v];

        ncount["All"] += ncount[v];

        eafsum["All"] += eafsum[v];

        ngenotyped["All"] += ngenotyped[v];
    }

    for(auto& w:p1sum)
    {
        auto v = w.first;

        eafsum[v] /= 2*ncount[v];

        double r2value = (sqrsum[v]-4.0*ncount[v]*eafsum[v]*eafsum[v])/(p1sum[v] + 4*p2sum[v] - 4.0*ncount[v]*eafsum[v]*eafsum[v]);

        if(eafsum[v] == 0 || eafsum[v]==1) r2value = -1;

        bool is_genotyped = r2value > 0.99 && double(ngenotyped[v])/ncount[v] > 0.94;

        result[v] = make_tuple(r2value, eafsum[v], ngenotyped[v], ncount[v], is_genotyped);
    }

    return result;
}

void ImputeclassDetailed::Filter(const vector<int>& include, Imputeclass& imputeclass)
{
    imputeclass.genotypes.resize(0);

    tie(imputeclass.rsnum,imputeclass.chromosome,imputeclass.position,imputeclass.ref,imputeclass.alt) = {rsnum,chromosome,position,ref,alt};

    double p1sum=0, p2sum=0, sqrsum=0, eafsum = 0;

    int ncount = 0;

    forc(i, probs)
    {
        if(include[i])
        {
            double w = 2*probs[i][2]+probs[i][1];

            sqrsum += Sqr(w);

            p1sum += probs[i][1];

            p2sum += probs[i][2];

            ++ncount;

            eafsum += w;

            imputeclass.genotypes.push_back(w);
        }
    }

    eafsum /= 2*ncount;

    imputeclass.eaf = eafsum;

    imputeclass.r2 = (sqrsum-4.0*ncount*eafsum*eafsum)/(p1sum + 4*p2sum - 4.0*ncount*eafsum*eafsum);

    if(eafsum == 0 || eafsum==1) imputeclass.r2 = -1;
}


std::istream& operator>>(std::istream& is, ImputeclassDetailed& v)
{
    is >> v.chromosome >> v.rsnum >> v.position >> v.ref >> v.alt >> std::ws;

    std::string linestr;

    getline(is, linestr);

    double p0, p1, p2, p1sum=0, p2sum=0, sqrsum=0, eafsum = 0;

    int ncount = 0;

    auto iter = linestr.begin();

    auto end_iter = linestr.end();

    v.genotypes.resize(0);

    v.probs.resize(0);

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
//            p1 /= den;

//            p2 /= den;

            w = (2*p2 + p1);

            sqrsum += Sqr(w);

            p1sum += p1;

            p2sum += p2;

            ++ncount;

            eafsum += w;

            v.genotypes.push_back(w);

            v.probs.push_back({p0,p1,p2});
        }
        else
        {
//            w = -1;

            v.genotypes.push_back(-1);

            v.probs.push_back({0,0,0});
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
