#include "Mlogit.h"

#include <MachineTimer.h>

// This initialises the logit model
// The null values are calculated so that they can be used as a starting point. The null likelihood is also calculated

Mlogit::Mlogit(const PhenotypeData& covdata)
{
    phenotypedata = covdata;

    nstudies = phenotypedata.study.maxCoeff()+1;

    ngroups = phenotypedata.status.maxCoeff();

    map<int, int> groups_chosen;

    forc(i, phenotypedata.status)
    {
        if(phenotypedata.status[i]!=-1)
        {
            groups_chosen[phenotypedata.status[i]]++;
        }
    }

    forl(i, ngroups+1)
    {
        if(groups_chosen[i]==0)
        {
            cout << "No valid outcomes for outcome " << i << " with " << ngroups+1 << " groups overall\n";

            exit(1);
        }
    }

    if(ngroups > 30)
    {
        cout << "Over 30 groups specified: might be an error in the phenotype file\n";

        exit(1);
    }

    xstartnull = VectorXd::Zero((phenotypedata.cov.rows() + nstudies)*ngroups);

    MatrixXd hessianorig;

    double lik;

    MatrixXd cov_sse;

    ConvertToSseMatrix(phenotypedata.cov, cov_sse);

    liknull = 0;

    LogitStats_sse(phenotypedata.cov, cov_sse, lik, xstartnull, hessianorig);

    LogitStats_sse(phenotypedata.cov, cov_sse, lik, xstartnull, hessianorig);

    liknull = lik;

    xstart.resize(xstartnull.size()+ngroups);

    int j=0;

    if(phenotypedata.cov.rows()==0)
    {
        forl(i, ngroups)
        {
            xstart[i]=0;
        }

        forc(i, xstartnull)
        {
            xstart[i+ngroups] = xstartnull[i];
        }
    }
    else
    {
        forc(i, xstartnull)
        {
            if(i%phenotypedata.cov.rows()==0 && i/phenotypedata.cov.rows() < ngroups)
            {
                xstart[j++]=0;
            }

            xstart[j++] = xstartnull[i];
        }
    }
}


// This calculates the logistic regrssion for the prs. It doesn't use sse instructions.

void Mlogit::LogitStats(const MatrixXd& cov, double& lik, VectorXd& x, MatrixXd& hessian) const
{
    if(cov.rows() > phenotypedata.cov.rows())
    {
        x = xstart;
    }
    else
    {
        x = xstartnull;
    }

    using  FunctorType = MlogitFunctorStudy;

    VectorXi status, study;

    FunctorType functor(phenotypedata.status, cov, phenotypedata.study, nstudies, cov.rows());

    MatrixXd invhessian;

    VectorXd value(x.size());

    double prev_lik = -1000000;

    lik =-200000;

    VectorXd change;

    while(abs(lik - prev_lik) > 1.0e-5)
    {
        prev_lik = lik;

        functor.CalcAll(x, lik, value, hessian);

        while(prev_lik > lik+1e-9 || !isfinite(lik))
        {
            change = 0.5*change;

            x += change;

            functor.CalcAll(x, lik, value, hessian);
        }

        invhessian = hessian.inverse();

        change = invhessian*value;

        x -= change;
    }

    hessian = invhessian;
}


// This calculates the logistic regression for the null model using sse instructions

void Mlogit::LogitStats_sse(const MatrixXd& cov, const MatrixXd& cov_sse, double& lrt, VectorXd& x, MatrixXd& hessian) const
{
    if(cov.rows() > phenotypedata.cov.rows())
    {
        x = xstart;
    }
    else
    {
        x = xstartnull;
    }

    VectorXi status, study;

    MatrixXd cov_analysis, cov_analysis_null;

    double lik;

    auto Func = [this](const VectorXi &status, const MatrixXd &cov, const VectorXi &study, double &lik, VectorXd &x, MatrixXd &hessian)
    {
        MaximiseData<LogitFunctorStudy_sse2>(status, cov, study, lik, x, hessian);
    };


    if(cov.rows() > phenotypedata.cov.rows())
    {
        Func(phenotypedata.status, cov_sse, phenotypedata.study, lik, x, hessian);

        lrt = 2*(lik - liknull);
    }
    else
    {
        Func(phenotypedata.status, cov_sse, phenotypedata.study, lrt, x, hessian);
    }
}


// This maximises the likelihood for the logistic regression model.

template<class T>
void Mlogit::MaximiseData(const VectorXi& status, const MatrixXd& cov, const VectorXi& study, double& lik, VectorXd& x, MatrixXd& hessian) const
{
    T functor(status, cov, study, nstudies, cov.rows()/2);

    MatrixXd invhessian;

    VectorXd value(x.size());

    double prev_lik = -1.0e20;

    lik = -1.0e15;

    VectorXd change;

    int max_runs = 0;

    VectorXd prev_x;

    bool reject_value = true;

    int ntries = 0;

    // This updates the likelihood using the newton step except when this decreases the likelihood in which case a reduced step is taken

    while((abs(lik - prev_lik) > 1.0e-5 && max_runs++ < 40000) || (reject_value && ntries<3))
    {
        reject_value = false;

        prev_lik = lik;

        functor.CalcAll(x, lik, value, hessian);

        for(int runnum=2;lik+1.0e-12 < prev_lik || !isfinite(lik);runnum++)
        {
            reject_value = true;

            ntries++;

            change = pow(0.7,runnum)*change;

            x = prev_x - change;

            functor.CalcAll(x, lik, value, hessian);
        }

        invhessian = hessian.inverse();

        change = invhessian*value;

        if(!isfinite(change.sum()))
        {
            hessian += 1.0e-7*MatrixXd::Identity(hessian.rows(), hessian.cols());

            invhessian = hessian.inverse();

            change = invhessian*value;
        }

        prev_x = x;

        x -= change;
    }

    hessian = invhessian;
}

// This calculates the logistic likelihood for the model parameters

double Mlogit::CalcFullLik(const MatrixXd& cov, VectorXd& param) const
{
    int ngroups = phenotypedata.status.maxCoeff();

    vector<double> p(ngroups);

    int nparam1 = cov.rows();

    double lik = 0;

    for (unsigned i=0;i<phenotypedata.status.size();i++)
    {
        if (phenotypedata.status[i]==-1) continue;

        const int offset_study = ngroups*nparam1 + ngroups*phenotypedata.study[i];

        double psum = 1;

        forl(k, ngroups)
        {
            int offset_k = k*nparam1;

            double paramsum = 0;

            for (int j=0;j<nparam1;j++)
            {
                paramsum += cov(j,i)*param[j+offset_k];
            }

            paramsum += param[offset_study+k];

            p[k] = exp(paramsum);

            psum += p[k];
        }

        forl(k, ngroups)
        {
            p[k] /= psum;
        }

        if(phenotypedata.status[i]==0)
        {
            lik -= log(psum);
        }
        else
        {
            lik += log(p[phenotypedata.status[i]-1]);
        }
    }

    return lik;
}
