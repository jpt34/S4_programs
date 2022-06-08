#ifndef MLOGIT_H
#define MLOGIT_H

#include "routines.h"
#include <Eigen/Dense>
#include <VectorclassNew3/vectorclass.h>
#include <VectorclassNew3/vectormath_exp.h>
#include <jlibrary.h>

using namespace std;
using namespace Eigen;

// This calculates logistic regression. It was orignally designed for multinomial logit but also can be used for logisitic regression.


class Mlogit
{
public:
    Mlogit(const PhenotypeData& covdata);
    void LogitStudyStats(const MatrixXd &cov, VectorXd& x, MatrixXd& hessian) const;
    int nstudies;
    int ngroups;
    PhenotypeData phenotypedata;
    VectorXd xstart;
    VectorXd xstartnull;
    double liknull;
    double CalcFullLik(const MatrixXd &cov, VectorXd &param) const;
    void LogitStats(const MatrixXd &cov, double &lik, VectorXd &x, MatrixXd &hessian) const;
    void LogitStats_sse(const MatrixXd &cov, const MatrixXd &cov_sse, double &lrt, VectorXd &x, MatrixXd &hessian) const;
    template<class T>
    void MaximiseData(const VectorXi &status, const MatrixXd &cov, const VectorXi &study, double &lik, VectorXd &x, MatrixXd &hessian) const;
};


// This is structure used to calculate logistic regression estimates

struct MlogitFunctorStudy
{
    const VectorXi& status;
    const MatrixXd& cov;
    const VectorXi& study;
    int nstudies;
    int ngroups;
    int nparam1;

    MlogitFunctorStudy(const VectorXi& statusvalue, const MatrixXd& covvalue, const VectorXi& studyvalue, int nstudiesvalue, int nparam1value): status(statusvalue),  cov(covvalue), study(studyvalue), nstudies(nstudiesvalue), nparam1(nparam1value)
    {
        ngroups = statusvalue.maxCoeff();
    }

    double CalcFullLik(const VectorXd& param) const
    {
        vector<double> p(ngroups);

        double lik = 0;

        for (unsigned i=0;i<status.size();i++)
        {
            const int offset_study = ngroups*nparam1 + ngroups*study[i];

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

            if(status[i]==0)
            {
                lik -= log(psum);
            }
            else
            {
                lik += log(p[status[i]-1]);
            }
        }

        return lik;
    }

    void CalcAll(const VectorXd& param, double& lik, VectorXd& value, MatrixXd& mat)
    {
        vector<double> p(ngroups);

        lik = 0;

        value.setZero(param.size());

        mat.setZero(param.size(), param.size());

        forc(i, status)
        {
            const int offset_study = ngroups*nparam1 + ngroups*study[i];

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

                if(paramsum > 100)
                {
                    paramsum = 100;
                }
                else if(paramsum < -100)
                {
                    paramsum = -100;
                }


                p[k] = exp(paramsum);

                psum += p[k];
            }

            forl(k, ngroups)
            {
                p[k] /= psum;
            }

            if(status[i]==0)
            {
                lik -= log(psum);
            }
            else
            {
                lik += log(p[status[i]-1]);
            }

            forl(k, ngroups)
            {
                for (int j=0;j<nparam1;j++)
                {
                    value[j+k*nparam1] -= cov(j,i)*((status[i]==k+1) - p[k]);
                }
            }

            forl(k, ngroups)
            {
                value[offset_study+k] -= ((status[i]==k+1) - p[k]);
            }

            forl(m, ngroups)
            {
                forl(n, m+1)
                {
                    double var = (n==m?p[m]*(1-p[m]):-p[m]*p[n]);

                    forl(j, nparam1)
                    {
                        forl(k, j+1)
                        {
                            mat(n*nparam1+k, m*nparam1+j) += cov(j,i)*cov(k,i)*var;
                        }
                    }

                    forl(k, nparam1)
                    {
                        mat(n*nparam1+k, offset_study+m) += cov(k,i)*var;
                    }

                    mat(offset_study+n, offset_study+m) += var;
                }
            }
        }

        forl(m, ngroups)
        {
            forl(n, m+1)
            {
                forl(j, nparam1)
                {
                    forl(k, j)
                    {
                        mat(n*nparam1+j, m*nparam1+k) = mat(n*nparam1+k, m*nparam1+j);
                    }
                }

                forl(j, nstudies)
                {
                    forl(k, nparam1)
                    {
                        mat(m*nparam1+k, ngroups*nparam1 + ngroups*j+n) = mat(n*nparam1+k, ngroups*nparam1 + ngroups*j+m);
                    }
                }
            }
        }

        forl(i, param.size())
        {
            forl(j, i)
            {
                mat(i,j) = mat(j,i);
            }
        }

        forl(i, param.size())
        {
            if(mat(i,i)==0) mat(i,i) = 1.0;   // to make sure hessian is invertible
        }
    }
};


// This is a structure used to calculate logistic regression estimates using sse instructions

struct LogitFunctorStudy_sse2
{
    const VectorXi& status;
    const MatrixXd& cov;
    const VectorXi& study;
    int nstudies;
    int nparam1;

    LogitFunctorStudy_sse2(const VectorXi& statusvalue, const MatrixXd& covvalue, const VectorXi& studyvalue, int nstudiesvalue, int nparam1value): status(statusvalue),  cov(covvalue), study(studyvalue), nstudies(nstudiesvalue), nparam1(nparam1value)
    {
    }

    int values() const
    {
        return (1+nstudies);
    }

    // This is used to update results for an individual that can't be calculated by the sse instructions as the total number of individuals is odd

    void UpdateIndividual(const VectorXd& param, int index, double& lik, VectorXd& value, MatrixXd& mat) const
    {
        double p;

        const int index2 = index*2;

        const int offset_study = nparam1 + study[index2];

        double psum = 1;

        double paramsum = 0;

        for (int j=0;j<nparam1;j++)
        {
            paramsum += cov(j,index)*param[j];
        }

        paramsum += param[offset_study];

        if(paramsum > 100)
        {
            paramsum = 100;
        }
        else if(paramsum < -100)
        {
            paramsum = -100;
        }


        p = exp(paramsum);

        psum += p;

        p /= psum;

        if(status[index2]==0)
        {
            lik -= log(psum);
        }
        else
        {
            lik += log(p);
        }

        for (int j=0;j<nparam1;j++)
        {
            value[j] -= cov(j,index)*(status[index2] - p);
        }

        value[offset_study] -= (status[index2] - p);

        double var = p*(1-p);

        forl(j, nparam1)
        {
            forl(k, j+1)
            {
                mat(k, j) += cov(j,index)*cov(k,index)*var;
            }
        }

        forl(k, nparam1)
        {
            mat(k, offset_study) += cov(k,index)*var;
        }

        mat(offset_study, offset_study) += var;

        forl(j, nparam1)
        {
            forl(k, j)
            {
                mat(j, k) = mat(k, j);
            }
        }

        forl(j, nstudies)
        {
            forl(k, nparam1)
            {
                mat(k, nparam1 + j) = mat(k, nparam1 + j);
            }
        }

        forl(i, param.size())
        {
            forl(j, i)
            {
                mat(i,j) = mat(j,i);
            }
        }
    }

    // This calculates the likelihood, derivative and hessian together
    // It is quicker to do this than to calculate them individually
    // The study vector is used separately so that many hessian values that will be 0 aren't calculated

    void CalcAll(const VectorXd& param_orig, double& lik_orig, VectorXd& value_orig, MatrixXd& mat_orig) const
    {
        Vec2d p;

        Vec2d lik = 0.0;

        vector<Vec2d> param(param_orig.size());

        forc(i, param_orig)
        {
            param[i] = Vec2d(param_orig[i]);
        }

        value_orig.resize(param_orig.size());

        vector<Vec2d> value(value_orig.size(), Vec2d(0.0));

        mat_orig.resize(param.size(), param.size());

        MatrixXd mat(param.size()*2, param.size());

        Map<Matrix<Vec2d,-1,-1>> mat_view((Vec2d*)&mat(0,0), param.size(), param.size());

        Map<Matrix<Vec2d,-1,-1>> cov_view((cov.rows()>0?(Vec2d*)&cov(0,0):(Vec2d*)&mat(0,0)), cov.rows()/2, status.size()/2);    // make sure last calculation doesn't use overlapped data

        mat.setZero(param.size()*2, param.size());

        Vec2db mask1(-1,0);

        Vec2db mask2(0,-1);

        Vec2d one(1);

        for (int i=0;i<cov_view.cols();i++)
        {
            const int offset_study1 = nparam1 + study[i*2];

            const int offset_study2 = nparam1 + study[i*2+1];

            Vec2d psum = 1;

            Vec2d status_vec(status[i*2], status[i*2+1]);

            Vec2d paramsum = 0;

            for (int j=0;j<nparam1;j++)
            {
                paramsum += cov_view(j,i)*param[j];
            }

            paramsum += (mask1 & param[offset_study1]);

            paramsum += (mask2 & param[offset_study2]);

            p = one/(one+exp(-paramsum));

            lik += log(select(status_vec==0,1-p,p));


            for (int j=0;j<nparam1;j++)
            {
                value[j] -= cov_view(j,i)*(status_vec - p);
            }

            value[offset_study1] -= mask1&(status_vec - p);

            value[offset_study2] -= mask2&(status_vec - p);

            Vec2d var = p*(1-p);

            forl(j, nparam1)
            {
                forl(k, j+1)
                {
                    mat_view(k, j) += cov_view(j,i)*cov_view(k,i)*var;
                }
            }

            forl(k, nparam1)
            {
                mat_view(k, offset_study1) += mask1&cov_view(k,i)*var;

                mat_view(k, offset_study2) += mask2&cov_view(k,i)*var;
            }

            mat_view(offset_study1, offset_study1) += mask1&var;

            mat_view(offset_study2, offset_study2) += mask2&var;
        }

        forl(j, nparam1)
        {
            forl(k, j)
            {
                mat_view(j, k) = mat_view(k, j);
            }
        }

        forl(j, nstudies)
        {
            forl(k, nparam1)
            {
                mat_view(k, nparam1 + j) = mat_view(k, nparam1 + j);
            }
        }


        forl(i, param.size())
        {
            forl(j, i+1)
            {
                mat_orig(i,j) = mat_orig(j,i) = horizontal_add(mat_view(j,i));
            }
        }

        forl(i, param.size())
        {
            if(mat_orig(i,i)==0) mat_orig(i,i) = 1.0;   // to make sure hessian is invertible

            value_orig[i] = horizontal_add(value[i]);
        }

        lik_orig = horizontal_add(lik);

        if(status.size()%2)
        {
            UpdateIndividual(param_orig, cov.cols()-1, lik_orig, value_orig, mat_orig);
        }

        forl(i, param.size())
        {
            if(mat(i,i)==0) mat(i,i) = 1.0;   // to make sure hessian is invertible
        }
    }
};


#endif // MLOGIT_H
