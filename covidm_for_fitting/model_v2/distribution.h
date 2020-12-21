// distribution.h

#ifndef DISTRIBUTION_H
#define DISTRIBUTION_H

#include "randomizer.h"
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_cdf.h>
#include <algorithm>
#include <string>
#include <iostream>

class Distribution
{
public:
    enum Type { Cauchy, LogNormal, Normal, Uniform, Beta, Exponential, Gamma, Rounded };
    static const std::string Types[];

    Distribution(Type t, std::vector<double> parameters, std::vector<double> trunc = {}, std::vector<double> shift = {}, std::vector<double> init = {});
    Distribution(std::string code);
    Distribution(std::vector<std::string> subcode, std::string code = "");

    double Random(Randomizer& R);
    double RandomInit(Randomizer& R);
    double LogProbability(double x);
    double LowerBound();
    double UpperBound();

    friend std::ostream& operator << (std::ostream& out, const Distribution& d)
    {
        out << Types[d.type] << "(" << d.a << ", " << d.b << ")";
        if (d.t_lp_adj != 0)
            out << " T [" << d.t0 << ", " << d.t1 << "] adj " << d.t_lp_adj;
        if (d.s0 != 0 || d.s1 != 1)
            out << " * " << d.s1 << " + " << d.s0;
        return out;
    }

private:
    static constexpr double ShoulderArea = 1e-20;

    void Init(Type t, std::vector<double> parameters, std::vector<double> trunc = {}, std::vector<double> shift = {}, std::vector<double> init = {});
    void SetTrunc(double xmin, double xmax);
    void SetFromCode(std::vector<std::string> subcode, std::string code);
    static std::vector<std::string> Unserialize(std::string s);
    static double RoundedUniformCDF(double min, double max, double shoulder, double x);

    Type type;
    double a, b;
    double t0, t1, t_lp_adj;
    double s0, s1;
    double i0, i1, i_amount;
};

#endif