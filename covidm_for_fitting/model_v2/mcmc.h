// mcmc.h

#ifndef MCMC_H
#define MCMC_H

#include <vector>
#include <string>
#include "parameters.h"
#include "randomizer.h"
#include "distribution.h"

// Functor for storing the results of optimization.
struct MCMCReporter
{
public:
    MCMCReporter(unsigned int iterations, unsigned int chains, std::vector<std::string>& param_names);

    void operator()(int TR, double LP, double CH, double LL, std::vector<double>& TH);

    unsigned int n_chains, n_samp;
    std::vector<int> trial, chain;
    std::vector<double> lp, ll;
    std::vector<std::vector<double>> theta;
    std::vector<std::string> pnames;
};

// Functor for calculating the likelihood of a model run.
struct Likelihood
{
public:
    Likelihood(Parameters& bp, unsigned int m_seed);
    Likelihood(vector<Parameters>& bp, vector<vector<unsigned int>>& xs, vector<double>& fp, unsigned int m_seed);

    double operator()(const std::vector<double>& theta);

private:
    Parameters base_parameters;
    unsigned int model_seed;

    bool multi;
    vector<Parameters> multi_base_parameters;
    vector<vector<unsigned int>> x_source;
    vector<double> fixed_params;
};

void DEMCMC_Priors(Randomizer& R, Likelihood& likelihood, MCMCReporter& report,
    int burn_in, int iterations, int n_chains, std::vector<Distribution>& priors,
    bool verbose = true, std::vector<std::string> param_names = std::vector<std::string>(),
    bool reeval_likelihood = false, bool in_parallel = false, int n_threads = -1, 
    bool classic_gamma = false, 
    std::vector<std::vector<double>> init = std::vector<std::vector<double>>(), int init_iter = 0);

void Optimize_Priors(Randomizer& R, Likelihood& likelihood, MCMCReporter& report, std::vector<Distribution>& priors,
    unsigned int maxeval, double ftol_abs, bool verbose, bool in_parallel, unsigned int n_threads);

void NelderMead_Priors(Randomizer& R, Likelihood& likelihood, MCMCReporter& report, std::vector<Distribution>& priors,
    vector<double> initial, bool chebyshev, bool adaptive, bool perturb);

#endif