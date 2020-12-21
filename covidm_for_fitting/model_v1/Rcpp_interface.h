// Rcpp_interface.h

// [[Rcpp::export]]
Rcpp::List cm_backend_simulate(Rcpp::List parameters, unsigned int n_run = 1, unsigned long int seed = 0)
{
    // Initialise parameters for this simulation
    Randomizer Rand(seed);
    Parameters covidm_parameters;
    SetParameters(covidm_parameters, parameters, Rand);

    Rcpp::List dynamics;
    Rcpp::List csvs;

    for (unsigned int r = 0; r < n_run; ++r)
    {
        // Run the simulation
        Parameters P = covidm_parameters;
        Reporter rep = RunSimulation(P, Rand);

        dynamics.push_back(rep.dynamics_df);
        csvs.push_back(rep.csv);
    }

    return Rcpp::List::create(
        Rcpp::Named("dynamics") = dynamics,
        Rcpp::Named("csv") = csvs
    );
}

// [[Rcpp::export]]
Rcpp::DataFrame cm_evaluate_distribution(string dist_code, unsigned int steps = 101, double xmin = 0, double xmax = -1)
{
    Distribution dist(dist_code);
    if (xmax < xmin)
    {
        xmin = dist.LowerBound();
        xmax = dist.UpperBound();
    }

    vector<double> x(steps, 0.);
    vector<double> p(steps, 0.);

    for (unsigned int s = 0; s < steps; ++s)
    {
        x[s] = xmin + ((xmax - xmin) / (steps - 1.)) * s;
        p[s] = exp(dist.LogProbability(x[s]));
    }

    Rcpp::DataFrame results = Rcpp::DataFrame::create(
        Rcpp::Named("x") = x,
        Rcpp::Named("p") = p
    );

    return results;
}


// Fitting interface

#include "mcmc.h"

class Caller
{
public:
    Caller(Rcpp::Function ff, Rcpp::List ep)
     : f(ff), p(ep), i(0) {}

    double operator()(const vector<double>& theta, int& obs)
    {
        (void)obs;
        return Rcpp::as<double>(f(theta, p, ++i));
    }

private:
    Rcpp::Function f;
    Rcpp::List p;
    unsigned int i;
};

class MCMCReporter
{
public:
    MCMCReporter(unsigned int iterations, unsigned int chains, vector<string>& param_names)
     : n_chains(chains), n_samp(iterations * n_chains), 
       trial(n_samp, 0), chain(n_samp, 0), lp(n_samp, 0.0), ll(n_samp, 0.0), 
       theta(n_samp, param_names.size()), pnames(param_names)
    {
    }

    void operator()(int TR, double LP, double CH, double LL, vector<double> TH, int& OBS)
    {
        (void)OBS;

        unsigned int i = TR * n_chains + CH;

        if (TR < 0)
            i = 0;

        trial[i] = TR;
        chain[i] = CH;
        lp[i] = LP;
        ll[i] = LL;
        for (unsigned int j = 0; j < TH.size(); ++j)
            theta(i, j) = TH[j];
    }

    Rcpp::DataFrame Contents()
    {
        Rcpp::DataFrame df = Rcpp::DataFrame::create(
            Rcpp::Named("trial") = trial,
            Rcpp::Named("lp") = lp,
            Rcpp::Named("chain") = chain,
            Rcpp::Named("ll") = ll
        );

        for (unsigned int j = 0; j < pnames.size(); ++j)
            df.push_back(theta.column(j), pnames[j]);

        return Rcpp::DataFrame(df);
    }

    Rcpp::List ListContents()
    {
        Rcpp::List li = Rcpp::List::create(
            Rcpp::Named("trial") = trial[0],
            Rcpp::Named("lp") = lp[0],
            Rcpp::Named("chain") = chain[0],
            Rcpp::Named("ll") = ll[0]
        );

        for (unsigned int j = 0; j < pnames.size(); ++j)
            li.push_back(theta(j, 0), pnames[j]);

        return Rcpp::List(li);
    }

private:
    unsigned int n_chains, n_samp;
    Rcpp::IntegerVector trial, chain;
    Rcpp::NumericVector lp, ll;
    Rcpp::NumericMatrix theta;
    vector<string> pnames;
};

// [[Rcpp::export]]
Rcpp::DataFrame cm_backend_mcmc(Rcpp::Function likelihood, Rcpp::List extra_params, Rcpp::List params_priors,
    int seed, unsigned int burn_in, unsigned int n_chains, unsigned int iterations, bool verbose, 
    bool reeval_likelihood, bool in_parallel, int n_threads)
{
    /// TODO implement seed
    (void) seed;

    vector<string> params = Rcpp::as<vector<string>>(params_priors.names());
    vector<Distribution> priors;

    for (unsigned int i = 0; i < params_priors.size(); ++i)
        priors.push_back(Distribution(Rcpp::as<string>(params_priors[i])));

    Randomizer rand;
    Caller lik(likelihood, extra_params);
    MCMCReporter rep(iterations, n_chains, params);

    DEMCMC_Priors<int>(rand, lik, rep, burn_in, iterations, n_chains, priors, verbose, params, 
        reeval_likelihood, in_parallel, n_threads, false);

    return rep.Contents();
}

// [[Rcpp::export]]
Rcpp::DataFrame cm_backend_mcmc_init(Rcpp::Function likelihood, Rcpp::List extra_params,
    Rcpp::List params_priors, Rcpp::NumericMatrix initial,
    int seed, unsigned int burn_in, unsigned int n_chains, unsigned int iterations, bool verbose, 
    bool reeval_likelihood, bool in_parallel, int n_threads)
{
    /// TODO implement seed
    (void) seed;

    vector<string> params = Rcpp::as<vector<string>>(params_priors.names());
    vector<Distribution> priors;

    for (unsigned int i = 0; i < params_priors.size(); ++i)
        priors.push_back(Distribution(Rcpp::as<string>(params_priors[i])));

    vector<vector<double>> init;

    if (initial.nrow() > 0)
    {
        for (unsigned int c = 0; c < n_chains; ++c)
        {
            if ((unsigned int)initial.nrow() <= c)
                throw logic_error("initial list must be large enough to supply all chains");
            Rcpp::NumericVector nv = initial.row(c); 
            vector<double> v = Rcpp::as<vector<double>>(nv);
            if (v.size() != (unsigned int)params_priors.size())
                throw logic_error("elements of initial list must be same size as number of parameters");
            init.push_back(v);
        }
    }

    Randomizer rand;
    Caller lik(likelihood, extra_params);
    MCMCReporter rep(iterations, n_chains, params);

    DEMCMC_Priors<int>(rand, lik, rep, burn_in, iterations, n_chains, priors, verbose, params, 
        reeval_likelihood, in_parallel, n_threads, true, -1, init);

    return rep.Contents();
}

// [[Rcpp::export]]
Rcpp::List cm_backend_optimize(Rcpp::Function likelihood, Rcpp::List extra_params, Rcpp::List params_priors,
    int seed, unsigned int maxeval, double ftol_abs, bool verbose)
{
    (void) seed;

    vector<string> params = Rcpp::as<vector<string>>(params_priors.names());
    vector<Distribution> priors;

    for (unsigned int i = 0; i < params_priors.size(); ++i)
        priors.push_back(Distribution(Rcpp::as<string>(params_priors[i])));

    Randomizer rand;
    Caller lik(likelihood, extra_params);
    MCMCReporter rep(1, 1, params);

    Optimize_Priors<int>(rand, lik, rep, priors, maxeval, ftol_abs, verbose);

    return rep.ListContents();
}

// [[Rcpp::export]]
Rcpp::NumericVector cm_backend_prior_sample(Rcpp::List params_priors)
{
    Randomizer rand;
    vector<string> params = Rcpp::as<vector<string>>(params_priors.names());
    Rcpp::NumericVector prior_sample;

    for (unsigned int i = 0; i < params_priors.size(); ++i)
    {
        Distribution dist(Rcpp::as<string>(params_priors[i]));
        prior_sample.push_back(dist.Random(rand), params[i]);
    }

    return prior_sample;
}
