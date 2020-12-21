// convenience.h

#ifndef CONVENIENCE_H
#define CONVENIENCE_H

#include <vector>

struct Parameters;
class Reporter;

// TODO perform rigorous error checking on inputs

// Helper functions for user code
std::vector<double> seq(double x0, double x1, double by = 1);

// binomial log density
double binom(double x, double size, double prob);

// negative binomial log density
double nbinom(unsigned int x, double mean, double size);

// beta binomial log density
double bbinom(double k, double n, double p, double a_plus_b);

// negative binomial log density with retrospective confirmation
double nbinom_gammaconf(unsigned int x, double mean, double size, double days_ago, double conf_delay_mean, double conf_delay_shape);

// normal log density
double norm(double x, double mean, double sd);

// skew normal log density
double skewnorm(double x, double xi, double omega, double alpha);

// beta density
double dbeta(double x, double alpha, double beta);

// construct a delay distribution following a gamma distribution with mean mu and shape parameter shape.
std::vector<double> delay_gamma(double mu, double shape, double t_max, double t_step, double mult = 1.);

// construct a delay distribution following a lognormal distribution with true mean mu and coefficient of variation cv.
std::vector<double> delay_lnorm(double mu, double cv, double t_max, double t_step);

// estimate the basic reproduction number
double estimate_R0(Parameters& P, Reporter& rep, double t, unsigned int p, unsigned int iter);

// estimate the effective reproduction number
double estimate_Rt(Parameters& P, Reporter& rep, double t, unsigned int p, unsigned int iter);

// clamp a number between two limits
double clamp(double x, double x0 = 0.0, double x1 = 1.0);

// smootherstep function
double smootherstep(double x0, double x1, double y0, double y1, double x);

#endif