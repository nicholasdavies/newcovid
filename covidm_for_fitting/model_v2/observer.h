// observer.h

#ifndef OBSERVER_H
#define OBSERVER_H

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppGSL)]]

#include <vector>
#include <iostream>
#include <fstream>
#include <stdexcept>
#include <algorithm>
#include <numeric>
#include <ctime>
#include <limits>
#include <omp.h>
#include <Rcpp.h>
using namespace std;

struct Parameters;
struct PopulationParameters;
class Reporter;

// Observer function
struct Observer
{
    Observer() : func("rnorm"), null(true) { } // TODO CHANGE THIS TO STH ELSE
    Observer& operator=(Rcpp::Function f) { func = f; null = false; return *this; }
    bool operator()(Parameters* parent, PopulationParameters& params, double t, Reporter& rep);
    Rcpp::Function func;
    bool null;
};

#endif