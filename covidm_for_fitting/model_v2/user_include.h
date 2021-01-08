// user_include.h

#ifndef USER_INCLUDE_H
#define USER_INCLUDE_H

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppGSL)]]

#include <vector>
#include <Rcpp.h>

struct Parameters;
class Reporter;
class Randomizer;

typedef void   (*cm_changes_func)       (const std::vector<double>&, Parameters&);
typedef double (*cm_loglikelihood_func) (const std::vector<double>&, Reporter&);
typedef bool   (*cm_observer_func)      (Parameters&, Randomizer&, Reporter&, double, std::vector<double>&);

#include "randomizer.h"
#include "convenience.h"
#include "parameters.h"
#include "reporter.h"

#endif
