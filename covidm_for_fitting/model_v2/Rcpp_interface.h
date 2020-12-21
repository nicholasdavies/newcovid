// Rcpp_interface.h

#ifndef RCPP_INTERFACE_H
#define RCPP_INTERFACE_H

#include "parameters.h"
#include "randomizer.h"
#include "reporter.h"

Reporter RunSimulation(Parameters& P, Randomizer& Rand, vector<double> x = vector<double>());

#endif
