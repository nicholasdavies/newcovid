#ifndef HELPER_H
#define HELPER_H

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppGSL)]]

#include <vector>
#include <stdexcept>
#include <numeric>

using namespace std;

struct Parameters;
class Randomizer;

//
// TIMING
//

// Clock functions
extern vector<double> ClockTimes;
extern double C0;

double Clock();
void StartClocking();
void ClockCheckpoint(unsigned int cp);
void ShowClockInfo();

//
// DISCRETE DISTRIBUTION
//

class MNApprox
{
public:
    static const unsigned int NVariants = 32;
    static const unsigned int VariantMask = 31;

    void Set(Parameters& P, Randomizer& Rand, vector<double>& p);


    void operator()(unsigned int N, vector<unsigned int>& out);

private:
    vector<vector<vector<unsigned int>>> x;
    vector<unsigned int> cycle;
};


// An arbitrary discrete distribution, used for maturation times in compartments
class Discrete
{
public:
    // Create distribution from set of unnormalised weights
    void operator=(std::vector<double> uw)
    {
        weights = uw;

        double total = accumulate(weights.begin(), weights.end(), 0.0);
        for (auto& w : weights)
            w /= total;

        storage.assign(weights.size(), 0);
    }

    // Calculate mean
    double Mean() const
    {
        double m = 0;
        for (unsigned int i = 0; i < weights.size(); ++i)
            m += i * weights[i];
        return m;
    }

    // Normalised weights
    std::vector<double> weights;

    // Storage for draws of integers
    std::vector<unsigned int> storage;

    // Fast approximate multinomial draws
    MNApprox mn_approx;
};


//
// MATRIX
//

struct Matrix
{
    Matrix()
     : nc(0)
    { }

    Matrix(double X, unsigned int nrow, unsigned int ncol)
     : x(nrow * ncol, X), nc(ncol)
    { }

    Matrix(vector<double> X, unsigned int nrow, unsigned int ncol)
     : x(X), nc(ncol)
    {
        if (x.size() != nrow * ncol)
            throw runtime_error("Improperly sized matrix.");
    }

    double& operator()(unsigned int i, unsigned int j)
    {
        return x[i * nc + j];
    }

    unsigned int NCol() { return nc; }
    unsigned int NRow() { return x.size() / nc; }

    vector<double> x;
    unsigned int nc;
};

#endif