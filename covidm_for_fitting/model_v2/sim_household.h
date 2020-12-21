// sim_household.h

#ifndef SIM_HOUSEHOLD_H
#define SIM_HOUSEHOLD_H

#include <vector>
#include <cstdint>
using namespace std;

//
// MODEL DYNAMICS
//

struct Parameters;
class Randomizer;
class Reporter;

// State for one individual
struct Individual
{
    // Individual properties
    uint32_t age, house, region, state, tminus;
    uint32_t group0, group1, group2, group3;
};

// Individual-based simulation of COVID-19 dynamics
class Households
{
public:
    Households(Parameters& P);

    // Run the model
    void Run(Parameters& P, Randomizer& Rand, Reporter& rep);

private:
    vector<Individual> pop;
    uint32_t n_houses;
    uint32_t n_regions;
    uint32_t n_groups;
    vector<double> x;
};

#endif