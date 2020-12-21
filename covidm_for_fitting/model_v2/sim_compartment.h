// sim_compartment.h

#ifndef SIM_COMPARTMENT_H
#define SIM_COMPARTMENT_H

#include <vector>
using namespace std;
class Compartment;

//
// MODEL DYNAMICS
//

struct Parameters;
class Randomizer;
class Reporter;

// A population of individuals, with SEI3HR dynamics.
class Population
{
public:
    // Construct a population with the specified size by age group; initially all uninfected
    Population(Parameters& P, unsigned int pindex);

    // Do seeding and calculate contagiousness
    void Contagiousness(Parameters& P, Randomizer& Rand, double t, vector<double>& contag, vector<double>& contag2);

    // Execute one time step's events
    void Tick(Parameters& P, Randomizer& Rand, double t, vector<double>& infec, vector<double>& infec2, Reporter& rep);

    // Print full population details
    void DebugPrint() const;

//private:
    vector<double> lambda;
    vector<double> lambda2;
    vector<double> N, S, R, R2, V;              // Total number, susceptible, recovered, recovered 2, vaccinated
    vector<Compartment> E, Ip, Ia, Is;          // Strain 1 exposed, presymptomatic, asymptomatic, symptomatic
    vector<Compartment> E2, Ip2, Ia2, Is2;      // Strain 2 exposed, presymptomatic, asymptomatic, symptomatic
    unsigned int seed_row, seed_row2;           // Which seed event is next
    unsigned int p;                             // Which population this is
    vector<vector<Compartment>> pc;             // User-specified process compartments, indexed by process id, then group
    vector<unsigned int> ni_out;                // Temporary storage
    vector<double> nd_out;                      // Temporary storage
    vector<double> pci;                         // Temporary storage
    vector<double> pco;                         // Temporary storage
};

// A metapopulation, containing multiple subpopulations.
class Metapopulation
{
public:
    Metapopulation(Parameters& P);

    // Execute one time step's events
    bool Tick(Parameters& P, Randomizer& Rand, double t, unsigned int ts, Reporter& rep);

    // Run the model
    void Run(Parameters& P, Randomizer& Rand, Reporter& rep, vector<double> x_fit = vector<double>());

//private:
    vector<vector<double>> contag, contag2;
    vector<vector<double>> infec, infec2;
    vector<Population> pops;
    vector<double> x;
};

#endif
