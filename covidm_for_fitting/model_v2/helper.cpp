// helper.cpp

#include "helper.h"

#include <gsl/gsl_sf_gamma.h>
#include <iostream>
#include "randomizer.h"
#include "parameters.h"

//
// TIMING
//

// Clock functions
vector<double> ClockTimes;
double C0;

double Clock()
{ 
    return double(clock()) / CLOCKS_PER_SEC;
}

void StartClocking()
{
    ClockTimes.clear();
    C0 = Clock();
}

void ClockCheckpoint(unsigned int cp)
{
    double C1 = Clock();
    if (cp >= ClockTimes.size()) ClockTimes.resize(cp + 1, 0.0);
    ClockTimes[cp] += C1 - C0;
    C0 = C1;
}

void ShowClockInfo()
{
    double total = 0;
    for (unsigned int i = 0; i < ClockTimes.size(); ++i)
    {
        cout << "Checkpoint " << i << ": " << ClockTimes[i] << "\n";
        total += ClockTimes[i];
    }
    cout << "Total time: " << total << " seconds.\n";
}




void MNApprox::Set(Parameters& P, Randomizer& Rand, vector<double>& p)
{
    x.assign(32, vector<vector<unsigned int>>(NVariants, vector<unsigned int>(p.size(), 0)));
    cycle.assign(NVariants, 0);
    for (unsigned int shift = 0; shift < 32; ++shift)
        for (unsigned int variant = 0; variant < NVariants; ++variant)
            gsl_ran_multinomial(Rand.GSL_RNG(), p.size(), 1 << shift, &p[0], &x[shift][variant][0]);
}


void MNApprox::operator()(unsigned int N, vector<unsigned int>& out)
{
    out.assign(out.size(), 0);
    for (unsigned int shift = 0; shift < 32; ++shift)
    {
        if ((N >> shift) & 1)
        {
            for (unsigned int i = 0; i < out.size(); ++i)
            {
                out[i] += x[shift][cycle[shift]][i];
            }
            cycle[shift] = (cycle[shift] + 1) & VariantMask;
        }
    }
}
