// sim_household.cpp

#include "sim_household.h"
#include "parameters.h"
#include "reporter.h"
#include "randomizer.h"
#include "user_defined.h"
#include <vector>
#include <iostream>
#include <fstream>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
using namespace std;

// Infection state
namespace State
{
    enum : uint32_t { S, E, Ip, Is, Ia, R };
}

// Cache of random discrete distribution
struct CacheDiscrete
{
    CacheDiscrete(gsl_rng* r, std::vector<double> prob, uint32_t n) 
     : i(0), size(n)
    {
        // Choose number of each item
        std::vector<uint32_t> each(prob.size(), 0);
        gsl_ran_multinomial(r, prob.size(), size, &prob[0], &each[0]);

        // Shuffle items
        draws.resize(size);
        for (unsigned int i = 0, j = 0; j < each.size(); i += each[j++])
            fill(draws.begin() + i, draws.begin() + i + each[j], j);
        gsl_ran_shuffle(r, &draws[0], draws.size(), sizeof(decltype(draws[0])));
    }

    uint32_t operator()()
    {
        if (++i >= size) 
            i = 0;
        return draws[i];
    }

    std::vector<uint32_t> draws;
    uint32_t i;
    uint32_t size;
};

// Cache of random continuous distribution
struct CacheContinuous
{
    template <typename Generator>
    CacheContinuous(Generator gen, uint32_t n) 
     : draws(n, 0.0), i(0), size(n)
    {
        for (unsigned int k = 0; k < n; ++k)
            draws[k] = gen();
    }

    uint32_t operator()()
    {
        if (++i >= size) 
            i = 0;
        return draws[i];
    }

    std::vector<double> draws;
    uint32_t i;
    uint32_t size;
};

// Household model implementation
Households::Households(Parameters& P)
 : n_houses(0), n_regions(0), n_groups(0)
{
    P.changes.Capture(P);
    // TODO specify in parameters
    ifstream fin("/Users/nick/Documents/household/ukhouse.txt");
    fin.ignore(numeric_limits<streamsize>::max(), '\n');

    Individual i;
    i.state = State::S;
    i.tminus = 0;
    i.group0 = i.group1 = i.group2 = i.group3 = 0;
    n_groups = 1;

    while (true)
    {
        fin >> i.age >> i.house >> i.region;
        fin >> i.group0 >> i.group1 >> i.group2 >> i.group3;

        if (!fin.good()) break;
        pop.push_back(i);

        if (n_houses < i.house + 1)
            n_houses = i.house + 1;

        if (n_regions < i.region + 1)
            n_regions = i.region + 1;
    }
    cout << "Read " << pop.size() << " individuals.\n";
}

void Households::Run(Parameters& P, Randomizer& Rand, Reporter& rep)
{
    gsl_rng* r = Rand.GSL_RNG();

    const uint32_t n_age_groups = P.pop[0].size.size();
    const uint32_t p = 0; // TEMP - population indicator

    // Simulation parameters
    vector<double> f(6 * n_age_groups, 0.0);
    copy(P.pop[0].fIp.begin(), P.pop[0].fIp.end(), f.begin() + 2 * n_age_groups);
    copy(P.pop[0].fIs.begin(), P.pop[0].fIs.end(), f.begin() + 3 * n_age_groups);
    copy(P.pop[0].fIa.begin(), P.pop[0].fIa.end(), f.begin() + 4 * n_age_groups);
    auto beta_h = [](double n) { return 0; }; /// TODO - no household transmission at the moment.

    // Construct travel distribution
    vector<CacheDiscrete> travel;
    for (uint32_t i = 0; i < n_regions; ++i)
        travel.push_back(CacheDiscrete(r, vector<double>(&P.travel(i, 0), &P.travel(i, n_regions)), 5131));

    CacheDiscrete dur_E (r, P.pop[p].dE.weights, 1111);
    CacheDiscrete dur_Ip(r, P.pop[p].dIp.weights, 1111);
    CacheDiscrete dur_Is(r, P.pop[p].dIs.weights, 1111);
    CacheDiscrete dur_Ia(r, P.pop[p].dIa.weights, 1111);

    vector<CacheDiscrete> y_trial;
    for (uint32_t ag = 0; ag < n_age_groups; ++ag)
        y_trial.push_back(CacheDiscrete(r, { 1. - P.pop[p].y[ag], P.pop[p].y[ag] }, 223));

    // Force of infection storage: effective proportion of individuals by household, by group, and by region & age, who are infectious
    vector<double> infn_h(n_houses, 0);
    vector<double> infn_g(n_groups, 0);
    vector<double> infp_c(n_regions * n_age_groups, 0);
    vector<double> infn_h_next(infn_h);
    vector<double> infn_g_next(infn_g);
    vector<double> infp_c_next(infp_c);

    // Incidence / prevalence storage by region
    vector<vector<uint32_t>>  inc(n_regions, vector<uint32_t>(6 * n_age_groups, 0));
    vector<vector<uint32_t>> prev(n_regions, vector<uint32_t>(6 * n_age_groups, 0));

    // Seed infection
    for (uint32_t seed = 0; seed < pop.size(); seed += 100000)
    {
        pop[seed].state = State::E;
        pop[seed].tminus = dur_E();
    }

    // Cache of exponential variates
    CacheContinuous exponential([&]() { return gsl_ran_exponential(r, 1.0 / P.time_step); }, 111);

    // Loop through each time step
    for (unsigned int ts = 0; ts < (1 + P.time1 - P.time0) / P.time_step; ++ts)
    {
        double t = P.time0 + ts * P.time_step;

        // Apply any changes to parameters for this time step
        P.changes.Apply(P, t);

        // Zero out prevalence and incidence storage at beginning of each day
        if (t == (int)t)
        {
            prev.assign(n_regions, vector<uint32_t>(6 * n_age_groups, 0));
            inc.assign (n_regions, vector<uint32_t>(6 * n_age_groups, 0));
        }

        // d is the distance, measured in units of force of infection per person per day, to the next individual to infect.
        double d = exponential();
        for (uint32_t i = 0; i < pop.size(); ++i)
        {
            // Infection
            double foi = P.pop[p].u[pop[i].age] * (
                beta_h(infn_h[pop[i].house]) + 
                infn_g[pop[i].group0] + infn_g[pop[i].group1] + infn_g[pop[i].group2] + infn_g[pop[i].group3] + 
                infp_c[pop[i].region * n_age_groups + pop[i].age]
            );
            d -= foi;

            if (d < 0)
            {
                d = exponential();
                if (pop[i].state == State::S)
                {
                    pop[i].state = State::E;
                    pop[i].tminus = dur_E();
                    ++inc[pop[i].region][pop[i].age * 6 + pop[i].state];
                }
            }

            // Movement between compartments (non-infection)
            if (pop[i].tminus == 0)
            {
                switch (pop[i].state)
                {
                    case State::E: // Either go to Ip or Ia
                        if (y_trial[pop[i].age]())
                        {
                            pop[i].state = State::Ip;
                            pop[i].tminus = dur_Ip();
                        }
                        else
                        {
                            pop[i].state = State::Ia;
                            pop[i].tminus = dur_Ia();
                        }
                        break;

                    case State::Ip:
                        pop[i].state = State::Is;
                        pop[i].tminus = dur_Is();
                        break;

                    case State::Is:
                        pop[i].state = State::R;
                        break;

                    case State::Ia:
                        pop[i].state = State::R;
                        break;

                    default:
                        break;
                }
                ++inc[pop[i].region][pop[i].age * 6 + pop[i].state];
            }
            --pop[i].tminus;

            // Contribution to next time step's FOI
            infn_h_next[pop[i].house] += f[pop[i].state * n_age_groups + pop[i].age];
            infn_g_next[pop[i].group0] += f[pop[i].state * n_age_groups + pop[i].age];
            infn_g_next[pop[i].group1] += f[pop[i].state * n_age_groups + pop[i].age];
            infn_g_next[pop[i].group2] += f[pop[i].state * n_age_groups + pop[i].age];
            infn_g_next[pop[i].group3] += f[pop[i].state * n_age_groups + pop[i].age];
            infp_c_next[travel[pop[i].region]() * n_age_groups + pop[i].age] += f[pop[i].state * n_age_groups + pop[i].age];

            // Track prevalence
            if (t == (int)t)
                ++prev[pop[i].region][pop[i].age * 6 + pop[i].state];
        }

        // Put household, group, and regional FOI in place for next time step
        swap(infn_h_next, infn_h);
        swap(infn_g_next, infn_g);
        for (uint32_t r = 0; r < n_regions; ++r)
            for (uint32_t a = 0; a < n_age_groups; ++a)
                infp_c_next[r * n_age_groups + a] /= P.pop[r].size[a];

        for (uint32_t r = 0; r < n_regions; ++r)
        {
            for (uint32_t a = 0; a < n_age_groups; ++a)
            {
                infp_c[r * n_age_groups + a] = 0;
                for (uint32_t b = 0; b < n_age_groups; ++b)
                    infp_c[r * n_age_groups + a] += P.pop[p].cm(a, b) * infp_c_next[r * n_age_groups + b];
            }
        }

        infn_h_next.assign(infn_h_next.size(), 0.);
        infn_g_next.assign(infn_g_next.size(), 0.);
        infp_c_next.assign(infp_c_next.size(), 0.);

        // Record prevalences
        if (t == (int)t)
        {
            for (uint32_t region = 0; region < n_regions; ++region)
            {
                for (uint32_t a = 0, i = 0; a < n_age_groups; ++a)
                {
                    rep(t, region, a, 0) = prev[region][i++];
                    rep(t, region, a, 1) = prev[region][i++];
                    rep(t, region, a, 2) = prev[region][i++];
                    rep(t, region, a, 3) = prev[region][i++];
                    rep(t, region, a, 4) = prev[region][i++];
                    rep(t, region, a, 5) = prev[region][i++];
                }
            }
        }
        // Record incidences
        if (t + P.time_step == int(t + P.time_step))
        {
            for (uint32_t region = 0; region < n_regions; ++region)
            {
                for (uint32_t a = 0; a < n_age_groups; ++a)
                {
                    uint32_t i = a * 6;
                    rep(t, region, a, 6) = inc[region][i + 3]; // cases
                    rep(t, region, a, 8) = inc[region][i + 4]; // subclinical
                }
            }
        }

        // Run observer at the last time step of each day.
        if (t + P.time_step == int(t + P.time_step))
            if (!CppObserver(P, Rand, rep, (int)t, x))
                return;
    }
}

// to do
// X multiple regions
// X age specific clinical proportion and susceptibility: age lookup table
// X age specific mixing within a region
// X travel to other regions
// age specific travel to other regions
