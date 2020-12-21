// sim_compartment.cpp

#include "sim_compartment.h"
#include "parameters.h"
#include "reporter.h"
#include "randomizer.h"
#include "user_defined.h"

//
// MODEL DYNAMICS
//

Population::Population(Parameters& P, unsigned int pindex)
 : seed_row(0), seed_row2(0), p(pindex)
{
    // Set up built-in compartments
    N  = P.pop[p].size;
    S  = N;
    R  = vector<double>(S.size(), 0.);
    R2 = vector<double>(S.size(), 0.);
    V  = vector<double>(S.size(), 0.);
    
    E  = vector<Compartment>(S.size());
    Ip = vector<Compartment>(S.size());
    Ia = vector<Compartment>(S.size());
    Is = vector<Compartment>(S.size());

    E2  = vector<Compartment>(S.size());
    Ip2 = vector<Compartment>(S.size());
    Ia2 = vector<Compartment>(S.size());
    Is2 = vector<Compartment>(S.size());

    // Initial immunity
    for (unsigned int a = 0; a < S.size(); ++a) {
        double imm = 0;
        if (P.deterministic) imm = S[a] * P.pop[p].imm0[a];
        else                 imm = (unsigned int)(S[a] * P.pop[p].imm0[a] + 0.5);
        S[a] -= imm;
        R[a] += imm;
    }

    // Set up user-specified processes
    unsigned int n_pc = 0;
    for (auto& p : P.processes)
        n_pc += p.ids.size();
    pc = vector<vector<Compartment>>(n_pc, vector<Compartment>(S.size()));
    pci = vector<double>(pc.size(), 0.);
    pco = vector<double>(pc.size(), 0.);
}

// Do seeding and calculate contagiousness
void Population::Contagiousness(Parameters& P, Randomizer& Rand, double t, vector<double>& contag, vector<double>& contag2)
{
    auto add = [&](unsigned int age, double n)
    {
        n = min(n, S[age]);
        S[age] -= n;
        E[age].Add(P, Rand, n, P.pop[p].dE);
    };

    auto add2 = [&](unsigned int age, double n)
    {
        n = min(n, S[age]);
        S[age] -= n;
        E2[age].Add(P, Rand, n, P.pop[p].dE);
    };

    // Do seeding
    while (seed_row < P.pop[p].seed_times.size() && t >= P.pop[p].seed_times[seed_row])
    {
        if (P.deterministic)
        {
            for (unsigned int a = 0; a < S.size(); ++a)
                add(a, P.pop[p].dist_seed_ages.weights[a]);
        }
        else
        {
            Rand.Multinomial(1, P.pop[p].dist_seed_ages.weights, P.pop[p].dist_seed_ages.storage);
            for (unsigned int a = 0; a < S.size(); ++a)
            {
                if (P.pop[p].dist_seed_ages.storage[a] == 1)
                {
                    add(a, 1);
                    break;
                }
            }
        }
        ++seed_row;
    }

    // Do seeding (2)
    while (seed_row2 < P.pop[p].seed_times2.size() && t >= P.pop[p].seed_times2[seed_row2])
    {
        if (P.deterministic)
        {
            for (unsigned int a = 0; a < S.size(); ++a)
                add2(a, P.pop[p].dist_seed_ages.weights[a]);
        }
        else
        {
            Rand.Multinomial(1, P.pop[p].dist_seed_ages.weights, P.pop[p].dist_seed_ages.storage);
            for (unsigned int a = 0; a < S.size(); ++a)
            {
                if (P.pop[p].dist_seed_ages.storage[a] == 1)
                {
                    add2(a, 1);
                    break;
                }
            }
        }
        ++seed_row2;
    }

    // Calculate contagiousness from this population
    for (unsigned int a = 0; a < contag.size(); ++a)
        contag[a]  = (N[a] == 0) ? 0 : (P.pop[p].fIp[a] * Ip[a].Size()  + P.pop[p].fIa[a] * Ia[a].Size()  + P.pop[p].fIs[a] * Is[a].Size() ) / N[a];
    for (unsigned int a = 0; a < contag2.size(); ++a)
        contag2[a] = (N[a] == 0) ? 0 : (P.pop[p].fIp[a] * Ip2[a].Size() + P.pop[p].fIa[a] * Ia2[a].Size() + P.pop[p].fIs[a] * Is2[a].Size()) / N[a];
}

// Execute one time step's events
void Population::Tick(Parameters& P, Randomizer& Rand, double t, vector<double>& infec, vector<double>& infec2, Reporter& rep)
{
    // Calculate force of infection in this compartment
    lambda.assign(infec.size(), 0.0);
    lambda2.assign(infec.size(), 0.0);
    for (unsigned int a = 0; a < lambda.size(); ++a)
    {
        for (unsigned int b = 0; b < lambda.size(); ++b)
        {
            lambda[a] += P.pop[p].u[a] * P.pop[p].cm(a,b) * infec[b];
            lambda2[a] += P.pop[p].u2[a] * P.pop[p].cm(a,b) * infec2[b];
        }
    }

    // Account for seasonality
    if (P.pop[p].season_A[0] != 0)
    {
        double f = 1.0 + P.pop[p].season_A[0] * cos(2. * M_PI * (t - P.pop[p].season_phi[0]) / P.pop[p].season_T[0]);
        for (unsigned int a = 0; a < lambda.size(); ++a)
        {
            lambda[a]  *= f;
            lambda2[a] *= f;
        }
    }

    // Account for importation
    // NOTE - assuming no importation of second strain.
    for (unsigned int a = 0; a < lambda.size(); ++a)
        lambda[a] += P.pop[p].omega[a];

    // Helpers
    auto multinomial = [&](double n, vector<double>& p, vector<double>& nd_out, vector<unsigned int>& ni_out) {
        nd_out.resize(p.size(), 0.);
        if (P.deterministic)
        {
            for (unsigned int i = 0; i < p.size(); ++i)
                nd_out[i] = n * p[i];
        }
        else
        {
            ni_out.resize(p.size(), 0);
            Rand.Multinomial(n, p, ni_out);
            for (unsigned int i = 0; i < p.size(); ++i)
                nd_out[i] = ni_out[i];
        }
    };

    auto poisson = [&](double l) {
        if (P.deterministic)
            return l;
        else
            return (double)Rand.Poisson(l);
    };

    auto binomial = [&](double n, double p) {
        if (P.deterministic)
            return n * p;
        else
            return (double)Rand.Binomial(n, p);
    };

    auto num = [&](double n) {
        if (P.deterministic)
            return n;
        else
            return round(n);
    };

    // Do state transitions and reporting for each age group
    for (unsigned int a = 0; a < lambda.size(); ++a)
    {
        // 0. Report prevalences
        if (t == (int)t)
        {
            // Built-in states
            rep(t, p, a, 0) = S[a];
            rep(t, p, a, 1) = E[a].Size();
            rep(t, p, a, 2) = Ip[a].Size();
            rep(t, p, a, 3) = Is[a].Size();
            rep(t, p, a, 4) = Ia[a].Size();
            rep(t, p, a, 5) = R[a];
            rep(t, p, a, 6) = E2[a].Size();
            rep(t, p, a, 7) = Ip2[a].Size();
            rep(t, p, a, 8) = Is2[a].Size();
            rep(t, p, a, 9) = Ia2[a].Size();
            rep(t, p, a, 10) = R2[a];
            rep(t, p, a, 11) = V[a];
            rep(t, p, a, 12) = lambda[a];
            rep(t, p, a, 13) = lambda2[a];

            // User-specified processes
            for (auto& process : P.processes)
            {
                for (unsigned int i = 0; i < process.p_cols.size(); ++i)
                    rep(t, p, a, process.p_cols[i]) = pc[process.p_ids[i]][a].Size();
            }

        }

        // 1. Built-in states

        // Vaccination and waning of natural immunity and vaccine protection
        // S -> V, V -> S, R -> S, R2 -> S
        double nS_V  = min(S[a],       num(P.pop[p].v[a] * S[a] / N[a] * P.time_step));
        double nV_S  = binomial(V[a],  1.0 - exp(-P.pop[p].wv[a]  * P.time_step));
        double nR_S  = binomial(R[a],  1.0 - exp(-P.pop[p].wn[a]  * P.time_step));
        double nR2_S = binomial(R2[a], 1.0 - exp(-P.pop[p].wn2[a] * P.time_step));

        S[a]  -= nS_V;
        V[a]  += nS_V;
        V[a]  -= nV_S;
        S[a]  += nV_S;
        R[a]  -= nR_S;
        S[a]  += nR_S;
        R2[a] -= nR2_S;
        S[a]  += nR2_S;

        // Infection rates
        // S -> E, S -> E2
        double lambda_12 = lambda[a] + lambda2[a];
        double nS_E12 = binomial(S[a], 1.0 - exp(-lambda_12 * P.time_step));
        double nS_E = binomial(nS_E12, lambda_12 > 0 ? lambda[a] / lambda_12 : 0);
        double nS_E2 = nS_E12 - nS_E;

        // Reinfection rates
        // R -> E, R -> E2, R2 -> E, R2 -> E2, V -> E, V -> E2
        double foi_r   = lambda[a]  * (1 - P.pop[p].pi_r[a]);
        double foi2_r  = lambda2[a] * (1 - P.pop[p].pi2_r[a]);
        double foi_r2  = lambda[a]  * (1 - P.pop[p].pi_r2[a]);
        double foi2_r2 = lambda2[a] * (1 - P.pop[p].pi2_r2[a]);
        double foi_v   = lambda[a]  * (1 - P.pop[p].ei_v[a]);
        double foi2_v  = lambda2[a] * (1 - P.pop[p].ei2_v[a]);
        
        double nR_E12 = binomial(R[a], 1.0 - exp(-(foi_r + foi2_r) * P.time_step));
        double nR_E = binomial(nR_E12, foi_r + foi2_r > 0 ? foi_r / (foi_r + foi2_r) : 0);
        double nR_E2 = nR_E12 - nR_E;
        double nR2_E12 = binomial(R2[a], 1.0 - exp(-(foi_r2 + foi2_r2) * P.time_step));
        double nR2_E = binomial(nR2_E12, foi_r2 + foi2_r2 > 0 ? foi_r2 / (foi_r2 + foi2_r2) : 0);
        double nR2_E2 = nR2_E12 - nR2_E;
        double nV_E12 = binomial(V[a], 1.0 - exp(-(foi_v + foi2_v) * P.time_step));
        double nV_E = binomial(nV_E12, foi_v + foi2_v > 0 ? foi_v / (foi_v + foi2_v) : 0);
        double nV_E2 = nV_E12 - nV_E;
        
        // Infection / reinfection changes
        S[a]  -= nS_E  + nS_E2;
        R[a]  -= nR_E  + nR_E2;
        R2[a] -= nR2_E + nR2_E2;
        V[a]  -= nV_E  + nV_E2;

        // Add new exposed. If dE2 has been supplied, use that as latent period for strain 2, otherwise use dE.
        E[a] .Add(P, Rand, nS_E  + nR_E  + nR2_E  + nV_E,  P.pop[p].dE);
        if (P.pop[p].dE2.weights.size() <= 1)
            E2[a].Add(P, Rand, nS_E2 + nR_E2 + nR2_E2 + nV_E2, P.pop[p].dE);
        else
            E2[a].Add(P, Rand, nS_E2 + nR_E2 + nR2_E2 + nV_E2, P.pop[p].dE2);

        // Strain 1
        // E -> Ip/Ia
        double nE_Ipa = E[a].Mature();
        double nE_Ip = binomial(nE_Ipa, P.pop[p].y[a]);
        double nE_Ia = nE_Ipa - nE_Ip;
        Ip[a].Add(P, Rand, nE_Ip, P.pop[p].dIp);
        Ia[a].Add(P, Rand, nE_Ia, P.pop[p].dIa);

        // Ip -> Is
        double nIp_Is = Ip[a].Mature();
        Is[a].Add(P, Rand, nIp_Is, P.pop[p].dIs);

        // Is -> R
        double nIs_R = Is[a].Mature();
        R[a] += nIs_R;

        // Ia -> R
        double nIa_R = Ia[a].Mature();
        R[a] += nIa_R;

        // Strain 2
        // E2 -> Ip2/Ia2
        double nE2_Ipa2 = E2[a].Mature();
        double nE2_Ip2 = binomial(nE2_Ipa2, P.pop[p].y2[a]);
        double nE2_Ia2 = nE2_Ipa2 - nE2_Ip2;
        Ip2[a].Add(P, Rand, nE2_Ip2, P.pop[p].dIp);
        Ia2[a].Add(P, Rand, nE2_Ia2, P.pop[p].dIa);

        // Ip2 -> Is2
        double nIp2_Is2 = Ip2[a].Mature();
        Is2[a].Add(P, Rand, nIp2_Is2, P.pop[p].dIs);

        // Is2 -> R2
        double nIs2_R2 = Is2[a].Mature();
        R2[a] += nIs2_R2;

        // Ia2 -> R2
        double nIa2_R2 = Ia2[a].Mature();
        R2[a] += nIa2_R2;

        // 2. User-specified processes
        fill(pco.begin(), pco.end(), -1.);

        for (auto& process : P.processes)
        {
            // Determine number of individuals entering the process
            double n_entering = 0.;
            switch (process.source_id)
            {
                case src_newE:
                    n_entering = nS_E + nR_E + nR2_E + nV_E; break;
                case src_newI:
                    n_entering = nE_Ipa; break;
                case src_newIp:
                    n_entering = nE_Ip; break;
                case src_newIs:
                    n_entering = nIp_Is; break;
                case src_newIa:
                    n_entering = nE_Ia; break;
                case src_newE2:
                    n_entering = nS_E2 + nR_E2 + nR2_E2 + nV_E2; break;
                case src_newI2:
                    n_entering = nE2_Ipa2; break;
                case src_newIp2:
                    n_entering = nE2_Ip2; break;
                case src_newIs2:
                    n_entering = nIp2_Is2; break;
                case src_newIa2:
                    n_entering = nE2_Ia2; break;
                case src_newEE2:
                    n_entering = nS_E + nR_E + nR2_E + nV_E + nS_E2 + nR_E2 + nR2_E2 + nV_E2; break;
                case src_newII2:
                    n_entering = nE_Ipa + nE2_Ipa2; break;
                case src_newIpIp2:
                    n_entering = nE_Ip + nE2_Ip2; break;
                case src_newIsIs2:
                    n_entering = nIp_Is + nIp2_Is2; break;
                case src_newIaIa2:
                    n_entering = nE_Ia + nE2_Ia2; break;
                default:
                    n_entering = pco[process.source_id];
                    if (n_entering < 0)
                        throw logic_error("Process sourced from unset user process. Have user processes been specified in the right order?");
                    break;
            }

            multinomial(n_entering, process.prob[a], nd_out, ni_out);

            // Seed and mature this process's compartments
            unsigned int c = 0;
            for (unsigned int compartment_id : process.ids)
            {
                if (compartment_id != Null)
                {
                    pc[compartment_id][a].Add(P, Rand, nd_out[c], process.delays[c]);
                    pci[compartment_id] = nd_out[c];
                    pco[compartment_id] = pc[compartment_id][a].Mature();
                }
                ++c;
            }
        }

        // 3. Report incidence / outcidence

        // User-specified processes
        for (auto& process : P.processes)
        {
            for (unsigned int i = 0; i < process.i_cols.size(); ++i)
                rep(t, p, a, process.i_cols[i]) += pci[process.i_ids[i]];
            for (unsigned int i = 0; i < process.o_cols.size(); ++i)
                rep(t, p, a, process.o_cols[i]) += pco[process.o_ids[i]];
        }
    }

    // Births, deaths, aging
    double Ntot = accumulate(N.begin(), N.end(), 0.0);
    for (unsigned int a = N.size() - 1; ; --a)
    {
        // Births
        double B = poisson(Ntot * (exp(P.pop[p].B[a] * P.time_step) - 1.));

        // Deaths (outer compartments)
        double death_prob = 1.0 - exp(-P.pop[p].D[a] * P.time_step);
        double DS   = binomial(S[a],  death_prob);
        double DR   = binomial(R[a],  death_prob);
        double DR2  = binomial(R2[a], death_prob);
        double DV   = binomial(V[a],  death_prob);

        // Changes
        N[a] += B;
        S[a] += B;
        S[a]  -= DS;
        R[a]  -= DR;
        R2[a] -= DR2;
        V[a]  -= DV;
        
        // Deaths (inner compartments)
        double DE   = E[a]  .RemoveProb(P, Rand, death_prob);
        double DIp  = Ip[a] .RemoveProb(P, Rand, death_prob);
        double DIa  = Ia[a] .RemoveProb(P, Rand, death_prob);
        double DIs  = Is[a] .RemoveProb(P, Rand, death_prob);
        double DE2  = E2[a] .RemoveProb(P, Rand, death_prob);
        double DIp2 = Ip2[a].RemoveProb(P, Rand, death_prob);
        double DIa2 = Ia2[a].RemoveProb(P, Rand, death_prob);
        double DIs2 = Is2[a].RemoveProb(P, Rand, death_prob);

        N[a]  -= DS + DR + DR2 + DV + DE + DIp + DIa + DIs + DE2 + DIp2 + DIa2 + DIs2;

        // Agings
        if (a != lambda.size() - 1)
        {
            double age_prob = 1.0 - exp(-P.pop[p].A[a] * P.time_step);
            double AS   = binomial(S[a],  age_prob);
            double AR   = binomial(R[a],  age_prob);
            double AR2  = binomial(R2[a], age_prob);
            double AV   = binomial(V[a],  age_prob);

            S[a]      -= AS;
            S[a + 1]  += AS;
            R[a]      -= AR;
            R[a + 1]  += AR;
            R2[a]     -= AR2;
            R2[a + 1] += AR2;
            V[a]      -= AV;
            V[a + 1]  += AV;
            double AE   = E[a]  .MoveProb(E   [a + 1], P, Rand, age_prob);
            double AIp  = Ip[a] .MoveProb(Ip  [a + 1], P, Rand, age_prob);
            double AIa  = Ia[a] .MoveProb(Ia  [a + 1], P, Rand, age_prob);
            double AIs  = Is[a] .MoveProb(Is  [a + 1], P, Rand, age_prob);
            double AE2  = E2[a] .MoveProb(E2  [a + 1], P, Rand, age_prob);
            double AIp2 = Ip2[a].MoveProb(Ip2 [a + 1], P, Rand, age_prob);
            double AIa2 = Ia2[a].MoveProb(Ia2 [a + 1], P, Rand, age_prob);
            double AIs2 = Is2[a].MoveProb(Is2 [a + 1], P, Rand, age_prob);

            N[a]      -= AS + AR + AR2 + AV + AE + AIp + AIa + AIs + AE2 + AIp2 + AIa2 + AIs2;
            N[a + 1]  += AS + AR + AR2 + AV + AE + AIp + AIa + AIs + AE2 + AIp2 + AIa2 + AIs2;
        }

        if (a == 0)
            break;
    }
}

// Print full population details
void Population::DebugPrint() const
{
    auto vecprint = [&](const vector<double>& vec, string name) {
        cout << name;
        for (auto& v : vec)
            cout << " " << v;
        cout << "\n";
    };

    auto comprint = [&](const vector<Compartment>& comp, string name) {
        cout << name;
        for (unsigned int c = 0; c < comp.size(); ++c) {
            cout << "element " << c << "\n";
            comp[c].DebugPrint();
        }
    };

    vecprint(lambda, "lambda");
    vecprint(lambda2, "lambda2");
    vecprint(N, "N");
    vecprint(S, "S");
    vecprint(R, "R");
    vecprint(V, "V");
    comprint(E, "E");
    comprint(Ip, "Ip");
    comprint(Ia, "Ia");
    comprint(Is, "Is");
    vecprint(R2, "R2");
    comprint(E2, "E2");
    comprint(Ip2, "Ip2");
    comprint(Ia2, "Ia2");
    comprint(Is2, "Is2");
    cout << "seed_row " << seed_row << " p " << p << "\n";
    cout << "seed_row2 " << seed_row2 << " p " << p << "\n";
    for (auto& c : pc)
        comprint(c, "User");

    cout << "\n\n";
}


Metapopulation::Metapopulation(Parameters& P)
{
    P.changes.Capture(P);

    for (unsigned int i = 0; i < P.pop.size(); ++i)
        pops.push_back(Population(P, i));
}

// Execute one time step's events
bool Metapopulation::Tick(Parameters& P, Randomizer& Rand, double t, unsigned int ts, Reporter& rep)
{
    // Apply any changes to parameters
    P.changes.Apply(P, t);

    unsigned int n_ages = P.pop[0].size.size();

    // Calculate contagiousness from each population
    // NOTE -- 'contag' subscripted first by j, then by a.
    // It's the effective number of infectious individuals FROM subpop j of age a.
    contag.assign(pops.size(), vector<double>(n_ages, 0.0));
    contag2.assign(pops.size(), vector<double>(n_ages, 0.0));
    for (unsigned int j = 0; j < pops.size(); ++j)
        pops[j].Contagiousness(P, Rand, t, contag[j], contag2[j]);

    // note -- 'infec' subscripted first by i, then by a
    // It's the effective number of infectious individuals who are CURRENTLY IN subpop i of age a.
    infec.assign(pops.size(), vector<double>(n_ages, 0.0));
    infec2.assign(pops.size(), vector<double>(n_ages, 0.0));
    for (unsigned int i = 0; i < pops.size(); ++i)
        for (unsigned int j = 0; j < pops.size(); ++j)
            for (unsigned int a = 0; a < n_ages; ++a)
            {
                infec[i][a]  += P.travel(j, i) * contag[j][a]  * (j != i ? P.pop[j].tau[a] : 1.0);
                infec2[i][a] += P.travel(j, i) * contag2[j][a] * (j != i ? P.pop[j].tau[a] : 1.0);
            }

    // Update populations
    //#pragma omp parallel for schedule(dynamic) reduction(&&:keep_going)
    for (unsigned int i = 0; i < pops.size(); ++i)
        pops[i].Tick(P, Rand, t, infec[i], infec2[i], rep);

    // Run observer at the last time step of each day.
    if (t + P.time_step == int(t + P.time_step))
        return CppObserver(P, Rand, rep, (int)t, x);

    return true;
}

void Metapopulation::Run(Parameters& P, Randomizer& Rand, Reporter& rep, vector<double> x_fit)
{
    x = x_fit;

    // Run simulation
    unsigned int time_steps = (1 + P.time1 - P.time0) / P.time_step;
    for (unsigned int ts = 0; ts < time_steps; ++ts)
    {
        if (!Tick(P, Rand, P.time0 + ts * P.time_step, ts, rep))
            break;
    }
}
