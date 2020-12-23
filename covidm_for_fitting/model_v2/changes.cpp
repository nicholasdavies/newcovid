// changes.cpp

#include "changes.h"
#include "parameters.h"
#include "helper.h"
#include <stdexcept>

#define PVector(x)   x
#define PDiscrete(x) x.weights

#define ParamCapture(param, get_vec) \
    else if (param_name == #param) { \
        if (pops.empty()) { \
            for (auto& pp : P.pop) { \
                param_ptr.push_back(&get_vec(pp.param)); \
                param_orig.push_back(get_vec(pp.param)); \
            } \
        } else { \
            for (auto& pi : pops) { \
                if (pi >= P.pop.size()) \
                    throw logic_error("Referring to nonexistent population in Change::Capture."); \
                param_ptr.push_back(&get_vec(P.pop[pi].param)); \
                param_orig.push_back(get_vec(P.pop[pi].param)); \
            } \
        } \
    }

// Construct a change impacting parameter pname in populations po of parameters P;
// apply value v with mode m at times t
Change::Change(Parameters& P, vector<unsigned int>& po, vector<unsigned int>& ru, string pname,
    Mode m, vector<double>& t, vector<vector<double>>& v)
 : mode(m), times(t), values(v), param_name(pname), current(-1), pops(po), runs(ru)
{
    if (times.size() != values.size())
        throw logic_error("Change: times and values must be same length.");
}

// Capture parameters to change
void Change::Capture(Parameters& P)
{
    current = -1;
    if (false) {}
    ParamCapture(dE,      PDiscrete)
    ParamCapture(dIp,     PDiscrete)
    ParamCapture(dIa,     PDiscrete)
    ParamCapture(dIs,     PDiscrete)
    ParamCapture(dE2,     PDiscrete)
    ParamCapture(contact, PVector)
    ParamCapture(contact_mult,    PVector)
    ParamCapture(contact_lowerto, PVector)
    ParamCapture(u,       PVector)
    ParamCapture(u2,   PVector)
    ParamCapture(fIp,     PVector)
    ParamCapture(fIa,     PVector)
    ParamCapture(fIs,     PVector)
    ParamCapture(y,       PVector)
    ParamCapture(y2,      PVector)
    ParamCapture(omega,   PVector)
    ParamCapture(tau,     PVector)
    ParamCapture(pi_r,    PVector)
    ParamCapture(pi_r2,   PVector)
    ParamCapture(pi2_r,   PVector)
    ParamCapture(pi2_r2,  PVector)
    ParamCapture(wn,      PVector)
    ParamCapture(wn2,     PVector)
    ParamCapture(v,       PVector)
    ParamCapture(wv,      PVector)
    ParamCapture(ei_v,    PVector)
    ParamCapture(ei2_v,   PVector)
    ParamCapture(ed_vi,   PVector)
    ParamCapture(ed_vi2,  PVector)
    ParamCapture(A,       PVector)
    ParamCapture(B,       PVector)
    ParamCapture(D,       PVector)
    ParamCapture(season_A,      PVector)
    ParamCapture(season_T,      PVector)
    ParamCapture(season_phi,    PVector)

    else if (param_name == "travel") {
        param_ptr.push_back(&P.travel.x);
        param_orig.push_back(P.travel.x);
    }
    else
        throw logic_error("Unsupported parameter change for " + param_name);
}

// Return true if parameters will change at time t
bool Change::Update(double t)
{
    bool refresh = false;
    while (current + 1 < (int)times.size() && t >= times[current + 1])
    {
        ++current;
        refresh = true;
    }
    return refresh;
}

// Apply parameter changes
void Change::Apply(double t)
{
    auto op = [&](vector<double>* ptr) {
        if (ptr->size() == values[current].size() || values[current].empty())
        {
            if (!values[current].empty())
            {
                for (unsigned int j = 0; j < ptr->size(); ++j)
                {
                    if (values[current][j] == -999)
                        continue;

                    switch (mode)
                    {
                        case Assign:    (*ptr)[j] = values[current][j]; break;
                        case Add:       (*ptr)[j] += values[current][j]; break;
                        case Multiply:  (*ptr)[j] *= values[current][j]; break;
                        case LowerTo:   (*ptr)[j] = min((*ptr)[j], values[current][j]); break;
                        case RaiseTo:   (*ptr)[j] = max((*ptr)[j], values[current][j]); break;
                        case Bypass:    break;
                    }
                }
            }
        }
        else
        {
            cout << current << "\n";
            cout << ptr->size() << "\n";
            cout << values[current].size() << "\n";
            throw logic_error("Change::Apply: change value has different size from parameter.");
        }
    };

    if (current >= 0)
        for (auto ptr : param_ptr)
            op(ptr);
}

// Reset all linked parameters to their original values
void Change::Reset()
{
    for (unsigned int i = 0; i < param_ptr.size(); ++i)
        *param_ptr[i] = param_orig[i];
}

// Capture parameters to change
void ChangeSet::Capture(Parameters& P)
{
    for (auto& c : ch)
        c.Capture(P);
}


// Apply any needed changes
void ChangeSet::Apply(Parameters& P, double t)
{
    bool refresh = false;

    for (auto& c : ch)
    {
        if (c.Update(t))
        {
            c.Reset();
            refresh = true;
        }
    }

    if (refresh)
    {
        for (auto& c : ch)
            c.Apply(t);
        for (auto& pp : P.pop) { // TODO -- slightly wasteful; Change could keep track if linked parameter sets need refreshing; then again in most cases probably needed
            pp.needs_recalc = true;
            pp.Recalculate();
        }
    }
}
