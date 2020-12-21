// changes.h

#ifndef CHANGES_H
#define CHANGES_H

#include <vector>
#include <string>
using namespace std;


// Forward declarations
struct Parameters;
struct PopulationParameters;

struct Change
{
    // Operation to apply with change
    enum Mode { Assign, Add, Multiply, LowerTo, RaiseTo, Bypass };

    // Construct a change impacting parameter pname in populations po of parameters P;
    // apply value v with mode m at times t
    Change(Parameters& P, vector<unsigned int>& po, vector<unsigned int>& ru, string pname, 
        Mode m, vector<double>& t, vector<vector<double>>& v);

    // Capture parameters to change
    void Capture(Parameters& P);

    // Return true if parameters will change at time t
    bool Update(double t);

    // Apply parameter changes
    void Apply(double t);

    // Reset all linked parameters to their original values
    void Reset();

    Mode mode;
    vector<double> times;
    vector<vector<double>> values;
    
    string param_name;
    int current;
    vector<unsigned int> pops;
    vector<unsigned int> runs;
    vector<vector<double>*> param_ptr;
    vector<vector<double>> param_orig;
};

class ChangeSet
{
public:
    vector<Change> ch;

    // Capture parameters to change
    void Capture(Parameters& P);

    // Apply any needed changes
    void Apply(Parameters& P, double t);
};

#endif
