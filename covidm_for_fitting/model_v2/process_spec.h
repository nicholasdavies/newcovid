//
// PROCESS
//

#ifndef PROCESS_SPEC_H
#define PROCESS_SPEC_H

#include "helper.h"

enum SourceID
{
    src_newE = 1000000,
    src_newI,
    src_newIp,
    src_newIs,
    src_newIa,
    
    src_newE2,
    src_newI2,
    src_newIp2,
    src_newIs2,
    src_newIa2,
    
    src_newEE2,
    src_newII2,
    src_newIpIp2,
    src_newIsIs2,
    src_newIaIa2
};

const unsigned int Null = 999999;
class Discrete;

struct ProcessSpec
{
    unsigned int source_id; // identifier for source compartment
    string source_name;     // name of the above
    string type;            // ignored for now - multinomial or dirichlet multinomial

    vector<string> names;       // names of sub-processes
    vector<unsigned int> ids;   // identifiers of sub-process compartments
    vector<string> report;      // reporting mode of sub-processes: empty, "i", "o", "p", or a combination of these
    vector<vector<double>> prob;// probability by group of entering each sub-process from the source above; indexed by group then by subprocess
    vector<Discrete> delays;    // delays for each sub-process

    vector<unsigned int> p_cols;
    vector<unsigned int> p_ids;
    vector<unsigned int> i_cols;
    vector<unsigned int> i_ids;
    vector<unsigned int> o_cols;
    vector<unsigned int> o_ids;     // prevalence, incidence, and outcidence data columns and sub-process identifiers for this process
};

#endif
