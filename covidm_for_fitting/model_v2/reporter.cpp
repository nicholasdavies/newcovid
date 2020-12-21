// reporter.cpp

#include "reporter.h"
#include "parameters.h"
#include <numeric>

Reporter::Reporter(Parameters& P)
 : t0(P.time0),
   n_times((unsigned int)(P.time1 - P.time0) + 1),
   n_populations(P.pop.size()),
   n_age_groups(P.pop[0].size.size()),
   col_names({ "S", "E", "Ip", "Is", "Ia", "R", "E2", "Ip2", "Is2", "Ia2", "R2", "V", "foi", "foi2" })
{
    if (P.time_step != 1. / P.report_every)
        throw logic_error("Reporter requires P.time_step = 1 / P.report_every.");

    // Create space for built-in compartments
    for (unsigned int c = 0; c < col_names.size(); ++c)
        data.push_back(vector<double>(n_times * n_populations * n_age_groups, 0.));

    // User-defined compartments
    for (auto& pr : P.processes)
    {
        // Reset existing references within parameters to columns
        // TODO this code should be integrated partly into Parameters, so there are no problems with multiple
        // instances of Reporter using the same set of Parameters and running into crosstalk issues.
        pr.p_cols.clear();
        pr.p_ids.clear();
        pr.i_cols.clear();
        pr.i_ids.clear();
        pr.o_cols.clear();
        pr.o_ids.clear();

        for (unsigned int j = 0; j < pr.names.size(); ++j)
        {
            for (unsigned int k = 0; k < pr.report[j].size(); ++k)
            {
                col_names.push_back(pr.names[j] + "_" + pr.report[j][k]);
                data.push_back(vector<double>(n_times * n_populations * n_age_groups, 0.));

                if (pr.report[j][k] == 'p') {
                    pr.p_cols.push_back(data.size() - 1);
                    pr.p_ids.push_back(pr.ids[j]);
                } else if (pr.report[j][k] == 'i') {
                    pr.i_cols.push_back(data.size() - 1);
                    pr.i_ids.push_back(pr.ids[j]);
                } else if (pr.report[j][k] == 'o') {
                    pr.o_cols.push_back(data.size() - 1);
                    pr.o_ids.push_back(pr.ids[j]);
                } else {
                    throw runtime_error("Unrecognized process report type " + string(1, pr.report[j][k]) + ".");
                }
            }
        }
    }
}


// Access data, summed over populations and groups
double Reporter::operator()(string compartment, double t, initializer_list<unsigned int> pops, initializer_list<unsigned int> groups)
{
    // Locate compartment
    auto cc = find(col_names.begin(), col_names.end(), compartment);
    if (cc == col_names.end()) {
        cout << "No such compartment " << compartment << "\n" << flush;
        throw logic_error("No such compartment: " + compartment);
    }
    unsigned int c = cc - col_names.begin();

    // Define populations and groups to loop over
    // TODO optimize...
    vector<unsigned int> pp = pops;
    vector<unsigned int> gg = groups;
    if (pp.empty())
    {
        pp.assign(n_populations, 0);
        iota(pp.begin(), pp.end(), 0);
    }
    if (gg.empty())
    {
        gg.assign(n_age_groups, 0);
        iota(gg.begin(), gg.end(), 0);
    }

    // Sum over populations and groups
    double x = 0.0;
    for (auto p : pp)
        for (auto g : gg)
            x += (*this)(t, p, g, c);

    return x;
}

// Save data to file
void Reporter::Save(string basename, unsigned long int seed)
{
    ofstream fout(basename + ".txt");

    // Output header
    fout << "# t0 " << t0 << " n_times " << n_times << " n_populations " << n_populations << " n_groups " << n_age_groups << " seed " << seed << "\n";
    for (unsigned int cd = 0; cd < col_names.size(); ++cd)
        fout << col_names[cd] << (cd == col_names.size() - 1 ? "" : "\t");

    if (obs.size() > 0)
        for (unsigned int co = 0; co < obs.size(); ++co)
            fout << "\tobs" << co;

    fout << "\n";

    // Output data
    for (unsigned int r = 0; r < data[0].size(); ++r)
    {
        for (unsigned int cd = 0; cd < data.size(); ++cd)
            fout << data[cd][r] << (cd == data.size() - 1 ? "" : "\t");

        if (obs.size() > 0)
            for (unsigned int co = 0; co < obs.size(); ++co)
                fout << "\t" << obs[co][r];

        fout << "\n";
    }
    fout.close();

    // Output csv
    if (!csv.empty())
    {
        ofstream fcsv(basename + ".csv");
        fcsv << csv;
    }
}
