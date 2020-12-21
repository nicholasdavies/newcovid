// observer.cpp

#include "observer.h"
#include "reporter.h"
#include "parameters.h"

bool Observer::operator()(Parameters* parent, PopulationParameters& pp, double t, Reporter& rep)
{
    if (null)
        return true;

    (void) parent; (void) pp; (void) t; (void) rep;
    return true;
    
    // using namespace Rcpp;
    // RObject observer_return = func(t, rep.dynamics_df);

    // if (is<List>(observer_return))
    // {
    //     List ret = as<List>(observer_return);
    //     if (ret.containsElementNamed("changes"))
    //     {
    //         if (is<List>(ret["changes"]))
    //         {
    //             List changes = as<List>(ret["changes"]);
    //             for (unsigned int i = 0; i < changes.size(); ++i)
    //             {
    //                 string name = as<string>(as<CharacterVector>(changes.names())[i]);
    //                 RObject value = as<RObject>(changes[i]);
    //                 pp.Set(parent, name, value);
    //             }
    //             pp.Recalculate();
    //         }
    //         else
    //         {
    //             throw logic_error("changes must be a list with named entries.");
    //         }
    //     }
    //     if (ret.containsElementNamed("print"))
    //     {
    //         string message = as<string>(ret["print"]);
    //         cout << message << '\n';
    //     }
    //     if (ret.containsElementNamed("csv"))
    //     {
    //         string csv = as<string>(ret["csv"]);
    //         rep.csv += csv;
    //         if (csv.back() != '\n')
    //             rep.csv += '\n';
    //     }
    //     if (ret.containsElementNamed("halt"))
    //     {
    //         bool halt = as<bool>(ret["halt"]);
    //         return !halt;
    //     }
    // }
    // else if (!observer_return.isNULL())
    // {
    //     throw runtime_error("Observer function must return either a list or NULL.");
    // }
    // return true;
}