# # v2/cm2_funcs.R
# # for building user defined functions as part of parameters
# 
# cm_compile_changes_func = function(code)
# {
#     include_covidm_types = paste0('#include "', normalizePath("./covidm_for_fitting/model_v2/user_include.h"), '"\n');
#     src = glue::glue(include_covidm_types, '
#     void user_changes_func(const std::vector<double>& x, Parameters& P)
#     {
#         (void) x; (void) P;
#         ${paste(code, collapse = "\n")}
#     }
#     
#     // [[Rcpp::export]]
#     Rcpp::XPtr<cm_changes_func> get_func()
#     {
#         return Rcpp::XPtr<cm_changes_func>(new cm_changes_func(&user_changes_func));
#     }
#     ', .sep = "\n", .open = "${", .close = "}");
#     
#     e = new.env();
#     Rcpp::sourceCpp(code = src, env = e);
#     fp = get("get_func", e)();
#     
#     return (list(src = src, fp = fp))
# }
# 
# cm_compile_loglikelihood_func = function(code)
# {
#     include_covidm_types = paste0('#include "', normalizePath("./covidm_for_fitting/model_v2/user_include.h"), '"\n');
#     src = glue::glue(include_covidm_types, '
#     double user_loglikelihood_func(const std::vector<double>& x, Reporter& dyn)
#     {
#         (void) x; (void) dyn;
#         double ll = 0.0;
#         ${paste(code, collapse = "\n")}
#         return ll;
#     }
#     
#     // [[Rcpp::export]]
#     Rcpp::XPtr<cm_loglikelihood_func> get_func()
#     {
#         return Rcpp::XPtr<cm_loglikelihood_func>(new cm_loglikelihood_func(&user_loglikelihood_func));
#     }
#     ', .sep = "\n", .open = "${", .close = "}");
#     
#     e = new.env();
#     Rcpp::sourceCpp(code = src, env = e);
#     fp = get("get_func", e)();
#     
#     return (list(src = src, fp = fp))
# }
# 
# cm_compile_observer_func = function(code)
# {
#     include_covidm_types = paste0('#include "', normalizePath("./covidm_for_fitting/model_v2/user_include.h"), '"\n');
#     src = glue::glue(include_covidm_types, '
#     bool user_observer_func(Parameters& P, Randomizer& R, Reporter& dyn, double t, std::vector<double>& x)
#     {
#         (void) P; (void) R; (void) dyn; (void) t; (void) x;
#         ${paste(code, collapse = "\n")}
#         return true;
#     }
#     
#     // [[Rcpp::export]]
#     Rcpp::XPtr<cm_observer_func> get_func()
#     {
#         return Rcpp::XPtr<cm_observer_func>(new cm_observer_func(&user_observer_func));
#     }
#     ', .sep = "\n", .open = "${", .close = "}");
#     
#     e = new.env();
#     Rcpp::sourceCpp(code = src, env = e);
#     fp = get("get_func", e)();
#     
#     return (list(src = src, fp = fp))
# }
