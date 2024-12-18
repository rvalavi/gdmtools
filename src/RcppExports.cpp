// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// par_cpp
Rcpp::NumericVector par_cpp(const Rcpp::NumericMatrix& rast_vals, const Rcpp::NumericMatrix& ref_vals, const Rcpp::NumericMatrix& samp_vals, const double intercept, int nthreads);
RcppExport SEXP _gdmtools_par_cpp(SEXP rast_valsSEXP, SEXP ref_valsSEXP, SEXP samp_valsSEXP, SEXP interceptSEXP, SEXP nthreadsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::NumericMatrix& >::type rast_vals(rast_valsSEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericMatrix& >::type ref_vals(ref_valsSEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericMatrix& >::type samp_vals(samp_valsSEXP);
    Rcpp::traits::input_parameter< const double >::type intercept(interceptSEXP);
    Rcpp::traits::input_parameter< int >::type nthreads(nthreadsSEXP);
    rcpp_result_gen = Rcpp::wrap(par_cpp(rast_vals, ref_vals, samp_vals, intercept, nthreads));
    return rcpp_result_gen;
END_RCPP
}
// presist_cpp
Rcpp::NumericVector presist_cpp(const Rcpp::NumericMatrix& rast_vals, const Rcpp::NumericMatrix& ref_vals, const Rcpp::NumericMatrix& cond_vals, const double intercept, const double power, int nthreads);
RcppExport SEXP _gdmtools_presist_cpp(SEXP rast_valsSEXP, SEXP ref_valsSEXP, SEXP cond_valsSEXP, SEXP interceptSEXP, SEXP powerSEXP, SEXP nthreadsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::NumericMatrix& >::type rast_vals(rast_valsSEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericMatrix& >::type ref_vals(ref_valsSEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericMatrix& >::type cond_vals(cond_valsSEXP);
    Rcpp::traits::input_parameter< const double >::type intercept(interceptSEXP);
    Rcpp::traits::input_parameter< const double >::type power(powerSEXP);
    Rcpp::traits::input_parameter< int >::type nthreads(nthreadsSEXP);
    rcpp_result_gen = Rcpp::wrap(presist_cpp(rast_vals, ref_vals, cond_vals, intercept, power, nthreads));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_gdmtools_par_cpp", (DL_FUNC) &_gdmtools_par_cpp, 5},
    {"_gdmtools_presist_cpp", (DL_FUNC) &_gdmtools_presist_cpp, 6},
    {NULL, NULL, 0}
};

RcppExport void R_init_gdmtools(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
