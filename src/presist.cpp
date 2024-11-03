// [[Rcpp::plugins(openmp)]]
// [[Rcpp::plugins("cpp11")]]
#include <Rcpp.h>
#include <cmath>
#include <vector>
#include "Lightweight_matrix.h"
#include "helper.h"
#ifdef _OPENMP
#include <omp.h>
#endif

using namespace Rcpp;


// [[Rcpp::export]]
Rcpp::NumericVector presist(
        const Rcpp::NumericMatrix &rast_vals,
        const Rcpp::NumericMatrix &ref_vals,
        const Rcpp::NumericMatrix &cond_vals,
        const double intercept,
        const double power = 0.25,
        int nthreads = -1)
{
    Lightweight_matrix<double> rast(rast_vals);
    Lightweight_matrix<double> refs(ref_vals);
    Lightweight_matrix<double> cond(cond_vals);

    int nr = rast.nrow();
    int nref = refs.nrow();
    // output vector
    std::vector<double> output(nr);

        // set the number of threads
    #ifdef _OPENMP
        if (nthreads < 1) nthreads = omp_get_max_threads();
        omp_set_num_threads(nthreads);
    #endif

    // define the variables for the total species persistence
    double numer = 0.0;
    double denom = 0.0;

    #pragma omp parallel for
    for (int i = 0; i < nr; i++)
    {
        double sim_cond = 0.0;
        double sim_pers = 0.0;

        // loop through the reference samples
        for (int j = 0; j < nref; j++)
        {
            auto eco_dist = dist(rast, refs, i, j);
            double sim_ij = similarity(intercept, eco_dist);
            sim_cond += (cond(j, 0) * sim_ij);
            sim_pers += sim_ij;
        }

        // dividing the sum similarity condition by the sum similarity pristine
        // and take to the power of z (here = 0.25) for species-area conversion
        double prop_pres = std::pow(sim_cond / sim_pers, power);
        double w = 1.0 / sim_pers;

        numer += (prop_pres * w);
        denom += w;

        #pragma omp critical
        output[i] = prop_pres;
    }

    Rcpp::Rcout << "The proportion of species persist: " << numer / denom << "\n";

    return Rcpp::wrap(output);
}

