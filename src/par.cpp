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
Rcpp::NumericVector par_cpp(
        const Rcpp::NumericMatrix &rast_vals,
        const Rcpp::NumericMatrix &ref_vals,
        const Rcpp::NumericMatrix &samp_vals,
        const double intercept,
        int nthreads = -1)
{
    Lightweight_matrix<double> rast(rast_vals);
    Lightweight_matrix<double> refs(ref_vals);
    Lightweight_matrix<double> samples(samp_vals);

    int nr = rast.nrow();
    int nref = refs.nrow();
    int nsam = samples.nrow();

    // output vector
    std::vector<double> output(nr);

    // set the number of threads
    #ifdef _OPENMP
        if (nthreads < 1) nthreads = omp_get_max_threads();
        omp_set_num_threads(nthreads);
    #endif

    #pragma omp parallel for
    for (int i = 0; i < nr; i++)
    {
        // loop through the reference samples
        double pa_dist = 0.0;
        for (int j = 0; j < nref; j++)
        {
            auto dist_1 = distance(rast, refs, i, j);
            pa_dist += similarity(intercept, dist_1);
        }

        // loop through the reference cell for whole region
        double reg_dist = 0.0;
        for (int j = 0; j < nsam; j++)
        {
            auto dist_2 = distance(rast, samples, i, j);
            reg_dist += similarity(intercept, dist_2);
        }

        #pragma omp critical
        output[i] = (pa_dist / static_cast<double>(nref)) / (reg_dist / static_cast<double>(nsam));
    }

    return Rcpp::wrap(output);
}

