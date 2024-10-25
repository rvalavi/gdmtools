// [[Rcpp::plugins(openmp)]]
// [[Rcpp::plugins("cpp11")]]
#include <Rcpp.h>
#include <cmath>
#include <vector>
#include "Lightweight_matrix.h"

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace Rcpp;

inline auto dist(
        const Lightweight_matrix<double> &x,
        const Lightweight_matrix<double> &y,
        int i,
        int j)
{
    auto out = 0.0;
    for (int k; k < x.ncol(); k ++)
    {
        out += std::abs(x(i, k) - y(j, k));
    }

    return out;
}

// [[Rcpp::export]]
Rcpp::NumericVector PAR(
        const Rcpp::NumericMatrix &rast_vals,
        const Rcpp::NumericMatrix &ref_vals,
        const std::vector<int> &samples,
        const double intercept,
        int num_threads = -1)
{
    Lightweight_matrix<double> rast(rast_vals);
    Lightweight_matrix<double> refs(ref_vals);

    int nr = rast.nrow();
    int nref = refs.nrow();
    int nsam = samples.size();

    // output vector
    std::vector<double> output(nr);

    // set the number of threads
    #ifdef _OPENMP
        if (num_threads < 1) num_threads = omp_get_max_threads();
        omp_set_num_threads(num_threads);
    #endif

    #pragma omp parallel for
    for (int i = 0; i < nr; i++)
    {
        // loop through the reference samples
        double pa_dist = 0.0;
        for (int j; j < nref; j++)
        {
            auto dist_1 = dist(rast, refs, i, j);
            pa_dist += std::exp(-1 * (intercept + dist_1));
        }

        // loop through the reference cell for whole region
        double reg_dist = 0.0;
        for (int j; j < nsam; j++)
        {
            auto dist_2 = dist(rast, refs, i, samples[j]);
            reg_dist += std::exp(-1 * (intercept + dist_2));
        }

        output[i] = (pa_dist / static_cast<double>(nref)) / (reg_dist / static_cast<double>(nsam));
    }

    return Rcpp::wrap(output);
}

