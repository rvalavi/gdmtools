#pragma once

#ifndef SIMILARITY_HELPER
#define SIMILARITY_HELPER

inline double distance(
        const Lightweight_matrix<double> &x,
        const Lightweight_matrix<double> &y,
        int i,
        int j)
{
    double dist = 0.0;
    for (int k = 0; k < x.ncol(); k ++)
    {
        dist += std::abs(x(i, k) - y(j, k));
    }

    return dist;
}

// calculate the similarity
inline double similarity(double x, double y)
{
    double sim = static_cast<double>(std::exp(-1.0 * (x + y)));

    return sim;
}

#endif /* SIMILARITY_HELPER */
