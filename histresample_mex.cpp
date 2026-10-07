// histresample_mex.cpp
// Multinomial resampling of a histogram by sequential binomial splitting.
// newhist = histresample_mex(inputhist, n, niter, seed)
//   inputhist : nbins vector of non-negative counts
//   n         : photons per resampled histogram
//   niter     : number of resampled histograms
//   seed      : RNG seed (pass randi(2^32-1) to tie it to MATLAB's rng)
// Returns nbins x niter double matrix. Cost is O(nbins * niter), independent of n.
//
// Compile: mex -O histresample_mex.cpp

#include "mex.h"
#include <random>
#include <vector>
#include <cstdint>

// Exact binomial sampler (Knuth TAOCP Vol 2, 3.4.1): split large n with a
// beta-distributed order statistic, then finish with Bernoulli trials.
// O(log n) per draw. (std::binomial_distribution in libstdc++ is biased
// for mid-range n, so it is not used.)
static long long binom_exact(long long n, double p, std::mt19937_64 &rng)
{
    std::uniform_real_distribution<double> unif(0.0, 1.0);
    long long k = 0;
    while (n > 32) {
        long long a = 1 + n / 2;
        long long b = n + 1 - a;
        std::gamma_distribution<double> ga((double) a, 1.0), gb((double) b, 1.0);
        double x1 = ga(rng), x2 = gb(rng);
        double x = x1 / (x1 + x2);   // a-th order statistic of n uniforms
        if (x >= p) {
            n = a - 1;
            p = p / x;
        } else {
            k += a;
            n = b - 1;
            p = (p - x) / (1 - x);
        }
    }
    for (long long i = 0; i < n; i++)
        if (unif(rng) < p) k++;
    return k;
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    if (nrhs != 4)
        mexErrMsgIdAndTxt("histresample_mex:nrhs", "Usage: histresample_mex(inputhist, n, niter, seed)");

    const double *hist = mxGetPr(prhs[0]);
    const mwSize nbins = mxGetNumberOfElements(prhs[0]);
    const long long n = (long long) mxGetScalar(prhs[1]);
    const mwSize niter = (mwSize) mxGetScalar(prhs[2]);
    const uint64_t seed = (uint64_t) mxGetScalar(prhs[3]);

    // Conditional probability of bin i given the photon is in bins i..end
    std::vector<double> pcond(nbins, 0.0);
    double tail = 0;
    for (mwSize i = nbins; i-- > 0;) {
        if (hist[i] < 0)
            mexErrMsgIdAndTxt("histresample_mex:neg", "inputhist must be non-negative");
        tail += hist[i];
        pcond[i] = tail > 0 ? hist[i] / tail : 0.0;
    }

    plhs[0] = mxCreateDoubleMatrix(nbins, niter, mxREAL);
    double *out = mxGetPr(plhs[0]);

    std::mt19937_64 rng(seed);
    for (mwSize k = 0; k < niter; k++) {
        double *col = out + k * nbins;
        long long remaining = n;
        for (mwSize i = 0; i < nbins && remaining > 0; i++) {
            if (pcond[i] <= 0) continue;
            long long draw;
            if (pcond[i] >= 1) {
                draw = remaining;
            } else {
                draw = binom_exact(remaining, pcond[i], rng);
            }
            col[i] = (double) draw;
            remaining -= draw;
        }
    }
}
