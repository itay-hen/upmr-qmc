// edd.h -- log|f[E_0,...,E_L]| for f(x) = exp(-beta x), the exponential divided difference
// used as the UPMR kernel. The sign of f[E_0..E_L] is (-1)^L.
//
// Method: with y_i = beta*(max_j E_j - E_i) >= 0,
//   f[E_0..E_L] = (-beta)^L exp(-beta max E) * g[y_0..y_L],  g(y) = exp(y),
// and g[y_0..y_L] = [exp(A)]_{0L} for the Opitz matrix A = diag(y) + N (N = superdiagonal of ones).
// We use the balanced matrix A_s = diag(y) + s N, for which [exp(A_s)]_{0L} = s^L [exp(A)]_{0L},
// and sum the Taylor series of exp(A_s) e_L. All matrix entries are nonnegative, so every term is
// nonnegative (no cancellation). Without balancing (s = 1) the components of the Taylor iterate
// span a range of order y_max^L, which exceeds double precision for L >~ 180 at beta >= 5; with
// s = max(1, y_max) the range is of order 2^k instead. The iterate is renormalized whenever its
// largest entry leaves [1e-150, 1e150], and the logarithm of the scale is accumulated separately.
#pragma once
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <vector>

static double log_edd(const std::vector<double>& E, double beta) {
    const int L = (int)E.size() - 1;
    if (L == 0) return -beta * E[0];
    double ymax = -1e300;
    for (double e : E) ymax = std::max(ymax, beta * e);
    std::vector<double> a(L + 1), v(L + 1, 0.0), t(L + 1);
    double amax = 0;
    for (int i = 0; i <= L; i++) { a[i] = ymax - beta * E[i]; amax = std::max(amax, a[i]); }
    const double s = std::max(1.0, amax);
    const double N = amax + s + 1.0;          // the series converges once k >> amax + s
    v[L] = 1.0;
    double logscale = 0.0, S = 0.0;
    for (int k = 1;; k++) {
        double mx = 0;
        for (int i = 0; i < L; i++) { t[i] = (a[i] * v[i] + s * v[i + 1]) / k; mx = std::max(mx, t[i]); }
        t[L] = a[L] * v[L] / k; mx = std::max(mx, t[L]);
        v.swap(t);
        if (mx > 1e150 || (mx > 0 && mx < 1e-150)) {
            for (double& x : v) x /= mx;
            logscale += std::log(mx); S /= mx; mx = 1.0;
        }
        if (k >= L) {
            S += v[0];
            if (k > 2 * N && mx < 1e-18 * S) break;
            if (k > 1000000) { fprintf(stderr, "log_edd: no convergence (q=%d); aborting\n", L); exit(2); }
        }
    }
    const double res = L * std::log(beta) - ymax + std::log(S) + logscale - L * std::log(s);
    if (!(S > 0) || !std::isfinite(S) || !std::isfinite(logscale) || !std::isfinite(res)) {
        fprintf(stderr, "log_edd: non-finite arithmetic (q=%d, S=%g, logscale=%g); aborting\n", L, S, logscale);
        exit(2);
    }
    return res;
}
