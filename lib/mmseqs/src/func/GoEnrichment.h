#ifndef MMSEQS_GOENRICHMENT_H
#define MMSEQS_GOENRICHMENT_H

#include <vector>
#include <algorithm>
#include <cmath>
#include <cstring>

namespace GoEnrichment {

// log of binomial coefficient C(n,k), via lgamma for numerical stability
static inline double logComb(int n, int k) {
    if (k < 0 || k > n) return -1e300;
    return lgamma(n + 1.0) - lgamma(k + 1.0) - lgamma(n - k + 1.0);
}

// One-sided hypergeometric p-value P(X >= k), X ~ Hyper(N, n, M)  [enrichment]
// N = background total, n = background with GO, M = study size, k = study with GO
static inline double hypergeomPval(int N, int n, int M, int k) {
    if (k <= 0 || n <= 0 || M <= 0 || N <= 0) return 1.0;
    if (k > n || k > M) return 0.0;
    double logDenom = logComb(N, M);
    double pval = 0.0;
    int xMax = std::min(n, M);
    for (int x = k; x <= xMax; ++x) {
        double logNum = logComb(n, x) + logComb(N - n, M - x);
        if (logNum > logDenom + 1e-9) { pval = 1.0; break; }
        pval += std::exp(logNum - logDenom);
        if (pval >= 1.0) { pval = 1.0; break; }
    }
    return pval;
}

// One-sided hypergeometric p-value P(X <= k), X ~ Hyper(N, n, M)  [depletion]
static inline double hypergeomPvalLower(int N, int n, int M, int k) {
    if (n <= 0 || M <= 0 || N <= 0) return 1.0;
    double logDenom = logComb(N, M);
    double pval = 0.0;
    int xMin = std::max(0, M - (N - n));
    for (int x = xMin; x <= k; ++x) {
        double logNum = logComb(n, x) + logComb(N - n, M - x);
        if (logNum > logDenom + 1e-9) { pval = 1.0; break; }
        pval += std::exp(logNum - logDenom);
        if (pval >= 1.0) { pval = 1.0; break; }
    }
    return pval;
}

// Benjamini-Hochberg FDR correction (returns adjusted p-values, same order as input)
static inline std::vector<double> bhCorrect(const std::vector<double>& pvals) {
    int m = (int)pvals.size();
    if (m == 0) return {};
    std::vector<std::pair<double, int>> sorted(m);
    for (int i = 0; i < m; ++i) sorted[i] = {pvals[i], i};
    std::sort(sorted.begin(), sorted.end());
    std::vector<double> adj(m, 1.0);
    double runMin = 1.0;
    for (int i = m - 1; i >= 0; --i) {
        double q = sorted[i].first * (double)m / (double)(i + 1);
        if (q < runMin) runMin = q;
        adj[sorted[i].second] = std::min(runMin, 1.0);
    }
    return adj;
}

// Parse GO integer IDs from tGoDbr entry (lines of "GO:NNNNNNN[ ...]")
static inline void parseGoIds(const char* data, size_t len, std::vector<int>& out) {
    const char* ptr = data;
    const char* end = data + len;
    while (ptr < end) {
        const char* nl = (const char*)memchr(ptr, '\n', end - ptr);
        size_t lineLen = nl ? (size_t)(nl - ptr) : (size_t)(end - ptr);
        if (lineLen >= 10 && ptr[0] == 'G' && ptr[1] == 'O' && ptr[2] == ':') {
            int id = 0;
            for (size_t j = 3; j < lineLen && ptr[j] >= '0' && ptr[j] <= '9'; ++j)
                id = id * 10 + (ptr[j] - '0');
            if (id > 0) out.push_back(id);
        }
        if (!nl) break;
        ptr = nl + 1;
    }
}

} // namespace GoEnrichment

#endif
