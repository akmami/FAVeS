// ============================================================================
// match-length.cpp
//
// Exact k-mer uniqueness analysis via the maximal-match-length distribution.
//
// For every genomic position i we compute L_i, the length of the longest
// substring starting at i that also occurs elsewhere in the genome (its
// "maximal match length"), using a suffix array + Longest Common Prefix (libdivsufsort).
//
// The k-mer at i is unique  <=>  k > L_i  (its length-k prefix occurs nowhere
// else). Hence the fraction of unique k-mers is exactly the CDF of L:
//
//     p(k) = Pr(L_i < k | position i can hold a length-k k-mer).
//
// One suffix-array pass therefore yields the ENTIRE uniqueness-vs-k curve for
// all k at once, plus the raw match-length (repeat-length) distribution and the
// base composition needed for the null model. This backs analyses 1-5 and 7:
//   (1) direct check of the p(k) = CDF(L) identity,
//   (2) crossover k_c = ln N / ln(1/Q),
//   (3) plateau = repeat content (fit 1 - a*S_rep(k) to the tail),
//   (4) cross-species comparison (run per genome),
//   (5) composition control (Q from A/C/G/T counts),
//   (7) repeat-length resolution (the match-length histogram itself).
//
// Single strand, forward orientation -- matching sketch/kmer.c, so p(15)/p(21)
// here reproduce the KMER path's unique_seeds_ratio.
//
// ============================================================================
#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <cstring>
#include <cmath>
#include <string>
#include <vector>
#include <fstream>
#include <iostream>
#include <algorithm>

#include <divsufsort.h>
#include <divsufsort64.h>


static inline uint8_t base_code(char ch) {
    switch (ch) {
        case 'A': case 'a': return 1;
        case 'C': case 'c': return 2;
        case 'G': case 'g': return 3;
        case 'T': case 't': return 4;
        default:            return 0;
    }
}

struct Result {
    uint64_t genome_len = 0;                // total residues seen (incl. N)
    uint64_t valid_len  = 0;                // A/C/G/T residues
    uint64_t base_counts[5] = {0,0,0,0,0};  // idx 1..4 used
    std::vector<uint64_t> cap_hist;         // cap_hist[c] = # positions with capacity exactly c
    std::vector<uint64_t> diff_unique;      // difference array over k for unique(k)
    std::vector<uint64_t> L_hist;           // capped histogram of L_i (for repeat lengths)
    uint64_t L_hist_overflow = 0;           // positions with L_i >= Lhist_cap
    uint64_t max_cap = 0;
};

// Separators are inserted between records and wherever a non-ACGT residue sits, so no match can cross a chromosome boundary or an N run.
static void read_fasta(const char *path, std::vector<uint8_t> &text, Result &R) {
    std::ifstream in(path);
    if (!in.good()) { std::cerr << "Error opening: " << path << "\n"; exit(1); }
    text.reserve(1u << 20);
    std::string line;
    bool have_record = false;
    while (std::getline(in, line)) {
        if (!line.empty() && line[0] == '>') {
            if (have_record) text.push_back(0); // break between chromosomes
            have_record = true;
            continue;
        }
        for (char ch : line) {
            if (ch == '\r') continue;
            R.genome_len++;
            uint8_t c = base_code(ch);
            text.push_back(c);            // 0 for N/other == break
            if (c) { R.valid_len++; R.base_counts[c]++; }
        }
    }
    text.push_back(0); // final sentinel
}

// capacity array: cap[i] = # consecutive valid bases starting at i
static void build_capacity(const std::vector<uint8_t> &text, std::vector<uint32_t> &cap) {
    size_t n = text.size();
    cap.assign(n, 0);
    uint32_t run = 0;
    for (size_t idx = n; idx-- > 0; ) {
        if (text[idx] == 0) run = 0;
        else                run = (run == UINT32_MAX ? run : run + 1);
        cap[idx] = run;
    }
}

// templated core: build SA/LCP with idx_t, accumulate histograms
template <typename idx_t>
static int core(const std::vector<uint8_t> &text, const std::vector<uint32_t> &cap,
                Result &R, int (*sufsort)(const sauchar_t*, idx_t*, idx_t),
                uint32_t Lhist_cap) {
    idx_t n = (idx_t)text.size();
    std::vector<idx_t> SA(n);
    if (sufsort(text.data(), SA.data(), n) != 0) {
        std::cerr << "suffix sort failed\n"; return 1;
    }
    // rank = inverse SA
    std::vector<idx_t> rank(n);
    for (idx_t r = 0; r < n; r++) rank[SA[r]] = r;

    // Kasai LCP: lcp[r] = lcp(SA[r-1], SA[r]); lcp[0] = 0.
    std::vector<idx_t> lcp(n, 0);
    idx_t h = 0;
    for (idx_t i = 0; i < n; i++) {
        if (rank[i] > 0) {
            idx_t j = SA[rank[i] - 1];
            while (i + h < n && j + h < n && text[i + h] == text[j + h]) h++;
            lcp[rank[i]] = h;
            if (h > 0) h--;
        } else {
            h = 0;
        }
    }

    // Per position: L_i = max(lcp[rank_i], lcp[rank_i + 1]); cap to segment end.
    R.cap_hist.assign(R.max_cap + 2, 0);
    R.diff_unique.assign(R.max_cap + 2, 0);
    R.L_hist.assign(Lhist_cap, 0);
    for (idx_t r = 0; r < n; r++) {
        idx_t pos = SA[r];
        uint32_t c = cap[pos];
        if (c == 0) continue;                 // separator position, no k-mer
        idx_t left  = (r > 0)     ? lcp[r]     : 0;
        idx_t right = (r + 1 < n) ? lcp[r + 1] : 0;
        idx_t rawL  = std::max(left, right);
        uint32_t L  = (uint32_t)std::min<idx_t>(rawL, (idx_t)c); // cap at segment
        // total(k): this position holds a k-mer for k in [1, c]
        R.cap_hist[c]++;
        // unique(k): unique for k in [L+1, c]
        if (L < c) {
            R.diff_unique[L + 1]++;
            R.diff_unique[c + 1]--;
        }
        // repeat-length (match-length) histogram
        if (L < Lhist_cap) R.L_hist[L]++; else R.L_hist_overflow++;
    }
    return 0;
}

int main(int argc, char **argv) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0]
                  << " <genome.fa> [shortname] [kmax=512] [Lhist_cap=2048]\n";
        return 1;
    }
    const char *path      = argv[1];
    const char *shortname = (argc > 2) ? argv[2] : "";
    uint32_t kmax      = (argc > 3) ? (uint32_t)atoi(argv[3]) : 512;
    uint32_t Lhist_cap = (argc > 4) ? (uint32_t)atoi(argv[4]) : 2048;

    Result R;
    std::vector<uint8_t> text;
    read_fasta(path, text, R);

    std::vector<uint32_t> cap;
    build_capacity(text, cap);
    for (uint32_t c : cap) if (c > R.max_cap) R.max_cap = c;

    // choose 32- vs 64-bit index by genome size
    int rc;
    if (text.size() < (size_t)INT32_MAX) {
        rc = core<saidx_t>(text, cap, R, &divsufsort, Lhist_cap);
    } else {
        rc = core<saidx64_t>(text, cap, R, &divsufsort64, Lhist_cap);
    }
    if (rc) return rc;

    // derive p(k) curve
    // total(k)  = sum_{c>=k} cap_hist[c]           (suffix sum)
    // unique(k) = prefix sum of diff_unique up to k
    uint32_t K = std::min<uint32_t>(kmax, R.max_cap);
    std::vector<uint64_t> total(K + 2, 0), uniq(K + 2, 0);
    // suffix sum for total
    uint64_t running = 0;
    for (uint32_t c = R.max_cap; c >= 1; c--) {
        running += R.cap_hist[c];
        if (c <= K) total[c] = running;
    }
    // prefix sum for unique
    uint64_t acc = 0;
    for (uint32_t k = 1; k <= K; k++) { acc += R.diff_unique[k]; uniq[k] = acc; }

    // composition / null-model constants
    double N = (double)R.valid_len;
    double fA = R.base_counts[1] / N, fC = R.base_counts[2] / N;
    double fG = R.base_counts[3] / N, fT = R.base_counts[4] / N;
    double Q  = fA*fA + fC*fC + fG*fG + fT*fT;
    double gc = fC + fG;
    double kc = std::log(N) / std::log(1.0 / Q);

    // emit JSON
    std::cout.setf(std::ios::fixed);
    std::cout << "{\n";
    std::cout << "\t\"genome\": \"" << path << "\",\n";
    std::cout << "\t\"short\": \"" << shortname << "\",\n";
    std::cout << "\t\"genome_len\": " << R.genome_len << ",\n";
    std::cout << "\t\"valid_len\": "  << R.valid_len  << ",\n";
    std::cout << "\t\"base_counts\": {\"A\": " << R.base_counts[1]
              << ", \"C\": " << R.base_counts[2]
              << ", \"G\": " << R.base_counts[3]
              << ", \"T\": " << R.base_counts[4] << "},\n";
    std::cout.precision(6);
    std::cout << "\t\"gc\": " << gc << ",\n";
    std::cout << "\t\"Q\": " << Q << ",\n";
    std::cout << "\t\"sigma_eff\": " << (1.0 / Q) << ",\n";
    std::cout << "\t\"k_c\": " << kc << ",\n";
    std::cout << "\t\"kmax\": " << K << ",\n";

    // uniqueness curve
    std::cout << "\t\"uniqueness\": [\n";
    for (uint32_t k = 1; k <= K; k++) {
        double p = total[k] ? (double)uniq[k] / (double)total[k] : 0.0;
        std::cout << "\t\t{\"k\": " << k
                  << ", \"total\": " << total[k]
                  << ", \"unique\": " << uniq[k]
                  << ", \"p_unique\": " << p << "}"
                  << (k == K ? "\n" : ",\n");
    }
    std::cout << "\t],\n";

    // match-length (repeat-length) histogram
    std::cout << "\t\"match_length_hist\": [\n";
    uint32_t last = Lhist_cap;
    while (last > 0 && R.L_hist[last - 1] == 0) last--;
    for (uint32_t L = 0; L < last; L++) {
        std::cout << "\t\t{\"L\": " << L << ", \"count\": " << R.L_hist[L] << "}"
                  << (L + 1 == last ? "\n" : ",\n");
    }
    std::cout << "\t],\n";
    std::cout << "\t\"match_length_overflow\": " << R.L_hist_overflow << "\n";
    std::cout << "}" << std::endl;
    return 0;
}
