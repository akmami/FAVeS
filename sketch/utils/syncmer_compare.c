/*
 * syncmer_compare.c
 *
 * Reads a FASTA file and, for each record, computes closed (s,k)-syncmers two
 * independent ways, then checks that they agree exactly:
 *
 *   1. FAVeS     : the real circular-buffer sketch_syncmers() from syncmer.c
 *   2. Original  : the reference algorithm from the "syncmer" repository, i.e.
 *                  WordToKmer / GetSubkmer / GetMinSubkmerPos_Hash (kmer.cpp) and
 *                  the SyncmerIndex2::CalcIsSyncmer decision (syncmerindex2.cpp),
 *                  reproduced verbatim below. (The original's C++ class itself pulls in
 *                  OpenMP + its option/utils framework, so its core is inlined here
 *                  rather than linked.)
 *
 * Both are forward-strand only, use murmur64, and A/C/G/T = 0/1/2/3.
 *
 * Build (from FAVeS/sketch):
 *     cc -O2 -o syncmer_compare syncmer_compare.c syncmer.c
 *
 * Run:
 *     ./syncmer_compare <file.fasta> [k] [s]        (defaults k=15 s=5)
 *
 * Exit code 0 = every syncmer / hash / position matched (ignoring windows that
 * contain an ambiguous base, see note below); 1 = a real mismatch was found.
 *
 * Note on ambiguous bases (N etc.): FAVeS resets on a non-ACGT base and never
 * emits a k-mer spanning it. The reference's WordToKmer returns UINT64_MAX for
 * such a window and still indexes it, so the two legitimately differ there.
 * Mismatches whose k-mer window contains a non-ACGT base are reported
 * separately and do NOT fail the run.
 */

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "syncmer.h"

/* ------------------------------------------------------------------ */
/* Reference core, reproduced verbatim from the syncmer repository.   */
/* ------------------------------------------------------------------ */

/* murmur.h */
static inline uint64_t murmur64(uint64_t h) {
    h ^= (h >> 33);
    h *= 0xff51afd7ed558ccdULL;
    h ^= (h >> 33);
    h *= 0xc4ceb9fe1a85ec53ULL;
    h ^= (h >> 33);
    return h;
}

/* alpha.cpp: g_CharToLetterNucleo (A/a=0 C/c=1 G/g=2 T/t=3, else >3) */
static int char_to_letter(unsigned char c) {
    switch (c) {
        case 'A': case 'a': return 0;
        case 'C': case 'c': return 1;
        case 'G': case 'g': return 2;
        case 'T': case 't': return 3;
        default:            return 0xff;
    }
}

/* kmer.cpp: WordToKmer */
static uint64_t word_to_kmer(const char *Seq, unsigned k) {
    uint64_t Kmer = 0;
    for (unsigned i = 0; i < k; ++i) {
        int Letter = char_to_letter((unsigned char)Seq[i]);
        if (Letter > 3)
            return UINT64_MAX;
        Kmer = (Kmer << 2) | (unsigned)Letter;
    }
    return Kmer;
}

/* kmer.cpp: GetSubkmer */
static uint64_t get_subkmer(uint64_t Kmer, unsigned k, unsigned s, unsigned Pos) {
    uint64_t Mer = 0;
    Kmer >>= 2 * (k - Pos - s);
    for (unsigned j = 0; j < s; ++j) {
        unsigned Letter = (unsigned)((Kmer >> (2 * j)) & 3);
        Mer = Mer | ((uint64_t)Letter << (2 * j));
    }
    return Mer;
}

/* kmer.cpp: GetMinSubkmerPos_Hash (leftmost min on ties: strict '<', Pos ascending) */
static unsigned min_subkmer_pos_hash(uint64_t Kmer, unsigned k, unsigned s) {
    unsigned BestPos = 0;
    uint64_t BestMer = murmur64(get_subkmer(Kmer, k, s, 0));
    for (unsigned Pos = 1; Pos + s <= k; ++Pos) {
        uint64_t Mer = murmur64(get_subkmer(Kmer, k, s, Pos));
        if (Mer < BestMer) {
            BestMer = Mer;
            BestPos = Pos;
        }
    }
    return BestPos;
}

/* syncmerindex2.cpp: SyncmerIndex2::CalcIsSyncmer, closed variant, ds=0, d=0 */
static int is_syncmer(uint64_t Kmer, unsigned k, unsigned s) {
    unsigned MinPos = min_subkmer_pos_hash(Kmer, k, s);
    return (MinPos == 0 || MinPos == k - s);
}

/* ------------------------------------------------------------------ */
/* Minimal FASTA reader                                               */
/* ------------------------------------------------------------------ */

typedef struct {
    char  *name;
    char  *seq;
    size_t len;
} Record;

/* Returns number of records, fills *recs (caller frees). */
static size_t read_fasta(const char *path, Record **out_recs) {
    FILE *f = fopen(path, "r");
    if (!f) {
        fprintf(stderr, "error: cannot open %s\n", path);
        exit(2);
    }

    Record *recs = NULL;
    size_t  nrec = 0, cap = 0;
    char   *seqbuf = NULL;
    size_t  seqlen = 0, seqcap = 0;
    char   *name = NULL;

    char line[1 << 16];
    while (fgets(line, sizeof(line), f)) {
        if (line[0] == '>') {
            /* flush previous record */
            if (name) {
                if (nrec == cap) {
                    cap = cap ? cap * 2 : 8;
                    recs = (Record *)realloc(recs, cap * sizeof(Record));
                }
                char *s = (char *)malloc(seqlen + 1);
                memcpy(s, seqbuf, seqlen);
                s[seqlen] = 0;
                recs[nrec].name = name;
                recs[nrec].seq  = s;
                recs[nrec].len  = seqlen;
                ++nrec;
            }
            /* start new record: strip '>' and trailing newline */
            size_t n = strlen(line);
            while (n && (line[n - 1] == '\n' || line[n - 1] == '\r'))
                line[--n] = 0;
            name = strdup(line + 1);
            seqlen = 0;
        } else if (name) {
            for (size_t i = 0; line[i]; ++i) {
                char c = line[i];
                if (c == '\n' || c == '\r' || c == ' ' || c == '\t')
                    continue;
                if (seqlen == seqcap) {
                    seqcap = seqcap ? seqcap * 2 : 1024;
                    seqbuf = (char *)realloc(seqbuf, seqcap);
                }
                seqbuf[seqlen++] = c;
            }
        }
    }
    /* flush last */
    if (name) {
        if (nrec == cap) {
            cap = cap ? cap * 2 : 8;
            recs = (Record *)realloc(recs, cap * sizeof(Record));
        }
        char *s = (char *)malloc(seqlen + 1);
        memcpy(s, seqbuf, seqlen);
        s[seqlen] = 0;
        recs[nrec].name = name;
        recs[nrec].seq  = s;
        recs[nrec].len  = seqlen;
        ++nrec;
    }

    free(seqbuf);
    fclose(f);
    *out_recs = recs;
    return nrec;
}

static int window_has_ambiguous(const char *seq, size_t i, unsigned k) {
    for (unsigned j = 0; j < k; ++j)
        if (char_to_letter((unsigned char)seq[i + j]) > 3)
            return 1;
    return 0;
}

/* ------------------------------------------------------------------ */

int main(int argc, char **argv) {
    if (argc < 2) {
        fprintf(stderr, "usage: %s <file.fasta> [k] [s]\n", argv[0]);
        return 2;
    }
    const char *path = argv[1];
    unsigned k = (argc > 2) ? (unsigned)atoi(argv[2]) : 15;
    unsigned s = (argc > 3) ? (unsigned)atoi(argv[3]) : 5;

    /* Reference valid domain for a closed syncmer (from CalcIsSyncmer asserts). */
    if (!(s >= 2 && s + 2 <= k && k <= 32)) {
        fprintf(stderr, "error: need s>=2 and s+2<=k and k<=32 (got k=%u s=%u)\n", k, s);
        return 2;
    }

    Record *recs = NULL;
    size_t  nrec = read_fasta(path, &recs);
    printf("file=%s  records=%zu  k=%u  s=%u\n\n", path, nrec, k, s);

    long tot_original = 0, tot_faves = 0;
    long real_mm = 0, ambig_mm = 0;
    long meta_bad = 0;

    for (size_t r = 0; r < nrec; ++r) {
        const char *seq = recs[r].seq;
        size_t L = recs[r].len;
        if (L < k) {
            printf("[%zu] %-30.30s len=%zu  (skipped, shorter than k)\n", r, recs[r].name, L);
            continue;
        }
        size_t K = L - k + 1;

        /* --- original side (coord order, forward strand) --- */
        unsigned char *is = (unsigned char *)calloc(K, 1);
        uint64_t      *h  = (uint64_t *)calloc(K, sizeof(uint64_t));
        long noriginal = 0;
        for (size_t i = 0; i < K; ++i) {
            uint64_t kmer = word_to_kmer(seq + i, k); /* UINT64_MAX on ambiguous */
            if (is_syncmer(kmer, k, s)) {
                is[i] = 1;
                h[i]  = murmur64(kmer);
                ++noriginal;
            }
        }

        /* --- FAVeS side --- */
        uint128_t *out = NULL;
        uint64_t   n = sketch_syncmers(seq, (int)L, (int)s, (int)k, (uint32_t)r, &out);
        unsigned char *fav_is = (unsigned char *)calloc(K, 1);
        uint64_t      *fav_h   = (uint64_t *)calloc(K, sizeof(uint64_t));
        for (uint64_t i = 0; i < n; ++i) {
            uint32_t pos    = __syncmer_get_index(out[i]);
            uint32_t rid    = __syncmer_get_reference_id(out[i]);
            uint32_t strand = __syncmer_get_strand(out[i]);
            if (rid != (uint32_t)r || strand != 0) {
                ++meta_bad;
                fprintf(stderr, "[%zu] bad meta at pos=%u: rid=%u strand=%u\n", r, pos, rid, strand);
            }
            if (pos < (uint32_t)K) {
                fav_is[pos] = 1;
                fav_h[pos]  = out[i].x;
            }
        }

        /* --- Compare --- */
        long rec_real = 0, rec_ambig = 0;
        for (size_t i = 0; i < K; ++i) {
            int disagree_set  = (is[i] != fav_is[i]);
            int disagree_hash = (is[i] && fav_is[i] && h[i] != fav_h[i]);
            if (disagree_set || disagree_hash) {
                if (window_has_ambiguous(seq, i, k)) ++rec_ambig;
                else {
                    ++rec_real;
                    if (rec_real <= 5)
                        fprintf(stderr,
                            "[%zu] MISMATCH coord=%zu  original(is=%d h=%016llx)  faves(is=%d h=%016llx)\n",
                            r, i, is[i], (unsigned long long)h[i],
                            fav_is[i], (unsigned long long)fav_h[i]);
                }
            }
        }

        tot_original  += noriginal;
        tot_faves     += (long)n;
        real_mm       += rec_real;
        ambig_mm      += rec_ambig;

        printf("[%zu] %-30.30s len=%zu  kmers=%zu  original=%ld faves=%llu  real_mm=%ld ambig_mm=%ld\n",
               r, recs[r].name, L, K, noriginal, (unsigned long long)n, rec_real, rec_ambig);

        free(is); free(h); free(fav_is); free(fav_h); free(out);
        free(recs[r].name); free(recs[r].seq);
    }
    free(recs);

    printf("\n---- summary ----\n");
    printf("total syncmers : original=%ld  faves=%ld\n", tot_original, tot_faves);
    printf("real mismatches (ACGT windows)      : %ld\n", real_mm);
    printf("ambiguous-window mismatches (ignored): %ld\n", ambig_mm);
    printf("bad FAVeS metadata (rid/strand)     : %ld\n", meta_bad);
    printf("RESULT: %s\n", (real_mm == 0 && meta_bad == 0) ? "MATCH" : "MISMATCH");

    return (real_mm == 0 && meta_bad == 0) ? 0 : 1;
}
