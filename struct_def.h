#ifndef __STRUCT_DEF_H__
#define __STRUCT_DEF_H__

#include <stdlib.h>
#include <stdint.h>
#include "pthread.h"
#include "utils/bvec.h"

#define __TOOL_SHORT_NAME__ "fvs"

#define __DEFAULT_BLEND_K__                     21
#define __DEFAULT_BLEND_w__                     11
#define __DEFAULT_BLEND_B__                     50
#define __DEFAULT_BLEND_N__                     3
#define __DEFAULT_SYNCMER_S__                   10
#define __DEFAULT_LCP_L__                       3
#define __DEFAULT_LCP_D__                       1
#define __DEFAULT_CONSENSUS_THRESHOLD__         10
#define __DEFAULT_CONSENSUS_FRAC_THRESHOLD__    0.75
#define __DEFAULT_RADIUS__                      4
#define __DEFAULT_THREAD_NUMBER__               4
#define __DEFAULT_SKETCH_CAPACITY__             40000000
#define __DEFAULT_VARIANT_CAPACITY__            10000000
#define __DEFAULT_DISTANCE_THRESHOLD__          200
#define __DEFAULT_QUEUE_SIZE__                  1024
#define __DEFAULT_BATCH_SIZE__                  512
#define __DEFAULT_PROGRESS_INTERVAL__           1

#define __DEFAULT_UPPER_BOUND_FREQUENCY__       2

#define __DEFAULT_PROGRESS__    0
#define __DEFAULT_VERBOSE__     0

// #define MATCH       2
#define MISMATCH    4
#define GAP         6
#define BAND        20     // +-20 bp band
#define MAX_LEN     512   // max extension length
#define MIN_SCORE   128
#define WFA_RATE_PERCENT 5ULL

#define abs_diff(x, y) ((x) < (y) ? (y) - (x) : (x) - (y))

#define to_uppercase_mask 0xDF


// -----------------------------------------------------------
// -----------------------------------------------------------
// SKETCHING RELATED DEFINITIONS
// -----------------------------------------------------------
// -----------------------------------------------------------
//
// Every sketching method emits its seeds into the same 128-bit record. The
// low 64 bits (`x`) pack the seed key together with its span; the split point
// between key and span differs per method and is captured at runtime by
// `kmer_shift` / `span_mask`. The high 64 bits (`y`) are laid out identically
// for all methods:  reference_id(32) | index(31) | strand(1).
//
// Because the accessors below take the layout by value (a plain shift and a
// plain mask), a runtime-selected shift compiles to the exact same single
// instruction a hard-coded constant would, so BLEND / MINIMIZER / SYNCMER stay
// as fast as the original BLEND-only build. LCP does not natively produce this
// record, so its cores are packed into the canonical layout once, at seed
// generation time (see sketch_lcp_seeds in faves.c).

#ifndef __UINT128_T__
#define __UINT128_T__
typedef struct {
    uint64_t x; // seed key + span
    uint64_t y; // reference_id(32) + index(31) + strand(1)
} uint128_t;
#endif

typedef enum {
    SKETCH_BLEND = 0,
    SKETCH_MINIMIZER,
    SKETCH_SYNCMER,
    SKETCH_LCP
} sketch_type_t;

typedef struct {
    sketch_type_t type;
    // per-method parameters
    int kmer_size;      // blend / minimizer / syncmer
    int window_size;    // blend / minimizer
    int blend_bits;     // blend
    int n_neighbors;    // blend
    int smer_size;      // syncmer
    int lcp_level;      // lcp
    int dct_count;      // lcp
    // layout derived from `type` (see sketch_args_finalize)
    int kmer_shift;
    uint64_t span_mask;
} sketch_args_t;

// Derive the (kmer_shift, span_mask) split from the selected method so that the accessors decode exactly what each sketch library encoded.
static inline void sketch_args_finalize(sketch_args_t *s) {
    switch (s->type) {
        case SKETCH_MINIMIZER:
        case SKETCH_SYNCMER:
            s->kmer_shift = 8;  s->span_mask = 0xFFULL;   break; // key(56) | span(8)
        case SKETCH_BLEND:
        case SKETCH_LCP:
        default:
            s->kmer_shift = 14; s->span_mask = 0x3FFFULL; break; // key(50) | span(14)
    }
}

// Seed field accessors. 
// `x`-derived fields take the runtime layout; 
// `y`-derived fields are identical for every method.
static inline uint64_t sketch_get_kmer(uint128_t s, const sketch_args_t *a)   { return s.x >> a->kmer_shift; }
static inline uint64_t sketch_get_length(uint128_t s, const sketch_args_t *a) { return s.x & a->span_mask; }
static inline uint64_t sketch_get_reference_id(uint128_t s)                   { return s.y >> 32; }
static inline uint64_t sketch_get_index(uint128_t s)                          { return (s.y & 0xFFFFFFFFULL) >> 1; }
static inline int      sketch_get_strand(uint128_t s)                         { return (int)(s.y & 1ULL); }

// -----------------------------------------------------------
// -----------------------------------------------------------
// MAIN PROGRAM RELATED STRUCTURES
// -----------------------------------------------------------
// -----------------------------------------------------------
typedef struct {
    char *header;
    char *chrom;
    uint64_t len;
} ref_seq_t;

typedef struct {
    char *fasta;
    char **fastq;
    int n_fastq;
    char *bed;
    int min_consensus;
    float min_consensus_frac;
    int use_consensus_frac;
    int radius;
    int n_threads;
    int progress;
    int verbose;
    uint64_t all_seq_len;       // length of all seqeunces combined
    sketch_args_t sketch;       // selected sketching method + its parameters
} params_t;

typedef struct {
    uint64_t total_uniq_count;  // per read
    uint64_t countains_unique;  // read-wise
    uint64_t only_non_unique;   // read-wise
    uint64_t no_seed;           // read-wise
    uint64_t total_seed_count;  // overall
    uint64_t mismatch_count;    // overall
    uint64_t neighbour_count;   // span
    uint64_t all_seq_len;       // length of all seqeunces combined
} stats_t;

// -----------------------------------------------------------
// -----------------------------------------------------------
// VARIATION RELATED STRUCTURES
// -----------------------------------------------------------
// -----------------------------------------------------------
static const uint64_t variation_index_mask = ((1ULL << 47) - 1) << 3;

#define __get_variation_chrom(x) ((x) >> 50)
#define __get_variation_index(x) (((x) & variation_index_mask) >> 3) 
#define __get_variation_type(x) ((x) & 7)

#define __set_variation_chrom(x, idx) ((x) |= ((uint64_t)(idx) << 50))
#define __set_variation_index(x, idx) ((x) |= ((uint64_t)(idx) << 3)) 
#define __set_variation_type(x, type) ((x) |= ((type) & 7))
#define __set_variation_bits(chrom, index, type) (((chrom) << 50) | ((index) << 3) | ((type) & 7))
#define __set_variation(x, chrom, index, type) ((x) = __set_variation_bits(chrom, index, type))

typedef uint64_t variation_t; // chrom (14 bits) + index (47 bits) + type (3 bits)

BVEC_INIT(var, variation_t)

// -----------------------------------------------------------
// -----------------------------------------------------------
// THREAD RELATED STRUCTURES
// -----------------------------------------------------------
// -----------------------------------------------------------
typedef struct {
    sketch_args_t sketch;     // selected sketching method + its parameters
    int radius;
    int n_threads;
    ref_seq_t *seqs;
    void *seeds;
    uint64_t seeds_len;
    void *index_table;
    stats_t stats;
    var_bvec_t *fv_variants;
} params_lite_t;

typedef struct {
    char *bases;
    int len;
} record_t;

typedef struct {
    record_t *records;
    int len;
} batch_records_t;

typedef struct {
    batch_records_t **buf;
    int capacity;
    int head, tail;
    int count;
    int done; // producer finished
    pthread_mutex_t lock;
    pthread_cond_t not_empty;
    pthread_cond_t not_full;
} job_queue_t;

typedef struct {
    job_queue_t *queue;
    params_lite_t *p;
} worker_ctx_t;


#endif