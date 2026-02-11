#ifndef __STRUCT_DEF_H__
#define __STRUCT_DEF_H__

#include <stdlib.h>
#include <stdint.h>
#include "pthread.h"
#include "utils/bvec.h"

#define __DEFAULT_BLEND_K__ 21
#define __DEFAULT_BLEND_w__ 11
#define __DEFAULT_BLEND_BITS__ 32
#define __DEFAULT_BLEND_NEIGHBOR_NUMBER__ 5
#define __DEFAULT_CONSENSUS_THRESHOLD__ 10
#define __DEFAULT_THREAD_NUMBER__ 4
#define __DEFAULT_SKETCH_CAPACITY__ 40000000
#define __DEFAULT_VARIANT_CAPACITY__ 10000000
#define __DEFAULT_DISTANCE_THRESHOLD__ 200
#define __DEFAULT_QUEUE_SIZE__ 1024
#define __DEFAULT_BATCH_SIZE__ 512
#define __DEFAULT_PROGRESS_INTERVAL__ 1

#define __DEFAULT_PROGRESS__ 0

#define MATCH       2
#define MISMATCH    -5
#define GAP         -3
#define BAND        20     // +-10 bp band
#define MAX_LEN     512   // max extension length
#define NEG_INF     -100000000
#define MIN_SCORE   64

#define abs_diff(x, y) ((x) < (y) ? (y) - (x) : (x) - (y))

#define to_uppercase_mask 0xDF
 

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
    char *fastq;
    char *bed;
    int k;
    int w;
    int blend_bits;
    int n_neighbors;
    int min_consensus;
    int n_threads;
    int progress;
} params_t;

typedef struct {
    uint64_t countains_unique;  // read-wise
    uint64_t only_non_unique;   // read-wise
    uint64_t no_seed;           // read-wise
    uint64_t total_seed_count;  // overall
    uint64_t mismatch_count;    // overall
    uint64_t neighbour_count;   // blend span
} stats_t;

// -----------------------------------------------------------
// -----------------------------------------------------------
// VARIATION RELATED STRUCTURES
// -----------------------------------------------------------
// -----------------------------------------------------------
static const uint64_t variation_index = ((1ULL << 52) - 1) << 3;

#define __get_variation_chrom(x) ((x) >> 55)
#define __get_variation_index(x) (((x) & variation_index) >> 3) 
#define __get_variation_type(x) ((x) & 7)

#define __set_variation_chrom(x, idx) ((x) |= ((uint64_t)(idx) << 55))
#define __set_variation_index(x, idx) ((x) |= ((uint64_t)(idx) << 3)) 
#define __set_variation_type(x, type) ((x) |= ((type) & 7))
#define __set_variation_bits(chrom, index, type) (((chrom) << 55) | ((index) << 3) | ((type) & 7))
#define __set_variation(x, chrom, index, type) ((x) = __set_variation_bits(chrom, index, type))

typedef uint64_t variation_t; // chrom (9 bits) + index (52 bits) + type (3 bits)

BVEC_INIT(var, variation_t)

// -----------------------------------------------------------
// -----------------------------------------------------------
// THREAD RELATED STRUCTURES
// -----------------------------------------------------------
// -----------------------------------------------------------
typedef struct {
    int k;
    int w;
    int blend_bits;
    int n_neighbors;
    int n_threads;
    ref_seq_t *seqs;
    void *fuzzy_seeds;
    uint64_t fuzzy_seeds_len;
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
    int  len;
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