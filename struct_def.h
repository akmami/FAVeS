#ifndef __STRUCT_DEF_H__
#define __STRUCT_DEF_H__

#include <stdlib.h>
#include <stdint.h>
#include "pthread.h"

#define BLEND_K_SHORT 21
#define BLEND_W_SHORT 11
#define BLEND_BITS_SHORT 32
#define BLEND_NEIGHBOR_NUMBER_SHORT 5
#define BLEND_K_HIFI 19
#define BLEND_W_HIFI 50
#define BLEND_BITS_HIFI 38
#define BLEND_NEIGHBOR_NUMBER_HIFI 5
#define THREAD_NUMBER 4
#define SKETCH_CAPACITY 40000000
#define DISTANCE_THRESHOLD 200
#define QUEUE_SIZE 128
#define BATCH_SIZE 1024

#define READ_COUNT_BREAKPOINT

#define PROGRESS 0

#define MATCH     2
#define MISMATCH -5
#define GAP      -3
#define BAND 20        // ±8 bp band
#define MAX_LEN 2048   // max extension length
#define NEG_INF -100000000
#define MIN_SCORE 64

#define abs_diff(x, y) (x < y ? y - x : x - y)


typedef struct {
    char *header;
    char *chrom;
    uint64_t len;
} ref_seq;

typedef struct {
    char *fasta;
    char *fastq;
    int k;
    int w;
    int blend_bits;
    int n_neighbors;
    int n_threads;
    int progress;
} params;

typedef struct {
    uint64_t found;
    uint64_t not_found;
    uint64_t mismatches;
    uint64_t total_len;
} stats_t;

typedef struct {
    int k;
    int w;
    int blend_bits;
    int n_neighbors;
    ref_seq *seqs;
    void *fuzzy_seeds;
    uint64_t unique_fuzzy_seeds_len;
    uint64_t relaxed_fuzzy_seeds_len;
    void *index_table;
    stats_t stats;
} params_lite;

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
    params_lite *p;
} worker_ctx_t;

#endif