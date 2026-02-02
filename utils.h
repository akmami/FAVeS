#ifndef __UTILS_H__
#define __UTILS_H__

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "struct_def.h"


static inline char rc_base(char b) {
    switch (b) {
        case 'A': return 'T';
        case 'C': return 'G';
        case 'G': return 'C';
        case 'T': return 'A';
        case 'a': return 'T';
        case 'c': return 'G';
        case 'g': return 'C';
        case 't': return 'A';
        default:  return 'N';
    }
}

static inline void reverse_complement(const char *in, char *out, int len) {
    for (int i = 0; i < len; i++) {
        out[len - i - 1] = rc_base(in[i]);
    }
}

void free_seqs(ref_seq_t **seqs, int seq_len);

int store_seqs(const char *path, ref_seq_t **seqs);

int banded_align_and_report(const char *ref, uint64_t ref_span, int ref_strand, const char *read, uint64_t read_span, int read_strand, uint64_t ref_pos, uint64_t ref_id, var_bvec_t *variants);
// -------------------------------------------
// --- THREADING HELPER FUNCS
// -------------------------------------------
void queue_init(job_queue_t *q, int capacity);

void queue_push(job_queue_t *q, batch_records_t *job);

batch_records_t *queue_pop(job_queue_t *q);

static inline batch_records_t *batch_records_create(int cap) {
    batch_records_t *b = malloc(sizeof(batch_records_t));
    b->records = (record_t *)malloc(sizeof(record_t) * cap);
    b->len = 0;
    return b;
}

static inline void batch_records_destroy(batch_records_t *b) {
    for (int i = 0; i < b->len; i++)
        free(b->records[i].bases);

    free(b->records);
    free(b);
}

#endif