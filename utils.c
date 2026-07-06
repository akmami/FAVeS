#include "utils.h"
#include <time.h>


#define abs_u64(x, y) (ref_span > read_span ? (ref_span - read_span) : (read_span - ref_span))

static __thread wavefront_aligner_t *wf_aligner_tls = NULL;

static wavefront_aligner_t *get_thread_aligner(void) {
    if (wf_aligner_tls) return wf_aligner_tls;
    wavefront_aligner_attr_t attr = wavefront_aligner_attr_default;
    
    attr.distance_metric = gap_linear;

    attr.linear_penalties.match = 0;
    attr.linear_penalties.mismatch = MISMATCH;
    attr.linear_penalties.indel = GAP;

    attr.alignment_scope = compute_alignment;
    attr.alignment_form.span = alignment_endsfree;
    attr.alignment_form.pattern_begin_free = BAND;
    attr.alignment_form.pattern_end_free = BAND;
    attr.alignment_form.text_begin_free = BAND;
    attr.alignment_form.text_end_free = BAND;

    attr.heuristic.strategy = wf_heuristic_banded_static;
    attr.heuristic.min_k = -BAND;
    attr.heuristic.max_k = +BAND;

    wf_aligner_tls = wavefront_aligner_new(&attr);
    return wf_aligner_tls;
}

void release_thread_aligner(void) {
    if (wf_aligner_tls) {
        wavefront_aligner_delete(wf_aligner_tls);
        wf_aligner_tls = NULL;
    }
}

static unsigned char seq_nt4_table[256] = {
	0, 1, 2, 3,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  3, 3, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  3, 3, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4
};

void free_seqs(ref_seq_t **seqs, int seq_len) {
    if (seq_len) {
        ref_seq_t *temp = *seqs;
        for (int i = 0; i < seq_len; i++) {
            if (temp[i].len) {
                free(temp[i].chrom);
            }
            if (temp[i].header) free(temp[i].header);
        }

        free(temp);
    }
}

int store_seqs(const char *path, ref_seq_t **seqs) {

    char *fai_path = (char *)malloc(strlen(path) + 5);
    if (fai_path == NULL) {
        fprintf(stderr, "[%s::err] Memory allocation failed\n", __TOOL_SHORT_NAME__);
        exit(-1);
    }
    sprintf(fai_path, "%s.fai", path);

    size_t line_cap = 1024;
    char *line = (char *)malloc(line_cap);
    int chrom_index = 0;

    FILE *fai = fopen(fai_path, "r");
    if (!fai) {
        fprintf(stderr, "[%s::err] Index file %s does not exist\n", __TOOL_SHORT_NAME__, fai_path);
        fprintf(stderr, "           samtools faidx %s\n", path);
        exit(-1);
    }

    while (fgets(line, line_cap, fai)) {
        chrom_index++;
    }

    if (chrom_index == 0) {
        fprintf(stderr, "[%s::err] Index file is empty.\n", __TOOL_SHORT_NAME__);
        exit(-1);
    }
    
    *seqs = (ref_seq_t *)malloc(sizeof(ref_seq_t) * chrom_index);
    if (!(*seqs)) {
        fprintf(stderr, "[%s::err] Couldn't allocate memory to ref sequences\n", __TOOL_SHORT_NAME__);
        exit(-1);
    }

    rewind(fai);
    chrom_index = 0;

    while (fgets(line, line_cap, fai)) {
        char *name, *length;

        // assign name
        char *saveptr;
        name = strtok_r(line, "\t", &saveptr);
        uint64_t name_len = strlen(name);
        (*seqs)[chrom_index].header = (char *)malloc(name_len+1);
        memcpy((*seqs)[chrom_index].header, name, name_len);
        (*seqs)[chrom_index].header[name_len] = '\0';

        // assign size and allocate in memory
        length = strtok_r(NULL, "\t", &saveptr);
        (*seqs)[chrom_index].len = strtol(length, NULL, 10);
        chrom_index++;
    }

    fclose(fai);
    free(fai_path);
    free(line);

    return chrom_index;
}

int banded_align_and_report(const char *ref, uint64_t ref_span, const char *read, uint64_t read_span, int read_strand, uint64_t ref_pos, uint64_t ref_id, var_bvec_t *variants) {

    if (abs_u64(ref_span, read_span) > BAND) return 0; // why bother?

    char read_buf[MAX_LEN] = {0};
    if (read_strand) reverse_complement(read, read_buf, read_span);
    else memcpy(read_buf, read, read_span);
    for (uint64_t k = 0; k < read_span; k++) read_buf[k] &= to_uppercase_mask;

    wavefront_aligner_t *wf = get_thread_aligner();
    int status = wavefront_align(wf, ref, (int)ref_span, read_buf, (int)read_span);
    if (status != WF_STATUS_ALG_COMPLETED) return 0;

    if (wf->align_status.score > (int)(4 * MAX_U64(2, CEIL_PERCENT_U64(MAX_U64(read_span, ref_span), WFA_RATE_PERCENT)) + 3)) {
        return 0;
    }

    char *ops = wf->cigar->operations;
    int end   = wf->cigar->end_offset;

    int pattern_pos = 0, text_pos = 0;
    int mismatches = 0;
    for (int i = wf->cigar->begin_offset; i < end; i++) {
        char op = ops[i];
        if (op == 'M') {
            pattern_pos++; text_pos++;
        } else if (op == 'X') {
            mismatches++;
            if (pattern_pos != 0 && text_pos != 0 &&
                pattern_pos != (int)ref_span - 1 && text_pos != (int)read_span - 1) {
                var_add(variants, __set_variation_bits(ref_id, ref_pos + pattern_pos, seq_nt4_table[(int)read_buf[text_pos]]));
            }
            pattern_pos++; text_pos++;
        } else if (op == 'I') {
            text_pos++;
        } else if (op == 'D') {
            pattern_pos++;
        }
    }

    return mismatches;
}

void queue_init(job_queue_t *q, int capacity) {
    q->buf = malloc(sizeof(batch_records_t *) * capacity);
    q->capacity = capacity;
    q->head = q->tail = q->count = 0;
    q->done = 0;

    pthread_mutex_init(&q->lock, NULL);
    pthread_cond_init(&q->not_empty, NULL);
    pthread_cond_init(&q->not_full, NULL);
}

void queue_push(job_queue_t *q, batch_records_t *records) {
    pthread_mutex_lock(&q->lock);

    while (q->count == q->capacity)
        pthread_cond_wait(&q->not_full, &q->lock);

    q->buf[q->tail] = records;
    q->tail = (q->tail + 1) % q->capacity;
    q->count++;

    pthread_cond_signal(&q->not_empty);
    pthread_mutex_unlock(&q->lock);
}

batch_records_t *queue_pop(job_queue_t *q) {
    pthread_mutex_lock(&q->lock);

    while (q->count == 0 && !q->done)
        pthread_cond_wait(&q->not_empty, &q->lock);

    if (q->count == 0 && q->done) {
        pthread_mutex_unlock(&q->lock);
        return NULL;
    }

    batch_records_t *records = q->buf[q->head];
    q->head = (q->head + 1) % q->capacity;
    q->count--;

    pthread_cond_signal(&q->not_full);
    pthread_mutex_unlock(&q->lock);

    return records;
}