#include "utils.h"


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
        fprintf(stderr, "[ERROR] Memory allocation failed\n");
        exit(-1);
    }
    sprintf(fai_path, "%s.fai", path);

    size_t line_cap = 1024;
    char *line = (char *)malloc(line_cap);
    int chrom_index = 0;

    FILE *fai = fopen(fai_path, "r");
    if (!fai) {
        fprintf(stderr, "[ERROR] Index file %s does not exist\n", fai_path);
        fprintf(stderr, "           samtools faidx %s\n", path);
        exit(-1);
    }

    while (fgets(line, line_cap, fai)) {
        chrom_index++;
    }

    if (chrom_index == 0) {
        fprintf(stderr, "[ERROR] Index file is empty.\n");
        exit(-1);
    }
    
    *seqs = (ref_seq_t *)malloc(sizeof(ref_seq_t) * chrom_index);
    if (!(*seqs)) {
        fprintf(stderr, "[ERROR] Couldn't allocate memory to ref sequences\n");
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
    
    int ref_len = ref_span;
    int read_len = read_span;

    if (abs(ref_len - read_len) > BAND) return 0; // why bother?

    char ref_buf[MAX_LEN];
    char read_buf[MAX_LEN];

    // Strand normalization
    memcpy(ref_buf, ref, ref_len);
    
    if (read_strand) reverse_complement(read, read_buf, read_len);
    else memcpy(read_buf, read, read_len);

    int dp[MAX_LEN + 1][MAX_LEN + 1];
    int bt[MAX_LEN + 1][MAX_LEN + 1];

    // init
    for (int i = 0; i <= ref_len; i++)
        for (int j = 0; j <= read_len; j++)
            dp[i][j] = NEG_INF;

    // free leading gaps within band
    for (int j = 0; j <= BAND && j <= read_len; j++)
        dp[0][j] = 0;

    for (int i = 0; i <= BAND && i <= ref_len; i++)
        dp[i][0] = 0;

    // DP
    for (int i = 1; i <= ref_len; i++) {

        int j_start = (i - BAND > 1) ? i - BAND : 1;
        int j_end   = (i + BAND < read_len) ? i + BAND : read_len;

        for (int j = j_start; j <= j_end; j++) {

            int best = NEG_INF, op = -1;

            // diagonal
            int diag = dp[i-1][j-1] + ((ref_buf[i-1] & to_uppercase_mask) == (read_buf[j-1] & to_uppercase_mask) ? MATCH : MISMATCH);
            best = diag; 
            op = 0; // diagonal movement

            int del = dp[i-1][j] + GAP;
            if (del > best) {
                best = del;
                op = 1; // downward movement
            }

            // insertion
            int ins = dp[i][j-1] + GAP;
            if (ins > best) {
                best = ins;
                op = 2; // right movement
            }

            dp[i][j] = best;
            bt[i][j] = op;
        }
    }

    // traceback
    int best = NEG_INF;
    int bi = ref_len, bj = read_len;

    // last row
    for (int j = 0; j <= read_len; j++) {
        if (dp[ref_len][j] > best) {
            best = dp[ref_len][j];
            bi = ref_len;
            bj = j;
        }
    }

    // last column
    for (int i = 0; i <= ref_len; i++) {
        if (dp[i][read_len] > best) {
            best = dp[i][read_len];
            bi = i;
            bj = read_len;
        }
    }

    if (best < MIN_SCORE) return 0;

    int i = bi, j = bj;

    int mismatches = 0;

    while (i > 0 && j > 0) {
        int op = bt[i][j];

        if (op == 0) {
            if ((ref_buf[i-1] & to_uppercase_mask) != (read_buf[j-1] & to_uppercase_mask)) {
                mismatches++;
                if (i-1 != 0 && j-1 != 0 && i != ref_len && j != read_len) {
                    var_add(variants, __set_variation_bits(ref_id, ref_pos + i - 1, seq_nt4_table[(int)read_buf[j-1]]));
                }
            }
            i--; j--;
        }
        else if (op == 1) {
            if (i-1 != 0 && j-1 != 0 && i != ref_len && j != read_len) {
                // do nothing for indels at this point
            }
            i--;
        }
        else {
            if (i-1 != 0 && j-1 != 0 && i != ref_len && j != read_len) {
                // do nothing for indels at this point
            }
            j--;
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