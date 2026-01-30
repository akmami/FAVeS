#include "utils.h"


void free_seqs(ref_seq **seqs, int seq_len) {
    if (seq_len) {
        ref_seq *temp = *seqs;
        for (int i = 0; i < seq_len; i++) {
            if (temp[i].len) {
                free(temp[i].chrom);
            }
            if (temp[i].header) free(temp[i].header);
        }

        free(temp);
    }
}

int store_seqs(const char *path, ref_seq **seqs) {

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
        fprintf(stderr, "[ERROR] Couldn't open %s\n", fai_path);
        exit(-1);
    }

    while (fgets(line, line_cap, fai)) {
        chrom_index++;
    }

    if (chrom_index == 0) {
        fprintf(stderr, "[ERROR] Index file is empty.\n");
        exit(-1);
    }
    
    *seqs = (ref_seq *)malloc(sizeof(ref_seq) * chrom_index);
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

int banded_align_and_report(const char *ref, uint64_t ref_span, const char *read, uint64_t read_span, uint64_t ref_pos, uint64_t ref_id, int ref_strand, int read_strand, uint64_t read_pos, int len) {
    
    int ref_len = ref_span;
    int read_len = read_span;

    if (abs(ref_len - read_len) > BAND) return 0; // why bother?

    char ref_buf[MAX_LEN];
    char read_buf[MAX_LEN];

    // Strand normalization
    if (ref_strand) reverse_complement(ref, ref_buf, ref_len);
    else memcpy(ref_buf, ref, ref_len);
    
    if (read_strand) reverse_complement(read, read_buf, read_len);
    else memcpy(read_buf, read, read_len);

    static int dp[MAX_LEN + 1][MAX_LEN + 1];
    static int bt[MAX_LEN + 1][MAX_LEN + 1];

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
            int diag = dp[i-1][j-1] + (ref_buf[i-1] == read_buf[j-1] ? MATCH : MISMATCH);
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
#ifdef __DEBUG__
    char aln_ref[2*MAX_LEN];
    char aln_read[2*MAX_LEN];
    char aln_mid[2*MAX_LEN];
    int aln_len = 0;
#endif
    while (i > 0 && j > 0) {
        int op = bt[i][j];

        if (op == 0) {
#ifdef __DEBUG__
            aln_ref[aln_len]  = ref_buf[i-1];
            aln_read[aln_len] = read_buf[j-1];
            aln_mid[aln_len]  = (ref_buf[i-1] == read_buf[j-1]) ? '|' : '*';
#endif
            if (ref_buf[i-1] != read_buf[j-1]) {
                mismatches++;
                if (i-1 != 0 && j-1 != 0 && i != ref_len && j != read_len) {
                    // printf("SNP\tREFID=%lu\tREADID=%lu\tPOS=%lu\tREF=%c\tALT=%c\n", ref_id, read_id, ref_pos + i - 1, ref_buf[i-1], read_buf[j-1]);
                }
            }
            i--; j--;
        }
        else if (op == 1) {
#ifdef __DEBUG__
            aln_ref[aln_len]  = ref_buf[i-1];
            aln_read[aln_len] = '-';
            aln_mid[aln_len]  = ' ';
            mismatches++;
#endif
            if (i-1 != 0 && j-1 != 0 && i != ref_len && j != read_len) {
                // printf("DEL\tREFID=%lu\tREADID=%lu\tPOS=%lu\tREF=%c\tALT=%c\n", ref_id, read_id, ref_pos + i - 1, ref_buf[i-1], '-');
            }
            i--;
        }
        else {
#ifdef __DEBUG__
            aln_ref[aln_len]  = '-';
            aln_read[aln_len] = read_buf[j-1];
            aln_mid[aln_len]  = ' ';
            mismatches++;
#endif
            if (i-1 != 0 && j-1 != 0 && i != ref_len && j != read_len) {
                // printf("INS\tREFID=%lu\tREADID=%lu\tPOS=%lu\tREF=%c\tALT=%c\n", ref_id, read_id, ref_pos + i - 1, '-', read_buf[j-1]);
            }
            j--;
        }

        // if (mismatches > 20) return 0;
#ifdef __DEBUG__
        aln_len++;
#endif
    }
#ifdef __DEBUG__
    if (mismatches) {
        // reverse alignment strings
        for (int x = 0; x < aln_len / 2; x++) {
            char t;

            t = aln_ref[x];
            aln_ref[x] = aln_ref[aln_len - 1 - x];
            aln_ref[aln_len - 1 - x] = t;

            t = aln_mid[x];
            aln_mid[x] = aln_mid[aln_len - 1 - x];
            aln_mid[aln_len - 1 - x] = t;

            t = aln_read[x];
            aln_read[x] = aln_read[aln_len - 1 - x];
            aln_read[aln_len - 1 - x] = t;
        }

        aln_ref[aln_len]  = '\0';
        aln_mid[aln_len]  = '\0';
        aln_read[aln_len] = '\0';

        // print alignment
        printf("\n");
        printf("REF : %s\n", aln_ref);
        printf("      %s\n", aln_mid);
        printf("READ: %s\n", aln_read);

        printf("--------------------------------------------------------------------------\n");

        // print seqs
        printf("REF:  %.*s (%d)\n", ref_len, ref, ref_len);
        printf("READ: %.*s (%d)\n\n", read_len, read, read_len);

        for (int i = 0; i < read_len; i++) {
            printf("'%c' ", read[i]);
        }
        printf(" read: %lu, Pos: %lu, Len: %d\n", read_id, read_pos, len);
    }
#endif
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