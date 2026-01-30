#include "faves.h"


KSEQ_INIT(gzFile, gzread)


static inline int cmp_fuzzy_seeds(const void *a, const void *b) {
    const uint128_t *x = (uint128_t *)a;
    const uint128_t *y = (uint128_t *)b;
    return (x->x >> 14) > (y->x >> 14) ? 1 : -1;
}

uint64_t process_fasta(params *p, uint128_t **fuzzy_seeds, uint64_t *fuzzy_seeds_len, map32_t **index_table, ref_seq **seqs, int *chrom_count) {

    *chrom_count = store_seqs(p->fasta, seqs);
    *fuzzy_seeds_len = 0;

    uint128_t *all_fuzzy_seeds;
    uint64_t all_fuzzy_seeds_len = 0;
    uint64_t all_fuzzy_seeds_cap = SKETCH_CAPACITY;

    all_fuzzy_seeds = (uint128_t *)malloc(sizeof(uint128_t) * all_fuzzy_seeds_cap);
    if (!all_fuzzy_seeds) {
        fprintf(stderr, "[ERROR] couldn't allocate array\n");
        exit(1);
    }

    gzFile fp = gzopen(p->fasta, "r");
    if (!fp) {
        fprintf(stderr, "[ERROR] cannot open reference %s\n", p->fasta);
        exit(1);
    }

    kseq_t *seq = kseq_init(fp);

    int chrom_index = 0;

    while (kseq_read(seq) >= 0) {
        const char *bases = seq->seq.s;
        int len = seq->seq.l;

        (*seqs[chrom_index]).chrom = (char *)malloc(len);
        memcpy((*seqs[chrom_index]).chrom, bases, len);

        uint128_t *chr_fuzzy_seeds;
        uint64_t chr_fuzzy_seeds_len = 0;
        chr_fuzzy_seeds_len = blend_sketch(bases, len, p->w, p->k, p->blend_bits, p->n_neighbors, chrom_index, &chr_fuzzy_seeds);

        if (chr_fuzzy_seeds_len) {
            if (all_fuzzy_seeds_cap <= all_fuzzy_seeds_len + chr_fuzzy_seeds_len) {
                while (all_fuzzy_seeds_cap <= all_fuzzy_seeds_len + chr_fuzzy_seeds_len) {
                    all_fuzzy_seeds_cap *= 2;
                }
                uint128_t *temp = (uint128_t *)realloc(all_fuzzy_seeds, sizeof(uint128_t) * all_fuzzy_seeds_cap);
                if (!temp) {
                    fprintf(stderr, "[ERROR] couldn't reallocate array\n");
                    exit(1);
                }
                all_fuzzy_seeds = temp;
            }
            
            memcpy(all_fuzzy_seeds + all_fuzzy_seeds_len, chr_fuzzy_seeds, sizeof(uint128_t) * chr_fuzzy_seeds_len);
            all_fuzzy_seeds_len += chr_fuzzy_seeds_len;

            free(chr_fuzzy_seeds);
        }

        chrom_index++;
    }
    // clean-up fasta reading 
    kseq_destroy(seq);
    gzclose(fp);

    // if no seeds are found, then no need to proceed
    if (!all_fuzzy_seeds_len) {
        free(all_fuzzy_seeds);
        return 0;
    }

    uint128_t *temp_all_fuzzy_seeds = (uint128_t *)malloc(sizeof(uint128_t) * all_fuzzy_seeds_len);
    if (!temp_all_fuzzy_seeds) {
        fprintf(stderr, "[ERROR] couldn't allocate array\n");
    }
    memcpy(temp_all_fuzzy_seeds, all_fuzzy_seeds, sizeof(uint128_t) * all_fuzzy_seeds_len);

    // get unique fuzzy_seeds
    qsort(all_fuzzy_seeds, all_fuzzy_seeds_len, sizeof(uint128_t), cmp_fuzzy_seeds);
    
    uint64_t unique_fuzzy_seeds_len = 0;

    for (uint64_t i = 0; i < all_fuzzy_seeds_len; ) {
        uint64_t j = i + 1;

        while (j < all_fuzzy_seeds_len && __blend_get_kmer(all_fuzzy_seeds[j]) == __blend_get_kmer(all_fuzzy_seeds[i])) {
            j++;
        }

        if (j == i + 1) {
            all_fuzzy_seeds[unique_fuzzy_seeds_len++] = all_fuzzy_seeds[i];
        }

        i = j;
    }
    
    map32_t *temp_index_table = *index_table;
    for (uint64_t i = 0; i < unique_fuzzy_seeds_len; i++) {
        khint_t k; int absent;
        k = map32_put(temp_index_table, __blend_get_kmer(all_fuzzy_seeds[i]), &absent);
        kh_val(temp_index_table, k) = i;
    }

    if (all_fuzzy_seeds) {
        *fuzzy_seeds = all_fuzzy_seeds;
    }

    uint64_t relaxed_fuzzy_seeds_len = unique_fuzzy_seeds_len;

    // todo: handle i=0 and i=all_fuzzy_seeds_len-1
    for (uint64_t i = 1; i < all_fuzzy_seeds_len - 1; i++) {
        khint_t k = map32_get(temp_index_table, __blend_get_kmer(temp_all_fuzzy_seeds[i]));
        
        if (k < kh_end(temp_index_table)) { // if unique
            continue;
        }

        k = map32_get(temp_index_table, __blend_get_kmer(temp_all_fuzzy_seeds[i-1]));

        if (k < kh_end(temp_index_table)) { // if left is unique
            all_fuzzy_seeds[relaxed_fuzzy_seeds_len++] = temp_all_fuzzy_seeds[i];
            continue;
        }

        k = map32_get(temp_index_table, __blend_get_kmer(temp_all_fuzzy_seeds[i+1]));

        if (k < kh_end(temp_index_table)) { // if right is unique
            all_fuzzy_seeds[relaxed_fuzzy_seeds_len++] = temp_all_fuzzy_seeds[i];
        }
    }

    // no need for temp_all_fuzzy_seeds any more
    free(temp_all_fuzzy_seeds);

    if (relaxed_fuzzy_seeds_len != unique_fuzzy_seeds_len) {
        qsort(all_fuzzy_seeds + unique_fuzzy_seeds_len, relaxed_fuzzy_seeds_len - unique_fuzzy_seeds_len, sizeof(uint128_t), cmp_fuzzy_seeds);
    }

    for (uint64_t i = unique_fuzzy_seeds_len; i < relaxed_fuzzy_seeds_len; ) {
        uint64_t j = i + 1;

        while (j < relaxed_fuzzy_seeds_len && __blend_get_kmer(all_fuzzy_seeds[j]) == __blend_get_kmer(all_fuzzy_seeds[i])) {
            j++;
        }

        khint_t k; int absent;
        k = map32_put(temp_index_table, __blend_get_kmer(all_fuzzy_seeds[i]), &absent);
        kh_val(temp_index_table, k) = i;

        i = j;
    }

    *fuzzy_seeds_len = relaxed_fuzzy_seeds_len;

    printf("[INFO] Processed reference (unique %lu (%.2f), relaxed %lu (%.2f), total %lu)\n", unique_fuzzy_seeds_len, (double)unique_fuzzy_seeds_len / (double)all_fuzzy_seeds_len, relaxed_fuzzy_seeds_len, (double)relaxed_fuzzy_seeds_len / (double)all_fuzzy_seeds_len, all_fuzzy_seeds_len);

    return unique_fuzzy_seeds_len;
}

void process_records(params *p, map32_t *index_table, uint128_t *fuzzy_seeds, uint64_t relaxed_fuzzy_seeds_len, ref_seq *seqs, uint64_t unique_fuzzy_seeds_len) {

    gzFile fp = gzopen(p->fastq, "r");
    if (!fp) {
        fprintf(stderr, "[ERROR] cannot open reads %s\n", p->fastq);
        exit(1);
    }

    kseq_t *seq = kseq_init(fp);

    time_t last_print = 0;

    // init threads
    job_queue_t queue;
    queue_init(&queue, QUEUE_SIZE);

    pthread_t threads[p->n_threads];
    worker_ctx_t ctx[p->n_threads];

    for (int i = 0; i < p->n_threads; i++) {
        ctx[i].queue = &queue;
        ctx[i].p = malloc(sizeof(params_lite) * p->n_neighbors);
        ctx[i].p->k = p->k;
        ctx[i].p->w = p->w;
        ctx[i].p->blend_bits = p->blend_bits;
        ctx[i].p->n_neighbors = p->n_neighbors;
        ctx[i].p->seqs = seqs;
        ctx[i].p->fuzzy_seeds = fuzzy_seeds;
        ctx[i].p->unique_fuzzy_seeds_len = unique_fuzzy_seeds_len;
        ctx[i].p->relaxed_fuzzy_seeds_len = relaxed_fuzzy_seeds_len;
        ctx[i].p->index_table = index_table;
        pthread_create(&threads[i], NULL, worker_process_record, (ctx + i));
    }

    uint64_t read_count = 0;


    batch_records_t *b = batch_records_create(BATCH_SIZE);

    while (kseq_read(seq) >= 0) {

        b->records[b->len].bases = strdup(seq->seq.s);
        b->records[b->len].len  = seq->seq.l;
        b->len++;

        read_count++;

        if (b->len == BATCH_SIZE) {
            queue_push(&queue, b);
            b = batch_records_create(BATCH_SIZE);
        }

        if (p->progress) {
            time_t now = time(NULL);
            if (now - last_print >= 1) {   // every 1 second
                printf("\r[INFO] Processed %ld records...", read_count);
                fflush(stdout);
                last_print = now;
            }
        }
    }

    if (b->len > 0) {
        queue_push(&queue, b);
    } else {
        batch_records_destroy(b);
    }

    kseq_destroy(seq);
    gzclose(fp);

    pthread_mutex_lock(&queue.lock);
    queue.done = 1;
    pthread_cond_broadcast(&queue.not_empty);
    pthread_mutex_unlock(&queue.lock);

    for (int i = 0; i < p->n_threads; i++) {
        pthread_join(threads[i], NULL);
    }
    
    free(queue.buf);
    pthread_mutex_destroy(&queue.lock);
    pthread_cond_destroy(&queue.not_empty);
    pthread_cond_destroy(&queue.not_full);

    stats_t stats = (stats_t){0ULL, 0ULL, 0ULL, 0ULL};
    for (int i = 0; i < p->n_threads; i++) {
        stats.found += ctx[i].p->stats.found;
        stats.not_found += ctx[i].p->stats.not_found;
        stats.mismatches += ctx[i].p->stats.mismatches;
        stats.total_len += ctx[i].p->stats.total_len;
        free(ctx[i].p);
    }

    if (p->progress) {
        printf("\n");
    }

    printf("[INFO] Found: %lu, not found: %lu, mismatches: %lu, total: %lu\n", stats.found, stats.not_found, stats.mismatches, stats.total_len);

    // TODO: sort reported variatns and take consensus and finalize report
}


void *worker_process_record(void *arg) {

    worker_ctx_t *ctx = (worker_ctx_t*)arg;

    int w = ctx->p->w;
    int k = ctx->p->k;
    int blend_bits = ctx->p->blend_bits;
    int n_neighbors = ctx->p->n_neighbors;

    ref_seq *seqs = ctx->p->seqs;
    uint128_t *fuzzy_seeds = (uint128_t *)ctx->p->fuzzy_seeds;
    uint64_t unique_fuzzy_seeds_len = ctx->p->unique_fuzzy_seeds_len;
    uint64_t relaxed_fuzzy_seeds_len = ctx->p->relaxed_fuzzy_seeds_len;
    map32_t *index_table = (map32_t *)ctx->p->index_table;

    ctx->p->stats = (stats_t){0ULL, 0ULL, 0ULL, 0ULL};
    ctx->p->stats.not_found = 0;
    ctx->p->stats.mismatches = 0;
    ctx->p->stats.total_len = 0;

    uint64_t found = 0, not_found = 0, mismatches = 0, total_len = 0;

    while (1) {
        batch_records_t *records = queue_pop(ctx->queue);
        if (!records) break;

        for (int i = 0; i < records->len; i++) {

            char *bases = records->records[i].bases;
            int len = records->records[i].len;

            uint128_t *temp_fuzzy_seeds;
            uint64_t temp_fuzzy_seeds_len = 0;
            temp_fuzzy_seeds_len = blend_sketch(bases, len, w, k, blend_bits, n_neighbors, 0, &temp_fuzzy_seeds);
            
            // store in set
            for (uint64_t i = 0; i < temp_fuzzy_seeds_len; i++) {
                khint_t k = map32_get(index_table, __blend_get_kmer(temp_fuzzy_seeds[i]));

                if (k < kh_end(index_table)) {

                    found++;

                    uint32_t seed_index = kh_val(index_table, k);

                    if (seed_index < unique_fuzzy_seeds_len) { // check if unique
                        // uint64_t ref_kmer = __blend_get_kmer(fuzzy_seeds[seed_index]);
                        uint64_t ref_span = __blend_get_length(fuzzy_seeds[seed_index]);
                        uint64_t ref_index = __blend_get_index(fuzzy_seeds[seed_index]);
                        uint64_t ref_id = __blend_get_reference_id(fuzzy_seeds[seed_index]);
                        int ref_strand = __blend_get_strand(fuzzy_seeds[seed_index]);

                        uint64_t read_span = __blend_get_length(temp_fuzzy_seeds[i]);
                        uint64_t read_index = __blend_get_index(temp_fuzzy_seeds[i]);
                        int read_strand = __blend_get_strand(temp_fuzzy_seeds[i]);

                        // get substrings
                        const char *ref = seqs[ref_id].chrom + ref_index;
                        const char *read = bases + read_index;

                        // align and report variants in alignment
                        int alignment_mismatches = banded_align_and_report(ref, ref_span, read, read_span, ref_index, ref_id, ref_strand, read_strand, read_index, len);
                        mismatches += alignment_mismatches;

                        continue;
                    }

                    if (i) { // check left
                        k = map32_get(index_table, __blend_get_kmer(temp_fuzzy_seeds[i-1]));

                        if (k < kh_end(index_table)) {

                            uint32_t left_seed_index = kh_val(index_table, k);

                            if (left_seed_index < unique_fuzzy_seeds_len) {

                                for (uint64_t temp_seed_index = seed_index; seed_index < relaxed_fuzzy_seeds_len; temp_seed_index++) {
                                    if (__blend_get_kmer(fuzzy_seeds[temp_seed_index]) != __blend_get_kmer(temp_fuzzy_seeds[i])) {
                                        break;
                                    }
                                    if (__blend_get_reference_id(fuzzy_seeds[temp_seed_index]) == __blend_get_reference_id(fuzzy_seeds[left_seed_index]) &&
                                        abs_diff(__blend_get_index(fuzzy_seeds[temp_seed_index]), __blend_get_index(fuzzy_seeds[left_seed_index])) < DISTANCE_THRESHOLD) {
                                        // found
                                        uint64_t ref_span = __blend_get_length(fuzzy_seeds[temp_seed_index]);
                                        uint64_t ref_index = __blend_get_index(fuzzy_seeds[temp_seed_index]);
                                        uint64_t ref_id = __blend_get_reference_id(fuzzy_seeds[temp_seed_index]);
                                        int ref_strand = __blend_get_strand(fuzzy_seeds[temp_seed_index]);

                                        uint64_t read_span = __blend_get_length(temp_fuzzy_seeds[i]);
                                        uint64_t read_index = __blend_get_index(temp_fuzzy_seeds[i]);
                                        int read_strand = __blend_get_strand(temp_fuzzy_seeds[i]);

                                        // get substrings
                                        const char *ref = seqs[ref_id].chrom + ref_index;
                                        const char *read = bases + read_index;

                                        // align and report variants in alignment
                                        int alignment_mismatches = banded_align_and_report(ref, ref_span, read, read_span, ref_index, ref_id, ref_strand, read_strand, read_index, len);
                                        mismatches += alignment_mismatches;
                                    }
                                }
                            }
                        }
                    } else if (i < temp_fuzzy_seeds_len -1) { // check right
                        k = map32_get(index_table, __blend_get_kmer(temp_fuzzy_seeds[i+1]));

                        if (k < kh_end(index_table)) {

                            uint32_t right_seed_index = kh_val(index_table, k);
                            
                            if (right_seed_index < unique_fuzzy_seeds_len) {

                                for (uint64_t temp_seed_index = seed_index; seed_index < relaxed_fuzzy_seeds_len; temp_seed_index++) {
                                    if (__blend_get_reference_id(fuzzy_seeds[temp_seed_index]) != __blend_get_kmer(temp_fuzzy_seeds[i])) {
                                        break;
                                    }
                                    if (__blend_get_reference_id(fuzzy_seeds[temp_seed_index]) == __blend_get_reference_id(fuzzy_seeds[right_seed_index]) &&
                                        abs_diff(__blend_get_index(fuzzy_seeds[temp_seed_index]), __blend_get_index(fuzzy_seeds[right_seed_index])) < DISTANCE_THRESHOLD) {
                                        // found
                                        uint64_t ref_span = __blend_get_length(fuzzy_seeds[temp_seed_index]);
                                        uint64_t ref_index = __blend_get_index(fuzzy_seeds[temp_seed_index]);
                                        uint64_t ref_id = __blend_get_reference_id(fuzzy_seeds[temp_seed_index]);
                                        int ref_strand = __blend_get_strand(fuzzy_seeds[temp_seed_index]);

                                        uint64_t read_span = __blend_get_length(temp_fuzzy_seeds[i]);
                                        uint64_t read_index = __blend_get_index(temp_fuzzy_seeds[i]);
                                        int read_strand = __blend_get_strand(temp_fuzzy_seeds[i]);

                                        // get substrings
                                        const char *ref = seqs[ref_id].chrom + ref_index;
                                        const char *read = bases + read_index;

                                        // align and report variants in alignment
                                        int alignment_mismatches = banded_align_and_report(ref, ref_span, read, read_span, ref_index, ref_id, ref_strand, read_strand, read_index, len);
                                        mismatches += alignment_mismatches;
                                    }
                                }
                            }
                        }
                    }
                } else {
                    not_found++;
                }
            }

            if (temp_fuzzy_seeds_len) free(temp_fuzzy_seeds);

            total_len += temp_fuzzy_seeds_len;
        }

        // cleanup
        batch_records_destroy(records);
    }

    ctx->p->stats.found = found;
    ctx->p->stats.not_found = not_found;
    ctx->p->stats.mismatches = mismatches;
    ctx->p->stats.total_len = total_len;
    return NULL;
}