#include "faves.h"


KSEQ_INIT(gzFile, gzread)

static char int2var[5] = {'A', 'C', 'G', 'T', 'N'};

static inline int cmp_variations(const void *a, const void *b) {
    uint64_t x = *(const uint64_t *)a;
    uint64_t y = *(const uint64_t *)b;
    return (x > y) - (x < y);
}

static inline int cmp_fuzzy_seeds(const void *a, const void *b) {
    const uint64_t x = (((uint128_t *)a)->x >> 14);
    const uint64_t y = (((uint128_t *)b)->x >> 14);
    return (x > y) - (x < y);
}

void process_fasta(params_t *p, uint128_t **fuzzy_seeds, uint64_t *fuzzy_seeds_len, map32_t **index_table, ref_seq_t **seqs, int *chrom_count) {

    *chrom_count = store_seqs(p->fasta, seqs);
    *fuzzy_seeds_len = 0;

    uint128_t *all_fuzzy_seeds;
    uint64_t all_fuzzy_seeds_len = 0;
    uint64_t all_fuzzy_seeds_cap = __DEFAULT_SKETCH_CAPACITY__;

    all_fuzzy_seeds = (uint128_t *)malloc(sizeof(uint128_t) * all_fuzzy_seeds_cap);
    if (!all_fuzzy_seeds) {
        fprintf(stderr, "[ERROR] couldn't allocate array\n");
        exit(1);
    }

    gzFile fp = gzopen(p->fasta, "r");
    if (!fp) {
        fprintf(stderr, "[ERROR] cannot open reference file %s\n", p->fasta);
        exit(1);
    }

    kseq_t *seq = kseq_init(fp);

    int chrom_index = 0;

    while (kseq_read(seq) >= 0) {
        const char *bases = seq->seq.s;
        int len = seq->seq.l;

        ((*seqs)[chrom_index]).chrom = (char *)malloc(len + 1);
        memcpy(((*seqs)[chrom_index]).chrom, bases, len);
        (*seqs)[chrom_index].chrom[len] = '\0';

        uint128_t *chr_fuzzy_seeds;
        uint64_t chr_fuzzy_seeds_len = 0;
        chr_fuzzy_seeds_len = blend_sb_sketch(bases, len, p->w, p->k, p->blend_bits, p->n_neighbors, chrom_index, &chr_fuzzy_seeds);

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
        return;
    }

    // find unique seeds
    map8_t *dup_map = map8_init();
    if (!dup_map) {
        fprintf(stderr, "[ERROR] couldn't allocate maps\n");
        exit(1);
    }

    for (uint64_t i = 0; i < all_fuzzy_seeds_len; i++) {
        uint32_t kmer = (uint32_t)__blend_get_kmer(all_fuzzy_seeds[i]);

        khint_t k = map8_get(dup_map, kmer);
        if (k == kh_end(dup_map)) {
            // first time seen
            int absent;
            k = map8_put(dup_map, kmer, &absent);
            kh_val(dup_map, k) = 0;
        } else {
            // seen before -> mark duplicate
            khint_t kd = map8_get(dup_map, kmer);
            if (kd != kh_end(dup_map)) {
                kh_val(dup_map, kd) = 1;
            } else {
                // should not happen, but be safe
                int absent;
                kd = map8_put(dup_map, kmer, &absent);
                kh_val(dup_map, kd) = 2;
            }
        }
    }

    map32_t *temp_index_table = *index_table;
    uint64_t unique_fuzzy_seeds_len = 0;
    double total_unique_seed_span = 0;

    for (uint64_t i = 0; i < all_fuzzy_seeds_len; i++) {
        uint32_t kmer = (uint32_t)__blend_get_kmer(all_fuzzy_seeds[i]);

        khint_t kd = map8_get(dup_map, kmer);
        if (kd == kh_end(dup_map)) continue; // safety

        if (kh_val(dup_map, kd) != 0) {
            // duplicate -> skip from unique list
            continue;
        }

        // appears exactly once; add to unique list and map kmer->unique_index
        khint_t k; int absent;
        k = map32_put(temp_index_table, kmer, &absent);
        kh_val(temp_index_table, k) = (uint32_t)i;

        unique_fuzzy_seeds_len++;
        total_unique_seed_span += __blend_get_length(all_fuzzy_seeds[i]);
    }

    // done with temp map
    map8_destroy(dup_map);

    *fuzzy_seeds = all_fuzzy_seeds;
    *fuzzy_seeds_len = all_fuzzy_seeds_len;

    printf("[INFO] Ref processed (uniq #: %lu - %.2f%c, seed #: %lu, span: %.2f)\n", unique_fuzzy_seeds_len, (double)unique_fuzzy_seeds_len / (double)all_fuzzy_seeds_len, '%', all_fuzzy_seeds_len, total_unique_seed_span / unique_fuzzy_seeds_len);
}

void process_records(params_t *p, map32_t *index_table, uint128_t *fuzzy_seeds, uint64_t fuzzy_seeds_len, ref_seq_t *seqs) {

    gzFile fp = gzopen(p->fastq, "r");
    if (!fp) {
        fprintf(stderr, "[ERROR] cannot open reads file %s\n", p->fastq);
        exit(1);
    }

    kseq_t *seq = kseq_init(fp);

    // init threads
    job_queue_t queue;
    queue_init(&queue, __DEFAULT_QUEUE_SIZE__);

    pthread_t threads[p->n_threads];
    worker_ctx_t ctx[p->n_threads];

    for (int i = 0; i < p->n_threads; i++) {
        ctx[i].queue = &queue;
        ctx[i].p = malloc(sizeof(params_lite_t));
        ctx[i].p->k = p->k;
        ctx[i].p->w = p->w;
        ctx[i].p->blend_bits = p->blend_bits;
        ctx[i].p->n_neighbors = p->n_neighbors;
        ctx[i].p->n_threads = p->n_threads;
        ctx[i].p->seqs = seqs;
        ctx[i].p->fuzzy_seeds = fuzzy_seeds;
        ctx[i].p->fuzzy_seeds_len = fuzzy_seeds_len;
        ctx[i].p->index_table = index_table;
        pthread_create(&threads[i], NULL, worker_process_record, (ctx + i));
    }

    uint64_t read_count = 0;

    batch_records_t *b = batch_records_create(__DEFAULT_BATCH_SIZE__);

    time_t program_start = 0, program_end, last_print = 0;
    double total_exec_sec = 0;

    program_start = time(NULL);

    while (kseq_read(seq) >= 0) {

        b->records[b->len].bases = strdup(seq->seq.s);
        b->records[b->len].len  = seq->seq.l;
        b->len++;

        read_count++;

        if (b->len == __DEFAULT_BATCH_SIZE__) {
            time_t queue_start = time(NULL);
            
            queue_push(&queue, b);
            b = batch_records_create(__DEFAULT_BATCH_SIZE__);

            total_exec_sec += time(NULL) - queue_start;
        }

        if (p->progress) {
            time_t now = time(NULL);
            if (now - last_print >= __DEFAULT_PROGRESS_INTERVAL__) {   // every 1 second
                printf("\r[INFO] Processed %ld records ...", read_count);
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

    program_end = time(NULL);

    if (p->progress) {
        printf("\n");
    }

    printf("[INFO] Elapsed time: %02d:%02d:%02d (%.2f%c idle)\n", (int)((program_end - program_start) / 3600), (int)((program_end - program_start) / 60) % 3600, ((int)(program_end - program_start) % 60), total_exec_sec / (program_end - program_start), '%');

    uint64_t total_var_size = 0;
    variation_t *variations = NULL;
    for (int i = 0; i < p->n_threads; i++) {
        total_var_size += ctx[i].p->fv_variants->size;
    }

    if (total_var_size) {
        variations = (variation_t *)malloc(sizeof(variation_t) * total_var_size);
        
        total_var_size = 0;
        for (int i = 0; i < p->n_threads; i++) {
            if (ctx[i].p->fv_variants->size) {
                memcpy(variations + total_var_size, ctx[i].p->fv_variants->array, sizeof(variation_t) * ctx[i].p->fv_variants->size);
            }
            total_var_size += ctx[i].p->fv_variants->size;
        }
    }

    stats_t stats = (stats_t){0ULL, 0ULL, 0ULL, 0ULL, 0ULL, 0ULL};
    for (int i = 0; i < p->n_threads; i++) {
        stats.countains_unique += ctx[i].p->stats.countains_unique;
        stats.only_non_unique += ctx[i].p->stats.only_non_unique;
        stats.no_seed += ctx[i].p->stats.no_seed;
        stats.total_seed_count += ctx[i].p->stats.total_seed_count;
        stats.mismatch_count += ctx[i].p->stats.mismatch_count;
        stats.neighbour_count += ctx[i].p->stats.neighbour_count;
        var_free(ctx[i].p->fv_variants);
        free(ctx[i].p);
    }

    printf("[INFO] query: uniq %lu, non-uniq %lu, noseed %lu | seeds %lu, mm %lu | span #: %lu\n", stats.countains_unique, stats.only_non_unique, stats.no_seed, stats.total_seed_count, stats.mismatch_count, stats.neighbour_count);

    uint64_t total_consensus_var = 0;

    if (total_var_size) {
        // sort reported variants
        qsort(variations, total_var_size, sizeof(variation_t), cmp_variations);

        // take consensus
        size_t i = 0;

        FILE *bed = fopen(p->bed, "w");
        if (!bed) {
            perror("fopen");
            exit(1);
        }

        while (i < total_var_size) {
            variation_t current = variations[i];
            int count = 1;

            // count how many times this variation repeats
            while ((i + count) < total_var_size &&
                __get_variation_chrom(variations[i+count]) == __get_variation_chrom(current) &&
                __get_variation_index(variations[i+count]) == __get_variation_index(current)) {
                count++;
            }

            // if consensus is strong enough, print it
            if (count > p->min_consensus) {
                total_consensus_var++;
                uint64_t chrom = __get_variation_chrom(current);
                uint64_t pos   = __get_variation_index(current);
                char ch_ref = seqs[chrom].chrom[pos] & to_uppercase_mask;
                char ch_alt = int2var[__get_variation_type(current)];

                // BED: 0-based start, 1-based end
                fprintf(bed,
                        "%s\t%lu\t%lu\t%c\t%c\t#=%d\n",
                        seqs[chrom].header,
                        pos,
                        pos + 1,
                        ch_ref,
                        ch_alt,
                        count);
                    }

            // move forward past this run
            i += count;
        }

        fclose(bed);

        // cleanup
        free(variations);
    }

    printf("[INFO] read: var %lu, snp: %lu\n", total_var_size, total_consensus_var);
}

void *worker_process_record(void *arg) {

    worker_ctx_t *ctx = (worker_ctx_t*)arg;

    int w = ctx->p->w;
    int k = ctx->p->k;
    int blend_bits = ctx->p->blend_bits;
    int n_neighbors = ctx->p->n_neighbors;

    ref_seq_t *seqs = ctx->p->seqs;
    uint128_t *fuzzy_seeds = (uint128_t *)ctx->p->fuzzy_seeds;
    uint64_t fuzzy_seeds_len = ctx->p->fuzzy_seeds_len;
    map32_t *index_table = (map32_t *)ctx->p->index_table;

    uint64_t countains_unique = 0, only_non_unique = 0, no_seed = 0, total_seed_count = 0, mismatch_count = 0, neighbour_count = 0;

    ctx->p->fv_variants = var_init(__DEFAULT_VARIANT_CAPACITY__ / ctx->p->n_threads);
    if (!ctx->p->fv_variants) {
        printf("Couldn't allocate variants array\n");
        return NULL;
    }

    var_bvec_t *variants = ctx->p->fv_variants;

    while (1) {
        batch_records_t *records = queue_pop(ctx->queue);
        if (!records) break;

        for (int record_index = 0; record_index < records->len; record_index++) {

            char *bases = records->records[record_index].bases;
            int len = records->records[record_index].len;

            uint128_t *temp_fuzzy_seeds;
            uint64_t temp_fuzzy_seeds_len = 0;
            temp_fuzzy_seeds_len = blend_sb_sketch(bases, len, w, k, blend_bits, n_neighbors, 0, &temp_fuzzy_seeds);
            
            int found_unique = 0;
            int last_valid_seed_index_read = -1;
            int last_valid_seed_index_ref = -1;
            
            for (uint64_t read_seed_index = 0; read_seed_index < temp_fuzzy_seeds_len; read_seed_index++) {
                khint_t k = map32_get(index_table, __blend_get_kmer(temp_fuzzy_seeds[read_seed_index]));

                if (k < kh_end(index_table)) {
                    uint32_t ref_seed_index = kh_val(index_table, k);

                    if (fuzzy_seeds_len <= ref_seed_index) {
                        fprintf(stderr, "[ERROR] Weird stuff is ongoing.\n");
                        abort();
                    }
                    
                    found_unique++;

                    for (int seed_iter_index = read_seed_index; last_valid_seed_index_read < seed_iter_index; seed_iter_index--) {
                        // uint64_t ref_kmer = __blend_get_kmer(fuzzy_seeds[ref_seed_index]);
                        uint64_t ref_span = __blend_get_length(fuzzy_seeds[ref_seed_index]);
                        uint64_t ref_index = __blend_get_index(fuzzy_seeds[ref_seed_index]);
                        uint64_t ref_id = __blend_get_reference_id(fuzzy_seeds[ref_seed_index]);
                        // int ref_strand = __blend_get_strand(fuzzy_seeds[ref_seed_index]);

                        uint64_t read_span = __blend_get_length(temp_fuzzy_seeds[seed_iter_index]);
                        uint64_t read_index = __blend_get_index(temp_fuzzy_seeds[seed_iter_index]);
                        int read_strand = __blend_get_strand(temp_fuzzy_seeds[seed_iter_index]);

                        // get substrings
                        const char *ref = seqs[ref_id].chrom + ref_index;
                        const char *read = bases + read_index;

                        // align and report variants in alignment
                        int alignment_mismatches = banded_align_and_report(ref, ref_span, read, read_span, read_strand, ref_index, ref_id, variants);
                        mismatch_count += alignment_mismatches;

                        neighbour_count++;
                    }
                    neighbour_count--;

                    last_valid_seed_index_ref = ref_seed_index;
                    last_valid_seed_index_read = read_seed_index;
                }
            }

            if (last_valid_seed_index_read != -1) {
                last_valid_seed_index_ref++;
                last_valid_seed_index_read++;

                while ((uint64_t)last_valid_seed_index_read < temp_fuzzy_seeds_len) {
                    // uint64_t ref_kmer = __blend_get_kmer(fuzzy_seeds[seed_index]);
                    uint64_t ref_span = __blend_get_length(fuzzy_seeds[last_valid_seed_index_ref]);
                    uint64_t ref_index = __blend_get_index(fuzzy_seeds[last_valid_seed_index_ref]);
                    uint64_t ref_id = __blend_get_reference_id(fuzzy_seeds[last_valid_seed_index_ref]);
                    // int ref_strand = __blend_get_strand(fuzzy_seeds[seed_index]);

                    uint64_t read_span = __blend_get_length(temp_fuzzy_seeds[last_valid_seed_index_read]);
                    uint64_t read_index = __blend_get_index(temp_fuzzy_seeds[last_valid_seed_index_read]);
                    int read_strand = __blend_get_strand(temp_fuzzy_seeds[last_valid_seed_index_read]);

                    // get substrings
                    const char *ref = seqs[ref_id].chrom + ref_index;
                    const char *read = bases + read_index;

                    // align and report variants in alignment
                    int alignment_mismatches = banded_align_and_report(ref, ref_span, read, read_span, read_strand, ref_index, ref_id, variants);
                    mismatch_count += alignment_mismatches;  
                    
                    last_valid_seed_index_read++;
                    last_valid_seed_index_ref++;
                }
            }

            if (temp_fuzzy_seeds_len) free(temp_fuzzy_seeds);

            total_seed_count += temp_fuzzy_seeds_len;

            if (found_unique) {
                countains_unique++;
            } else if (temp_fuzzy_seeds_len) {
                only_non_unique++;
            } else {
                no_seed++;
            }
        }

        // cleanup
        batch_records_destroy(records);
    }

    ctx->p->stats.countains_unique = countains_unique;
    ctx->p->stats.only_non_unique = only_non_unique;
    ctx->p->stats.no_seed = no_seed;
    ctx->p->stats.mismatch_count = mismatch_count;
    ctx->p->stats.total_seed_count = total_seed_count;
    ctx->p->stats.neighbour_count = neighbour_count;
    return NULL;
}