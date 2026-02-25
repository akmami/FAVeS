#include "faves.h"


KSEQ_INIT(gzFile, gzread)

static char int2var[5] = {'A', 'C', 'G', 'T', 'N'};

BVEC_INIT(uint64, uint64_t)

#define align(seqs, read, ref_seed, read_seed, mismatch_count) { \
    /* uint64_t ref_kmer = __sketch_get_kmer(ref_seed); */ \
    uint64_t ref_span = __sketch_get_length(ref_seed); \
    uint64_t ref_index = __sketch_get_index(ref_seed); \
    uint64_t ref_id = __sketch_get_reference_id(ref_seed); \
    /* int ref_strand = __sketch_get_strand(read_seed); */ \
    uint64_t read_span = __sketch_get_length(read_seed); \
    uint64_t read_index = __sketch_get_index(read_seed); \
    int read_strand = __sketch_get_strand(read_seed); \
    /* get substrings */ \
    const char *ref = seqs[ref_id].chrom + ref_index; \
    const char *read = bases + read_index; \
    /* align and report variants in alignment */ \
    int alignment_mismatches = banded_align_and_report(ref, ref_span, read, read_span, read_strand, ref_index, ref_id, variants); \
    mismatch_count += alignment_mismatches; \
}

#define pack_id(x, y) ((((x) & 0xFFFFFFFF) | ((y) << 32)))
#define unpack_ref(x) ((uint32_t)((x) & 0xFFFFFFFF))
#define unpack_read(x) ((uint32_t)((x) >> 32))

static inline int cmp_variations(const void *a, const void *b) {
    uint64_t x = *(const uint64_t *)a;
    uint64_t y = *(const uint64_t *)b;
    return (x > y) - (x < y);
}

static inline int cmp_seeds(const void *a, const void *b) {
    const uint64_t x = __sketch_get_kmer(*(uint128_t *)a);
    const uint64_t y = __sketch_get_kmer(*(uint128_t *)b);
    return (x > y) - (x < y);
}

void process_fasta(params_t *p, uint128_t **seeds, uint64_t *seeds_len, map32_t **index_table, ref_seq_t **seqs, int *chrom_count) {

    *chrom_count = store_seqs(p->fasta, seqs);
    *seeds_len = 0;

    uint128_t *all_seeds;
    uint64_t all_seeds_len = 0;
    uint64_t all_seeds_cap = __DEFAULT_SKETCH_CAPACITY__;

    all_seeds = (uint128_t *)malloc(sizeof(uint128_t) * all_seeds_cap);
    if (!all_seeds) {
        fprintf(stderr, "[%s::err] couldn't allocate array\n", __TOOL_SHORT_NAME__);
        exit(1);
    }

    gzFile fp = gzopen(p->fasta, "r");
    if (!fp) {
        fprintf(stderr, "[%s::err] cannot open reference file %s\n", __TOOL_SHORT_NAME__, p->fasta);
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

        uint128_t *chr_seeds;
        uint64_t chr_seeds_len = 0;
        chr_seeds_len = blend_sb_sketch(bases, len, p->w, p->k, p->blend_bits, p->n_neighbors, chrom_index, &chr_seeds);
        
        if (chr_seeds_len) {
            if (all_seeds_cap <= all_seeds_len + chr_seeds_len) {
                while (all_seeds_cap <= all_seeds_len + chr_seeds_len) {
                    all_seeds_cap *= 2;
                }
                uint128_t *temp = (uint128_t *)realloc(all_seeds, sizeof(uint128_t) * all_seeds_cap);
                if (!temp) {
                    fprintf(stderr, "[%s::err] couldn't reallocate array\n", __TOOL_SHORT_NAME__);
                    exit(1);
                }
                all_seeds = temp;
            }
            
            memcpy(all_seeds + all_seeds_len, chr_seeds, sizeof(uint128_t) * chr_seeds_len);
            all_seeds_len += chr_seeds_len;

            free(chr_seeds);
        }

        chrom_index++;
    }
    // clean-up fasta reading 
    kseq_destroy(seq);
    gzclose(fp);

    // if no seeds are found, then no need to proceed
    if (!all_seeds_len) {
        free(all_seeds);
        return;
    }

    // find unique seeds
    map8_t *dup_map = map8_init();
    if (!dup_map) {
        fprintf(stderr, "[%s::err] couldn't allocate maps\n", __TOOL_SHORT_NAME__);
        exit(1);
    }

    for (uint64_t i = 0; i < all_seeds_len; i++) {
        uint32_t kmer = (uint32_t)__sketch_get_kmer(all_seeds[i]);

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
    uint64_t unique_seeds_len = 0;
    double total_unique_seed_span = 0;

    for (uint64_t i = 0; i < all_seeds_len; i++) {
        uint32_t kmer = (uint32_t)__sketch_get_kmer(all_seeds[i]);

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

        unique_seeds_len++;
        total_unique_seed_span += __sketch_get_length(all_seeds[i]);
    }

    // done with temp map
    map8_destroy(dup_map);

    *seeds = all_seeds;
    *seeds_len = all_seeds_len;

    printf("[%s:ref] Processed (uniq #: %lu - %.2f%c, seed #: %lu, span avg: %.2f)\n", __TOOL_SHORT_NAME__, unique_seeds_len, (double)unique_seeds_len / (double)all_seeds_len, '%', all_seeds_len, total_unique_seed_span / unique_seeds_len);
}

void process_records(params_t *p, map32_t *index_table, uint128_t *seeds, uint64_t seeds_len, ref_seq_t *seqs) {

    gzFile fp = gzopen(p->fastq, "r");
    if (!fp) {
        fprintf(stderr, "[%s::err] cannot open reads file %s\n", __TOOL_SHORT_NAME__, p->fastq);
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
        ctx[i].p->seeds = seeds;
        ctx[i].p->seeds_len = seeds_len;
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
                printf("\r[%s:live] Processed %ld records ...", __TOOL_SHORT_NAME__, read_count);
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

    printf("[%s:time] elapsed: %02d:%02d:%02d (%.2f%c idle)\n", __TOOL_SHORT_NAME__, (int)((program_end - program_start) / 3600), (int)((program_end - program_start) / 60) % 3600, ((int)(program_end - program_start) % 60), total_exec_sec / (program_end - program_start), '%');

    // cleanup unnesesary indices
    map32_destroy(index_table);
    if (seeds_len) free(seeds);

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

    stats_t stats = (stats_t){0ULL, 0ULL, 0ULL, 0ULL, 0ULL, 0ULL, 0ULL};
    for (int i = 0; i < p->n_threads; i++) {
        stats.countains_unique += ctx[i].p->stats.countains_unique;
        stats.only_non_unique += ctx[i].p->stats.only_non_unique;
        stats.no_seed += ctx[i].p->stats.no_seed;
        stats.total_seed_count += ctx[i].p->stats.total_seed_count;
        stats.total_uniq_count += ctx[i].p->stats.total_uniq_count;
        stats.mismatch_count += ctx[i].p->stats.mismatch_count;
        stats.neighbour_count += ctx[i].p->stats.neighbour_count;
        var_free(ctx[i].p->fv_variants);
        free(ctx[i].p);
    }

    printf("[%s:records] query: uniq %lu, non-uniq %lu, no-seed %lu | seed #: %lu, uniq #: %lu, mismatch #: %lu | radious #: %lu\n", __TOOL_SHORT_NAME__, stats.countains_unique, stats.only_non_unique, stats.no_seed, stats.total_seed_count, stats.total_uniq_count, stats.mismatch_count, stats.neighbour_count);

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

    printf("[%s:records] read: var %lu, snp: %lu\n", __TOOL_SHORT_NAME__, total_var_size, total_consensus_var);
}

void *worker_process_record(void *arg) {

    worker_ctx_t *ctx = (worker_ctx_t*)arg;

    int w = ctx->p->w;
    int k = ctx->p->k;
    int blend_bits = ctx->p->blend_bits;
    int n_neighbors = ctx->p->n_neighbors;

    ref_seq_t *seqs = ctx->p->seqs;
    uint128_t *seeds = (uint128_t *)ctx->p->seeds;
    uint64_t seeds_len = ctx->p->seeds_len;
    map32_t *index_table = (map32_t *)ctx->p->index_table;

    uint64_bvec_t *ids = uint64_init(64);

    uint64_t countains_unique = 0, only_non_unique = 0, no_seed = 0, total_seed_count = 0, total_uniq_count = 0, mismatch_count = 0, neighbour_count = 0;

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

            uint128_t *temp_seeds;
            uint64_t temp_seeds_len = 0;
            temp_seeds_len = blend_sb_sketch(bases, len, w, k, blend_bits, n_neighbors, 0, &temp_seeds);
            
            int last_read_anchor_idx = -1;
            int last_ref_anchor_idx = -1;

            uint64_clear(ids);
            
            for (uint64_t read_anchor_idx = 0; read_anchor_idx < temp_seeds_len; read_anchor_idx++) {
                khint_t k = map32_get(index_table, __sketch_get_kmer(temp_seeds[read_anchor_idx]));

                if (k == kh_end(index_table)) continue;

                uint64_t ref_anchor_idx = kh_val(index_table, k);

                if (seeds_len <= ref_anchor_idx) {
                    fprintf(stderr, "[%s:err] Weird stuff is ongoing.\n", __TOOL_SHORT_NAME__);
                    abort();
                }

                uint64_add(ids, pack_id(ref_anchor_idx, read_anchor_idx));
            }
             
            assert(ids->size <= temp_seeds_len);

            total_uniq_count += ids->size;

            for (uint32_t unique_ids_index = 0; unique_ids_index < ids->size; unique_ids_index++) {

                uint32_t ref_anchor_idx = unpack_ref(ids->array[unique_ids_index]);
                uint32_t read_anchor_idx = unpack_read(ids->array[unique_ids_index]);

                uint32_t i = 0;

                for (int seed_iter_index = read_anchor_idx; i <= ref_anchor_idx && last_read_anchor_idx < seed_iter_index; seed_iter_index--) {
                    align(seqs, read, seeds[ref_anchor_idx - i], temp_seeds[seed_iter_index], mismatch_count)
                    neighbour_count++;
                    i++;
                }
                neighbour_count--;

                last_ref_anchor_idx = ref_anchor_idx;
                last_read_anchor_idx = read_anchor_idx;
            }

            if (last_read_anchor_idx != -1) {
                last_ref_anchor_idx++;
                last_read_anchor_idx++;

                while ((uint64_t)last_read_anchor_idx < temp_seeds_len) {
                    align(seqs, read, seeds[last_ref_anchor_idx], temp_seeds[last_read_anchor_idx], mismatch_count)                    
                    last_read_anchor_idx++;
                    last_ref_anchor_idx++;
                }
            }

            if (temp_seeds_len) free(temp_seeds);

            total_seed_count += temp_seeds_len;

            if (ids->size) {
                countains_unique++;
            } else if (temp_seeds_len) {
                only_non_unique++;
            } else {
                no_seed++;
            }
        }

        // cleanup
        batch_records_destroy(records);
    }

    uint64_free(ids);

    ctx->p->stats.countains_unique = countains_unique;
    ctx->p->stats.only_non_unique = only_non_unique;
    ctx->p->stats.no_seed = no_seed;
    ctx->p->stats.total_seed_count = total_seed_count;
    ctx->p->stats.total_uniq_count = total_uniq_count;
    ctx->p->stats.mismatch_count = mismatch_count;
    ctx->p->stats.neighbour_count = neighbour_count;
    return NULL;
}