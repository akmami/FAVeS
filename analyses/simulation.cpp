#include <stdint.h>
#include <math.h>
#include <iostream>
#include "../sketch/blend.h"

#define MAX_SEQ_LEN 10000
#define READ_LEN 10000
#define MAX_SKETCH 10000

#define READ_COUNT 10000
#define READ_COUNT_BREAKPOINT 1000


typedef struct {
    uint32_t max_seq_len;
    uint32_t read_len;
    uint32_t max_sketch;
    int k;
    int w;
    int b;
    int n;
    int progress;
} params_t;

typedef struct {
    uint128_t *anchors;
    uint64_t count;
} sketch_t;

typedef struct {
    uint64_t ref_pos;
    uint64_t read_pos;
} match;

int sample_geometric(double p_extend) {
    int len = 1;
    while (((double)rand() / RAND_MAX) < p_extend)
        len++;
    return len;
}

/* ============================
   SKETCH FUNCTION
   ============================ */

void sketch_blend_sequence(const char *seq, int len, params_t *p, sketch_t *sk) {

    uint128_t *fuzzy_seeds;
    uint64_t fuzzy_seeds_len = 0;

    fuzzy_seeds_len = sketch_blend(seq, len, p->w, p->k, p->b, p->n, 0, &fuzzy_seeds);

    sk->anchors = fuzzy_seeds;
    sk->count = fuzzy_seeds_len;
}

/* ============================
   Evaluation
   ============================ */

int cmp_match(const void *a, const void *b) {
    const match *x = (match *)a;
    const match *y = (match *)b;
    if (x->ref_pos != y->ref_pos) return x->ref_pos - y->ref_pos;
    return x->read_pos - y->read_pos;
}

void evaluate(const sketch_t *ref, const sketch_t *read, int *TP, int *FP, int *FN) {

    /* ----------------------------
       1. Collect candidate matches
       ---------------------------- */

    if (ref->count == 0) {
        *TP = 0;
        *FP = 0;
        *FN = ref->count;
        return;
    }

    uint64_t capacity = 2 * ref->count;
    match *matches = (match *)malloc(sizeof(match) * 2 * ref->count);
    uint64_t m = 0;

    for (uint64_t i = 0; i < ref->count; i++) {
        for (uint64_t j = 0; j < read->count; j++) {
            if (__blend_get_kmer(ref->anchors[i]) == __blend_get_kmer(read->anchors[j])) {
                matches[m].ref_pos = __blend_get_index(ref->anchors[i]);
                matches[m++].read_pos = __blend_get_index(read->anchors[j]);

                if (m == capacity) {
                    capacity *= 2;
                    matches = (match *)realloc(matches, capacity * sizeof(match));

                    if (!matches) {
                        fprintf(stderr, "realloc failed\n");
                        exit(-1);
                    }
                }
            }
        }
    }

    if (m == 0) {
        *TP = 0;
        *FP = 0;
        *FN = ref->count;
        free(matches);
        return;
    }

    /* ----------------------------
       2. Sort by reference position
       ---------------------------- */

    qsort(matches, m, sizeof(match), cmp_match);

    /* ----------------------------
       3. DP chaining
       ---------------------------- */

    int *dp = (int *)malloc(sizeof(int) * m);
    int best = 0;

    for (uint64_t i = 0; i < m; i++) {
        dp[i] = 1;
        for (uint64_t j = 0; j < i; j++) {
            if (matches[j].ref_pos < matches[i].ref_pos && matches[j].read_pos < matches[i].read_pos) {

                if (dp[j] + 1 > dp[i])
                    dp[i] = dp[j] + 1;
            }
        }
        if (dp[i] > best)
            best = dp[i];
    }

    /* ----------------------------
       4. Metrics
       ---------------------------- */

    *TP = best;
    *FP = ref->count - best;
    *FN = read->count - best;

    free(dp);
    free(matches);
}

void print_metrics(int read_len, int tp, int fp, int fn) {
    double precision = tp / (double)(tp + fp + 1e-9);
    double recall = tp / (double)(tp + fn + 1e-9);
    double f1 = 2 * precision * recall / (precision + recall + 1e-9);

    std::cout << "read_len: " << read_len << ',' << std::endl;
    std::cout << "tp" << tp << ',' << std::endl;
    std::cout << "fp" << fp << ',' << std::endl;
    std::cout << "fn" << fn << ',' << std::endl;
    std::cout << "precision" << precision << ',' << std::endl;
    std::cout << "recall" << recall << ',' << std::endl;
    std::cout << "f1" << f1 << ',' << std::endl;
}

/* ============================
   Main experiment
   ============================ */

int main(int argc, char **argv) {

    if (argc < 5) {
        std::cout << "Usage: ./simulation <k> <w> <b> <n>" << std::endl;
        return -1;
    }

    params_t p;
    p.k = atoi(argv[1]);
    p.w = atoi(argv[2]);
    p.b = atoi(argv[3]);
    p.n = atoi(argv[4]);
    p.max_seq_len = MAX_SEQ_LEN;
    p.read_len = READ_LEN;
    p.max_sketch = MAX_SKETCH;

    double error_rates[] = {0.001, 0.005, 0.01, 0.02};
    int n_rates = sizeof(error_rates) / sizeof(double);

    std::cout << '\t' << "{" << std::endl;

    std::cout << '\t' << '\t' << "\"k\": " << p.k << ',' << std::endl;
    std::cout << '\t' << '\t' << "\"w\": " << p.w << ',' << std::endl;
    std::cout << '\t' << '\t' << "\"b\": " << p.b << ',' << std::endl;
    std::cout << '\t' << '\t' << "\"n\": " << p.n << ',' << std::endl;

    std::cout << '\t' << '\t' << "\"stats\": [" << std::endl;

    for (int e = 0; e < n_rates; e++) {

        char maf_name[256];
        snprintf(maf_name, sizeof(maf_name), "../data/sim/sim.%d.maf", (int)(error_rates[e] * 1000));

        FILE *maf = fopen(maf_name, "r");
        if (!maf) {
            continue;
        }

        char *line = (char *)malloc(p.max_seq_len * 2 + 1);
        char *ref = (char *)malloc(p.max_seq_len + 1);
        char *read = (char *)malloc(p.max_seq_len * 2 + 1);

        int line_size = p.max_seq_len * 2 + 1;

        double tp = 0, fp = 0, fn = 0, precision = 0, recall = 0, f1 = 0, read_len = 0;    
        int read_count = 0;

        while (fgets(line, line_size, maf)) {
                    
            /* header */
            if (line[0] != '@') continue;

            /* read ref line */
            if (!fgets(line, line_size, maf)) break;
            line[strcspn(line, "\n")] = '\0';
            strncpy(ref, line, p.max_seq_len);
            ref[p.max_seq_len] = '\0';
            /* '+' separator */
            if (!fgets(line, line_size, maf)) break;
            /* read sequence */
            if (!fgets(line, line_size, maf)) break;
            line[strcspn(line, "\n")] = '\0';
            strncpy(read, line, p.max_seq_len * 2);
            read[p.max_seq_len * 2] = '\0';

            int temp_read_len = strlen(read);

            sketch_t ref_sk, read_sk;
            ref_sk.count = read_sk.count = 0;
            sketch_blend_sequence(ref, p.max_seq_len, &p, &ref_sk);
            sketch_blend_sequence(read, temp_read_len, &p, &read_sk);
                    

            int temp_TP, temp_FP, temp_FN;
            evaluate(&ref_sk, &read_sk, &temp_TP, &temp_FP, &temp_FN);

            double temp_precision = temp_TP / (double)(temp_TP + temp_FP + 1e-9);
            double temp_recall = temp_TP / (double)(temp_TP + temp_FN + 1e-9);
            double temp_f1 = 2 * temp_precision * temp_recall / (temp_precision + temp_recall + 1e-9);
                
            tp += temp_TP;
            fp += temp_FP;
            fn += temp_FN;
            precision += temp_precision;
            recall += temp_recall;
            f1 += temp_f1;
            read_len += temp_read_len;
            read_count++;

            // print_metrics(read_len, TP, FP, FN); // legacy :')

            if (ref_sk.count) free(ref_sk.anchors);
            if (read_sk.count) free(read_sk.anchors);
        }

        tp = tp / read_count;
        fp = fp / read_count;
        fn = fn / read_count;
        precision = precision / read_count;
        recall = recall / read_count;
        f1 = f1 / read_count;
        read_len = read_len / read_count;

        std::cout << '\t' << '\t' << '\t' << "{ "; 
        std::cout << "\"error_rate\": " << error_rates[e] << ", ";
        std::cout << "\"read_count\": " << read_count << ", ";
        std::cout << "\"read_len\": " << read_len << ", ";
        std::cout << "\"tp\": " << tp << ", ";
        std::cout << "\"fp\": " << fp << ", ";
        std::cout << "\"fn\": " << fn << ", ";
        std::cout << "\"precision\": " << precision << ", ";
        std::cout << "\"recall\": " <<  recall << ", ";
        std::cout << "\"f1\": " << f1;
        
        if (e + 1 == n_rates)
            std::cout << " }" << std::endl; 
        else
            std::cout << " }," << std::endl; 
        
        fclose(maf);

        free(line);
        free(ref);
        free(read);
    }

    std::cout << '\t' << '\t' << "]" << std::endl;

    std::cout << '\t' << "}," << std::endl;

    return 0;
}
