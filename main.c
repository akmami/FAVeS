#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include "struct_def.h"
#include "opt_parser.h"
#include "utils.h"
#include "faves.h"


int main(int argc, char **argv) {

    struct timespec ts_start, ts_end;
    clock_gettime(CLOCK_MONOTONIC, &ts_start);

    // parse arguments
    params_t p;
    parse_args(argc, argv, &p);

    // process reference and index
    map32_t *index_table = map32_init();
    uint128_t *seeds;
    uint64_t seeds_len;
    ref_seq_t *seqs;
    int seq_count = 0;

    process_fasta(&p, &seeds, &seeds_len, &index_table, &seqs, &seq_count);

    // process reads and call variants. at the end, it will free seqeunces and index
    process_records(&p, index_table, seeds, seeds_len, seqs);

    // cleanup
    free_seqs(&seqs, seq_count);

    clock_gettime(CLOCK_MONOTONIC, &ts_end);
    double total_s = (ts_end.tv_sec - ts_start.tv_sec) + (ts_end.tv_nsec - ts_start.tv_nsec) / 1e9;
    fprintf(stderr, "[%s::time] total: %.3f s\n", __TOOL_SHORT_NAME__, total_s);

    return 0;
}
