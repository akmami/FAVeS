#include <stdlib.h>
#include <time.h>
#include "struct_def.h"
#include "opt_parser.h"
#include "utils.h"
#include "faves.h"


int main(int argc, char **argv) {

    params p;
    init_params(&p);
    parse_args(argc, argv, &p);

    fprintf(stderr,
        "Params:\n"
        "  fasta       = %s\n"
        "  fastq       = %s\n"
        "  k           = %d\n"
        "  w           = %d\n"
        "  blend-bits  = %d\n"
        "  n-neighbors = %d\n",
        p.fasta, p.fastq, p.k, p.w,
        p.blend_bits, p.n_neighbors
    );

    map32_t *index_table = map32_init();
    uint128_t *fuzzy_seeds;
    uint64_t relaxed_fuzzy_seeds_len, unique_fuzzy_seeds_len;
    ref_seq *seqs;
    int seq_count = 0;
    unique_fuzzy_seeds_len = process_fasta(&p, &fuzzy_seeds, &relaxed_fuzzy_seeds_len, &index_table, &seqs, &seq_count);

    process_records(&p, index_table, fuzzy_seeds, relaxed_fuzzy_seeds_len, seqs, unique_fuzzy_seeds_len);

    // cleanup
    map32_destroy(index_table);
    free_seqs(&seqs, seq_count);

    if (relaxed_fuzzy_seeds_len) free(fuzzy_seeds);

    return 0;
}
