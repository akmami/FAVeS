#include <stdlib.h>
#include "struct_def.h"
#include "opt_parser.h"
#include "utils.h"
#include "faves.h"


int main(int argc, char **argv) {

    params_t p;
    parse_args(argc, argv, &p);

    fprintf(stderr,
        "Params:\n"
        "   fasta       = %s\n"
        "   fastq       = %s\n"
        "   k           = %d\n"
        "   w           = %d\n"
        "   blend-bits  = %d\n"
        "   n-neighbors = %d\n"
        "   consensus   = %d\n"
        "   threads     = %d\n",
        p.fasta, p.fastq, p.k, p.w,
        p.blend_bits, p.n_neighbors, p.min_consensus, p.n_threads
    );

    map32_t *index_table = map32_init();
    uint128_t *fuzzy_seeds;
    uint64_t fuzzy_seeds_len;
    ref_seq_t *seqs;
    int seq_count = 0;
    process_fasta(&p, &fuzzy_seeds, &fuzzy_seeds_len, &index_table, &seqs, &seq_count);

    process_records(&p, index_table, fuzzy_seeds, fuzzy_seeds_len, seqs);

    // cleanup
    map32_destroy(index_table);
    free_seqs(&seqs, seq_count);

    if (fuzzy_seeds_len) free(fuzzy_seeds);

    return 0;
}
