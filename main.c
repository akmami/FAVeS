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
    uint128_t *seeds;
    uint64_t seeds_len;
    ref_seq_t *seqs;
    int seq_count = 0;
    
    // process reference and index
    process_fasta(&p, &seeds, &seeds_len, &index_table, &seqs, &seq_count);

    // process reads and call variants. at the end, it will free seqeunces and index
    process_records(&p, index_table, seeds, seeds_len, seqs);

    // cleanup
    free_seqs(&seqs, seq_count);

    return 0;
}
