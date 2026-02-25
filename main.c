#include <stdlib.h>
#include "struct_def.h"
#include "opt_parser.h"
#include "utils.h"
#include "faves.h"


int main(int argc, char **argv) {

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

    return 0;
}
