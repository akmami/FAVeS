#include <stdlib.h>
#include <time.h>
#include <zlib.h>
#include "utils/khashl.h"
#include "utils/kseq.h"
#include "struct_def.h"
#include "utils.h"
#include "sketch/blend.h"


KHASHL_MAP_INIT(KH_LOCAL, map32_t, map32, uint64_t, uint32_t, kh_hash_uint64, kh_eq_generic)


uint64_t process_fasta(params_t *p, uint128_t **fuzzy_seeds, uint64_t *fuzzy_seeds_len, map32_t **index_table, ref_seq_t **seqs, int *chrom_count);

void process_records(params_t *p, map32_t *index_table, uint128_t *fuzzy_seeds, uint64_t relaxed_fuzzy_seeds_len, ref_seq_t *seqs, uint64_t unique_fuzzy_seeds_len);

void *worker_process_record(void *arg);
