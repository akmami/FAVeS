#include "kmer.h"


static unsigned char seq_nt4_table[256] = {
	0, 1, 2, 3,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  3, 3, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  3, 3, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4
};

uint64_t sketch_kmers(const char *str, int len, int kmer_size, uint32_t reference_id, uint128_t **kmers) {

	if (len <= 0) return 0;

	uint64_t kmer_bit_mask = (1ULL << (2 * kmer_size)) - 1;
	uint64_t current_kmer = 0;
	int char_index, kmer_len, kmer_span = 0;

	// legacy code
	assert(len > 0 && (kmer_size > 0 && kmer_size <= 28)); 

	kmer_len = 0;

	uint64_t kmers_cap = len - kmer_size + 1;
	*kmers = (uint128_t *)malloc(sizeof(uint128_t) * kmers_cap);

	uint128_t *tmp_kmers = *kmers;
	uint64_t kmers_len = 0;

	for (char_index = 0; char_index < len; ++char_index) {
		int current_char = seq_nt4_table[(uint8_t)str[char_index]];

		if (current_char < 4) { 
			int strand;
			kmer_span = kmer_len + 1 < kmer_size ? kmer_len + 1 : kmer_size;
			current_kmer = (current_kmer << 2 | current_char) & kmer_bit_mask;    // forward k-mer
			kmer_len++;	// reverse k-mer
			if (kmer_len >= kmer_size && kmer_span < 256) {
				tmp_kmers[kmers_len].x = current_kmer << 8 | kmer_span;
				tmp_kmers[kmers_len++].y = (uint64_t)reference_id << 32 | (uint32_t)(char_index - kmer_size + 1) << 1 | strand;
			}
		} else {
			kmer_len = 0, kmer_span = 0;
		}

	}
	
	return kmers_len;
}
