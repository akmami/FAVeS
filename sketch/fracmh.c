#include "fracmh.h"

// same nt table as sketch/kmer.c (a/A c/C g/G t/T -> 0..3, else 4 == break)
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

// SplitMix64 finalizer -- content-only hash (depends solely on the k-mer bits)
static inline uint64_t hash64(uint64_t x) {
	x ^= x >> 30; x *= 0xbf58476d1ce4e5b9ULL;
	x ^= x >> 27; x *= 0x94d049bb133111ebULL;
	x ^= x >> 31;
	return x;
}

uint64_t sketch_fracmh(const char *str, int len, int kmer_size, uint32_t mod,
                       uint32_t reference_id, uint128_t **kmers) {

	if (len <= 0) return 0;
	if (mod == 0) mod = 1;

	uint64_t kmer_bit_mask = (1ULL << (2 * kmer_size)) - 1;
	uint64_t current_kmer = 0;
	int char_index, kmer_len, kmer_span = 0;

	assert(len > 0 && (kmer_size > 0 && kmer_size <= 28));

	kmer_len = 0;

	uint64_t kmers_cap = (uint64_t)len - kmer_size + 1;
	*kmers = (uint128_t *)malloc(sizeof(uint128_t) * kmers_cap);

	uint128_t *tmp_kmers = *kmers;
	uint64_t kmers_len = 0;

	for (char_index = 0; char_index < len; ++char_index) {
		int current_char = seq_nt4_table[(uint8_t)str[char_index]];

		if (current_char < 4) {
			kmer_span = kmer_len + 1 < kmer_size ? kmer_len + 1 : kmer_size;
			current_kmer = (current_kmer << 2 | current_char) & kmer_bit_mask;
			kmer_len++;
			if (kmer_len >= kmer_size && kmer_span < 256) {
				// FracMinHash gate: keep iff hash(kmer) % mod == 0
				if (hash64(current_kmer) % mod == 0) {
					tmp_kmers[kmers_len].x = current_kmer << 8 | kmer_span;
					tmp_kmers[kmers_len++].y = (uint64_t)reference_id << 32
						| (uint32_t)(char_index - kmer_size + 1) << 1 | 0;
				}
			}
		} else {
			kmer_len = 0, kmer_span = 0;
		}
	}

	return kmers_len;
}
