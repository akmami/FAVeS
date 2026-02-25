#include "minimizer.h"


#define validate(type, array, size, capacity) { \
	if (size == capacity) {                                               	\
		capacity *= 2;                                                    	\
		type *temp = (type *)realloc(array, capacity * sizeof(type));	\
		if (!temp) { free(array); return 0; }                              	\
		array = temp;                                                    	\
	} 																		\
}

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

static inline uint64_t hash64(uint64_t key, uint64_t mask) {
	key = (~key + (key << 21)) & mask; // key = (key << 21) - key - 1;
	key = key ^ key >> 24;
	key = ((key + (key << 3)) + (key << 8)) & mask; // key * 265
	key = key ^ key >> 14;
	key = ((key + (key << 2)) + (key << 4)) & mask; // key * 21
	key = key ^ key >> 28;
	key = (key + (key << 31)) & mask;
	return key;
}

uint64_t sketch_minimizers(const char *str, int len, int window, int kmer_size, int only_symmetric, uint32_t reference_id, uint128_t **minimizers) {

	if (len <= 0) return 0;

	uint64_t shift = 2 * (kmer_size - 1);
	uint64_t kmer_bit_mask = (1ULL << (2 * kmer_size)) - 1;
	uint64_t current_kmer[2] = {0, 0};
	int char_index, helper_index, buf_pos, min_pos, kmer_len, kmer_span = 0;
	uint128_t buf[256];
	uint128_t min = { UINT64_MAX, UINT64_MAX };

	// legacy code
	assert(len > 0 && (window > 0 && window < 256) && (kmer_size > 0 && kmer_size <= 28)); 
	memset(buf, 0xFF, window * 16);

	kmer_len = 0;
	buf_pos = 0;
	min_pos = 0;

	uint64_t kmers_cap = 3 * len / window;
	*minimizers = (uint128_t *)malloc(sizeof(uint128_t) * kmers_cap);

	uint128_t *kmers = *minimizers;
	uint64_t kmers_len = 0;

	for (char_index = 0; char_index < len; ++char_index) {
		int current_char = seq_nt4_table[(uint8_t)str[char_index]];
		uint128_t current_info = { UINT64_MAX, UINT64_MAX };

		if (current_char < 4) { 
			int strand;
			kmer_span = kmer_len + 1 < kmer_size? kmer_len + 1 : kmer_size;
			current_kmer[0] = (current_kmer[0] << 2 | current_char) & kmer_bit_mask;    // forward k-mer
			current_kmer[1] = (current_kmer[1] >> 2) | (3ULL^current_char) << shift; 	// reverse k-mer
			if (only_symmetric && current_kmer[0] == current_kmer[1]) continue;
			strand = current_kmer[0] < current_kmer[1] ? 0 : 1;
			kmer_len++;
			if (kmer_len >= kmer_size && kmer_span < 256) {
				current_info.x = hash64(current_kmer[strand], kmer_bit_mask) << 8 | kmer_span;
				current_info.y = (uint64_t)reference_id << 32 | (uint32_t)(char_index - kmer_len + 1) << 1 | strand;
			}
		} else {
            if (min.x != UINT64_MAX) {
				validate(uint128_t, kmers, kmers_len, kmers_cap)
				kmers[kmers_len++] = min;
			}

			kmer_len = 0, kmer_span = 0;
            min = (uint128_t){ UINT64_MAX, UINT64_MAX };
			continue;
		}

		buf[buf_pos] = current_info;
		
		if (kmer_len == window + kmer_size - 1 && min.x != UINT64_MAX) { 
			for (helper_index = buf_pos + 1; helper_index < window; ++helper_index) {
				if (min.x == buf[helper_index].x && buf[helper_index].y != min.y) {
					validate(uint128_t, kmers, kmers_len, kmers_cap)
					kmers[kmers_len++] = buf[helper_index];
				}
			}
			for (helper_index = 0; helper_index < buf_pos; ++helper_index) {
				if (min.x == buf[helper_index].x && buf[helper_index].y != min.y) {
					validate(uint128_t, kmers, kmers_len, kmers_cap)
					kmers[kmers_len++] = min;
				}
			}
		}
		
		if (current_info.x <= min.x) {
			if (kmer_len >= window + kmer_size && min.x != UINT64_MAX) {
				validate(uint128_t, kmers, kmers_len, kmers_cap)
				kmers[kmers_len++] = min;
			}
			min = current_info, min_pos = buf_pos;
		} else if (buf_pos == min_pos) {
			if (kmer_len >= window + kmer_size - 1 && min.x != UINT64_MAX) {
				validate(uint128_t, kmers, kmers_len, kmers_cap)
				kmers[kmers_len++] = min;
			}
			for (helper_index = buf_pos + 1, min.x = UINT64_MAX; helper_index < window; ++helper_index) { 
				if (min.x >= buf[helper_index].x) {
					min = buf[helper_index];
					min_pos = helper_index; 
				}
			}
			for (helper_index = 0; helper_index <= buf_pos; ++helper_index) {
				if (min.x >= buf[helper_index].x) {
					min = buf[helper_index];
					min_pos = helper_index;
				}
			}
			if (kmer_len >= window + kmer_size - 1 && min.x != UINT64_MAX) {
				for (helper_index = buf_pos + 1; helper_index < window; ++helper_index) { 
					if (min.x == buf[helper_index].x && min.y != buf[helper_index].y) {
						validate(uint128_t, kmers, kmers_len, kmers_cap)
						kmers[kmers_len++] = buf[helper_index];
					}
				}
				for (helper_index = 0; helper_index <= buf_pos; ++helper_index) {
					if (min.x == buf[helper_index].x && min.y != buf[helper_index].y) {
						validate(uint128_t, kmers, kmers_len, kmers_cap)
						kmers[kmers_len++] = buf[helper_index];
					}
				}
			}
		}

		if (++buf_pos == window) {
			buf_pos = 0;
		}
	}
	if (min.x != UINT64_MAX) {
		validate(uint128_t, kmers, kmers_len, kmers_cap)
		kmers[kmers_len++] = min;
	}

	if (kmers_len == 0) { free(kmers); *minimizers = NULL; }
	else *minimizers = kmers;
	
	return kmers_len;
}