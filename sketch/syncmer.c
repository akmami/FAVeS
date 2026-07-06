#include "syncmer.h"

static inline uint64_t murmur64(uint64_t h) {
    h ^= (h >> 33);
    h *= 0xff51afd7ed558ccdULL;
    h ^= (h >> 33);
    h *= 0xc4ceb9fe1a85ec53ULL;
    h ^= (h >> 33);
    return h;
}

static inline void push_output(uint128_t **out, uint64_t *n, uint64_t *cap, uint128_t v) {
    if (*n == *cap) {
        *cap = (*cap < 16) ? 16 : (*cap << 1);
        uint128_t *p = (uint128_t *)realloc(*out, (*cap) * sizeof(uint128_t));
        if (p == NULL) {
            exit(1);
        }
        *out = p;
    }
    (*out)[(*n)++] = v;
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

uint64_t sketch_syncmers(const char *str, int len, int smer_size, int kmer_size, uint32_t reference_id, uint128_t **syncmers) {
    *syncmers = NULL;

    if (str == NULL || len <= 0) {
        return 0;
    }

    if (smer_size <= 0 || kmer_size <= 0 || smer_size > kmer_size) {
        return 0;
    }

    if (smer_size > 32 || kmer_size > 32) {
        return 0;
    }

    if (len < kmer_size) {
        return 0;
    }

    const size_t s = (size_t)smer_size;
    const size_t k = (size_t)kmer_size;
    const size_t smers_per_kmer = k - s + 1;

    uint64_t *smer_hash = (uint64_t *)malloc(smers_per_kmer * sizeof(uint64_t));
    size_t *smer_pos  = (size_t *)malloc(smers_per_kmer * sizeof(size_t));

    if (smer_hash == NULL || smer_pos == NULL) {
        free(smer_hash);
        free(smer_pos);
        exit(1);
    }

    uint128_t *kmers = NULL;
    uint64_t kmers_len = 0;
    uint64_t kmers_cap = 0;

    uint64_t fwd_kmer = 0;
    uint64_t fwd_smer = 0;

    const uint64_t kmer_mask = (k == 32) ? UINT64_MAX : ((1ULL << (2ULL * k)) - 1ULL);
    const uint64_t smer_mask = (s == 32) ? UINT64_MAX : ((1ULL << (2ULL * s)) - 1ULL);

    size_t valid_run = 0;

    for (size_t char_index = 0; char_index < (size_t)len; char_index++) {
        uint64_t current_char = seq_nt4_table[(uint8_t)str[char_index]];

        if (current_char == 4) {
            valid_run = 0;
            fwd_kmer = 0;
            fwd_smer = 0;
            continue;
        }

        fwd_kmer = ((fwd_kmer << 2) | current_char) & kmer_mask;
        fwd_smer = ((fwd_smer << 2) | current_char) & smer_mask;

        valid_run++;

        if (valid_run >= s) {
            size_t smer_start = char_index - s + 1;
            size_t idx = smer_start % smers_per_kmer;

            smer_hash[idx] = murmur64(fwd_smer);
            smer_pos[idx] = smer_start;
        }

        if (valid_run >= k) {
            size_t kmer_start = char_index - k + 1;
            size_t first_smer_pos = kmer_start;
            size_t last_smer_pos = char_index - s + 1;

            uint64_t min_hash = UINT64_MAX;
            size_t min_pos = SIZE_MAX;

            for (size_t j = 0; j < smers_per_kmer; j++) {
                size_t pos = smer_pos[j];

                if (pos < first_smer_pos || pos > last_smer_pos) {
                    continue;
                }

                if (smer_hash[j] < min_hash ||
                    (smer_hash[j] == min_hash && pos < min_pos)) {
                    min_hash = smer_hash[j];
                    min_pos = pos;
                }
            }

            size_t min_offset = min_pos - kmer_start;

            if (min_offset == 0 || min_offset == k - s) {
                uint64_t strand = 0;

                uint128_t info;
                info.x = murmur64(fwd_kmer);
                info.y = ((uint64_t)reference_id << 32) | (((uint64_t)(uint32_t)kmer_start) << 1) | strand;
                push_output(&kmers, &kmers_len, &kmers_cap, info);
            }
        }
    }

    free(smer_hash);
    free(smer_pos);

    *syncmers = kmers;
    return kmers_len;
}