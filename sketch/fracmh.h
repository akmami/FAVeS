#ifndef __FRACMH__
#define __FRACMH__

#ifdef __cplusplus
extern "C" {
#endif

#include <stdlib.h>
#include <stdint.h>
#include <assert.h>
#include <string.h>

// same bit layout as sketch/kmer.h
#define __fracmh_get_kmer(seed) 		((seed).x >> 8)
#define __fracmh_get_length(seed) 		((seed).x & 0xFF)
#define __fracmh_get_reference_id(seed) (((seed).y >> 32))
#define __fracmh_get_index(seed) 		(((seed).y & 0xFFFFFFFF) >> 1)
#define __fracmh_get_strand(seed) 		(((seed).y & 1))

#ifndef __UINT128_T__
#define __UINT128_T__
typedef struct {
	uint64_t x; // kmer(56) + span(8)
	uint64_t y; // reference_id(32) + index(31) + strand(1)
} uint128_t;
#endif

/**
 * FracMinHash k-mers: keep every k-mer whose 64-bit hash is congruent to 0
 * modulo `mod`, giving an expected retained fraction (density) of ~ 1/mod.
 *
 * Selection is CONTENT-KEYED: all identical k-mers share one hash, so every
 * copy of a repeated k-mer is retained or dropped together. This is the
 * control case for the sketching-invariance analysis -- uniqueness among the
 * retained set equals full-k-mer uniqueness (no Poisson thinning of repeats).
 *
 * @param str          sequence
 * @param len          length of sequence
 * @param kmer_size    k-mer size
 * @param mod          keep k-mer iff hash64(kmer) % mod == 0  (density ~ 1/mod)
 * @param reference_id reference ID embedded into the seed
 * @param kmers        output array (malloc'd inside)
 *
 * @note seed encoding matches sketch/kmer.c:
 *       kmer.x = hash << 8 | span
 *       kmer.y = ref  << 32 | pos << 1 | strand
 */
uint64_t sketch_fracmh(const char *str, int len, int kmer_size, uint32_t mod,
                       uint32_t reference_id, uint128_t **kmers);

#ifdef __cplusplus
}
#endif

#endif
