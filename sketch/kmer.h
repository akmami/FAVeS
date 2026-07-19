#ifndef __KMERS__
#define __KMERS__

#ifdef __cplusplus
extern "C" {
#endif

#include <stdlib.h>
#include <stdint.h>
#include <assert.h>
#include <string.h>

#define _get_kmer(seed) 		((seed).x >> 8)
#define _get_length(seed) 		((seed).x & 0xFF)
#define _get_reference_id(seed) (((seed).y >> 32))
#define _get_index(seed) 		(((seed).y & 0xFFFFFFFF) >> 1)
#define _get_strand(seed) 		(((seed).y & 1))

#ifndef __UINT128_T__
#define __UINT128_T__
typedef struct {
	uint64_t x; // kmer(56) + span(8)
	uint64_t y; // reference_id(32) + index(31) + strand(1)
} uint128_t;
#endif


/**
 * Find all kmers on a DNA sequence
 *
 * @param str    		 sequence
 * @param len    		 length of sequence
 * @param kmer_size      k-mer size
 * @param reference_id   reference ID; will be embeded to uint128_t struct
 * @param kmers          the array of kmers being processed and found
 * 
 * @note kmer is encoded as:
 * 		kmer.x = hash << 8 | span
 * 		kmer.y = ref << 32 | pos << 1 | strand 
 */
uint64_t sketch_kmers(const char *str, int len, int kmer_size, uint32_t reference_id, uint128_t **kmers);

#ifdef __cplusplus
}
#endif

#endif