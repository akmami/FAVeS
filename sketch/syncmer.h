#ifndef __SYNCMERS__
#define __SYNCMERS__

#ifdef __cplusplus
extern "C" {
#endif

#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <limits.h>

#define __syncmer_get_kmer(seed) 			((seed).x >> 8)
#define __syncmer_get_length(seed) 			((seed).x & 0xFF)
#define __syncmer_get_reference_id(seed)  	(((seed).y >> 32))
#define __syncmer_get_index(seed) 			(((seed).y & 0xFFFFFFFF) >> 1)
#define __syncmer_get_strand(seed) 			(((seed).y & 1))

#ifndef __UINT128_T__
#define __UINT128_T__
typedef struct {
	uint64_t x; // kmer(56) + span(8)
	uint64_t y; // reference_id(32) + index(31) + strand(1)
} uint128_t;
#endif


/**
 * Find symmetric (s, k)-syncmer on a DNA sequence
 *
 * @param str    		 sequence
 * @param len    		 length of sequence
 * @param smer_size      find a syncmer that has smer either at the beggining or at the end.
 * @param kmer_size      k-mer size
 * @param reference_id   reference ID; will be embeded to syncmer struct
 * @param syncmers       the array of syncmers being processed and found
 * 
 * @note kmer is encoded as:
 * 		syncmer.x = hash << 8 | span
 * 		syncmer.y = ref << 32 | pos << 1 | strand 
 */
uint64_t sketch_syncmers(const char *str, int len, int smer_size, int kmer_size, uint32_t reference_id, uint128_t **syncmers);

#ifdef __cplusplus
}
#endif

#endif