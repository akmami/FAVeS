#ifndef __MINIMIZERS__
#define __MINIMIZERS__

#ifdef __cplusplus
extern "C" {
#endif

#include <stdlib.h>
#include <stdint.h>
#include <assert.h>
#include <string.h>

#define __minimizer_get_kmer(seed) 		    ((seed).x >> 14)
#define __minimizer_get_length(seed) 		((seed).x & 0x3FFF)
#define __minimizer_get_reference_id(seed)  (((seed).y >> 32))
#define __minimizer_get_index(seed) 		(((seed).y & 0xFFFFFFFF) >> 1)
#define __minimizer_get_strand(seed) 		(((seed).y & 1))

#ifndef __UINT128_T__
#define __UINT128_T__
typedef struct {
	uint64_t x; // kmer(50) + span(14)
	uint64_t y; // reference_id(32) + index(31) + strand(1)
} uint128_t;
#endif


/**
 * Find symmetric (w, k)-minimizers on a DNA sequence
 *
 * @param str    		 sequence
 * @param len    		 length of sequence
 * @param window      	 find a minimizer for every `window` consecutive k-mers
 * @param kmer_size      k-mer size
 * @param only_symmetric include only symmetric kmers
 * @param reference_id   reference ID; will be embeded to minimizer struct
 * 
 * @note kmer is encoded as:
 * 		minimizer.x = hash
 * 		minimizer.y = pos << 32 | len 
 * 		minimizer.z = ref << 1 | strand
 */
uint64_t sketch_minimizers(const char *str, int len, int window, int kmer_size, int only_symmetric, uint32_t reference_id, uint128_t **minimizers);

#ifdef __cplusplus
}
#endif

#endif