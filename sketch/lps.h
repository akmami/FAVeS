/**
 * @file lps.h
 * @brief Defines the lps struct and its associated methods for handling locally consistent
 * parsing (LCP) of strings.
 *
 * The lps struct is responsible for performing LCP operations on input strings, constructing
 * cores, and supporting functionalities like parsing, compression (deepening). It includes 
 * methods for reading and writing the data to files as well as deepening the LCP to higher 
 * levels of compression.
 *
 * The lps struct leverages various helper classes like core, encoding, and hash to manage the
 * string data and its compressed forms. Additionally, it supports both standard and
 * reverse-complement parsing for specialized string handling in bioinformatics and other fields.
 *
 * Key functionalities include:
 * - Parsing an input string or file to extract LCP cores.
 * - Performing multi-level compression of LCP cores (deepening).
 * - Saving and loading LCP cores from files.
 * - Calculating memory usage of the constructed LCP structure.
 *
 * @namespace lcp
 * @struct lps
 *
 * @note Destructor handles clean-up of allocated memory for cores.
 *
 * @author Akmuhammet Ashyralyyev
 * @version 1.0
 * @date 2024-09-14
 *
 */

#ifndef LPS_H
#define LPS_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h>
#include <stdlib.h>
#include <inttypes.h>
#include <stdint.h> 
#include <string.h>
#include <math.h>

#define CONSTANT_FACTOR         1.5
#define DCT_ITERATION_COUNT     1

#define minimum(a, b) ((a) < (b) ? (a) : (b))

typedef uint32_t ubit_size;
typedef uint32_t ulabel;

#ifdef LARGE_GENOME
struct core {
    ubit_size bit_size;
    ulabel label;
    uint64_t bit_rep;
    uint64_t start;
    uint64_t end;
};
#else 
struct core {
    ubit_size bit_size;
    ulabel label;
    uint64_t bit_rep;
    uint32_t start;
    uint32_t end;
};
#endif

struct lps {
    int level;
    int size;
    struct core *cores;
};

/**
 * @brief Constructs an lps object from a string.
 * 
 * @param lps_ptr The `lps` object that will be initialized
 * @param str The input string to be parsed.
 * @param len The length of the string to be parsed.
 */
void init_lps(struct lps *lps_ptr, const char *str, int len);

/**
 * @brief Constructs an lps object from a string.
 * 
 * @param lps_ptr The `lps` object that will be initialized
 * @param str The input string to be parsed.
 * @param len The length of the string to be parsed.
 * @param offset The length of the offset in which each index will be shifted.
 */
void init_lps_offset(struct lps *lps_ptr, const char *str, int len, uint64_t offset);

/**
 * @brief Constructs an lps object from a string, with reverse complement
 * transformation.
 * 
 * @param lps_ptr The `lps` object that will be initialized
 * @param str The input string to be parsed.
 * @param len The length of the string to be parsed.
 */
void init_lps2(struct lps *lps_ptr, const char *str, int len);
/**
 * @brief Initializes an lps object by reading its contents from a binary file.
 *
 * This function reads the level and size of the `lps` object from the provided file, 
 * allocates memory for the cores array if necessary, and then reads each core's data. 
 * If any error occurs during reading, it will print an error message and terminate 
 * the program or return control, depending on error-handling strategy.
 *
 * @param lps_ptr The `lps` object that will be initialized
 * @param in File pointer to the binary file containing the serialized lps data.
 * 
 * @note The caller must ensure that the `FILE *in` is a valid and open binary file 
 *       for reading. The function will allocate memory for `lps_ptr->cores` if 
 *       `lps_ptr->size > 0`. The caller is responsible for freeing this memory later.
 */
void init_lps3(struct lps *lps_ptr, FILE *in);

/**
 * @brief Constructs an lps object from a string, using split and merge paradigm.
 * The give string will be roughly split into the length of `chunk_size`
 * and processed independently. The final cores will be merged into a single array.
 * 
 * @param lps_ptr The `lps` object that will be initialized
 * @param str The input string to be parsed.
 * @param len The length of the string to be parsed.
 * @param chunk_size The length of the chunks to be processed.
 */
void init_lps4(struct lps *lps_ptr, const char *str, int len, int lcp_level, int chunk_size);

/**
 * @brief Destructor for the lps object. Frees dynamically allocated memory for cores.
 */
void free_lps(struct lps *lps_ptr);

/**
 * @brief Calculates and returns the memory size used by the `lps` structure.
 *
 * @return The memory size (in bytes) used by the `lps` structure.
 */
int64_t lps_memsize(const struct lps *lps_ptr);

/**
 * @brief Deepens the compression level of the LPS structure by one level using a custom DCT iteration count.
 *
 * This function performs DCT compression on the existing cores, parses the remaining cores to find
 * the cores of the next level, removes obsolete cores, updates the size of the LPS structure, and
 * increments the current compression level.
 *
 * The `dct_iteration_count` parameter controls how many DCT iterations are applied before parsing
 * the remaining cores into the next level.
 *
 * @param lps_ptr The `lps` object to deepen.
 * @param dct_iteration_count The number of DCT iterations to perform.
 * @return 1 if the structure was successfully deepened, 0 otherwise.
 */
int lps_deepen1_dct_iters(struct lps *lps_ptr, int dct_iteration_count);

/**
 * @brief Deepens the compression level of the LPS structure by one level using the default DCT iteration count.
 *
 * This function compresses the existing cores and finds new cores for the next compression level
 * using `DCT_ITERATION_COUNT`. It is a convenience wrapper around `lps_deepen1_dct_iters`.
 *
 * @param lps_ptr The `lps` object to deepen.
 * @return 1 if the structure was successfully deepened, 0 otherwise.
 */
int lps_deepen1(struct lps *lps_ptr);

/**
 * @brief Deepens the compression level of the LPS structure to a specific level using a custom DCT iteration count.
 *
 * This function repeatedly deepens the LPS structure one level at a time until the requested
 * `lcp_level` is reached or no further deepening is possible. Each deepening step performs DCT
 * compression using the provided `dct_iteration_count`, then parses the remaining cores into the
 * next compression level.
 *
 * If the requested level is less than or equal to the current level, no deepening is performed.
 *
 * @param lps_ptr The `lps` object to deepen.
 * @param lcp_level The target compression level to deepen to.
 * @param dct_iteration_count The number of DCT iterations to perform at each level.
 * @return 1 if deepening was performed or attempted, 0 if the target level is not greater than the current level.
 */
int lps_deepen_dct_iters(struct lps *lps_ptr, int lcp_level, int dct_iteration_count);

/**
 * @brief Deepens the compression level of the LPS structure to a specific level using the default DCT iteration count.
 *
 * This function repeatedly deepens the LPS structure one level at a time until the requested
 * `lcp_level` is reached or no further deepening is possible. It uses `DCT_ITERATION_COUNT` for
 * each DCT compression step and is a convenience wrapper around `lps_deepen_dct_iters`.
 *
 * @param lps_ptr The `lps` object to deepen.
 * @param lcp_level The target compression level to deepen to.
 * @return 1 if deepening was performed or attempted, 0 if the target level is not greater than the current level.
 */
int lps_deepen(struct lps *lps_ptr, int lcp_level);

/**
 * @brief Outputs the representation of a `lcp` pointer.
 *
 * This function iterates over the cores in the `cores` array and outputs them.
 * It prints each cores as bit representation, after initially printing the LCP level.
 *
 * @param lps_ptr A pointer to the `lps` object to be output.
 */
void print_lps(const struct lps *lps_ptr);

#ifdef __cplusplus
}
#endif

#endif