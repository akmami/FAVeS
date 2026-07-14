#include "lps.h"


static unsigned char alphabet[256] = {
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

// struct Core part
static inline uint32_t rotl32(uint32_t x, int r) {
    return (x << r) | (x >> (32 - r));
}

// Specialization of MurmurHash3_32(data, 16, 42) for four 32-bit words. Bit-identical to the general routine.
static inline uint32_t hash4_u32(uint32_t w0, uint32_t w1, uint32_t w2, uint32_t w3) {
    const uint32_t c1 = 0xcc9e2d51, c2 = 0x1b873593;
    uint32_t h1 = 42u, k1;

    k1 = w0; k1 *= c1; k1 = rotl32(k1,15); k1 *= c2; h1 ^= k1; h1 = rotl32(h1,15); h1 = h1*5 + 0xe6546b64;
    k1 = w1; k1 *= c1; k1 = rotl32(k1,15); k1 *= c2; h1 ^= k1; h1 = rotl32(h1,15); h1 = h1*5 + 0xe6546b64;
    k1 = w2; k1 *= c1; k1 = rotl32(k1,15); k1 *= c2; h1 ^= k1; h1 = rotl32(h1,15); h1 = h1*5 + 0xe6546b64;
    k1 = w3; k1 *= c1; k1 = rotl32(k1,15); k1 *= c2; h1 ^= k1; h1 = rotl32(h1,15); h1 = h1*5 + 0xe6546b64;

    h1 ^= 16u;                 /* len */
    h1 ^= h1 >> 16; h1 *= 0x85ebca6b; h1 ^= h1 >> 13; h1 *= 0xc2b2ae35; h1 ^= h1 >> 16;
    return h1;
}

void init_core1(struct core *cr, const char *begin, uint64_t distance, uint64_t start_index, uint64_t end_index) {
    cr->start = start_index;
    cr->end = end_index;
    cr->label = 0;
    cr->label |= ((distance-2) << 6);
    cr->label |= (alphabet[(int)(*begin)] << 4);
    cr->label |= (alphabet[(int)(*(begin+distance-2))] << 2);
    cr->label |= (alphabet[(int)(*(begin+distance-1))]);
    cr->bit_rep = 0x8000000000000000 | cr->label;
    cr->bit_size = 2 * distance;
}

void init_core2(struct core *cr, const char *begin, uint64_t distance, uint64_t start_index, uint64_t end_index) {
    cr->start = start_index;
    cr->end = end_index;
    cr->label = 0;
    cr->label |= ((distance-2) << 6);
    cr->label |= ((3ULL ^ alphabet[(int)(*(begin))]) << 4);
    cr->label |= ((3ULL ^ alphabet[(int)(*(begin-distance+2))]) << 2);
    cr->label |= ((3ULL ^ alphabet[(int)(*(begin-distance+1))]));
    cr->bit_rep = 0x8000000000000000 | cr->label;
    cr->bit_size = 2 * distance;
}

void init_core3(struct core *cr, struct core *begin, uint64_t distance) {
    cr->start = begin->start;
    cr->end = (begin+distance-1)->end;
    cr->bit_rep = 0;
    cr->bit_size = 0;

    for (struct core *it=begin; it<begin+distance; it++) {
        cr->bit_size += it->bit_size;
    }

    int index = 0;
    for (struct core *it = begin+distance-1; begin <= it && index + it->bit_size <= 64; it--) {
        cr->bit_rep |= (it->bit_rep << index);
        index += it->bit_size;
    }

    cr->bit_rep = 0x7FFFFFFFFFFFFFFF & cr->bit_rep;
    cr->bit_size = minimum(cr->bit_size, 63);
    cr->label = hash4_u32(begin->label, (begin + distance - 2)->label, (begin + distance - 1)->label, (uint32_t)(distance - 2));
}

void init_core4(struct core *cr, ubit_size bit_size, uint64_t bit_rep, ulabel label, uint64_t start, uint64_t end) {
    cr->bit_size = bit_size;
    cr->bit_rep = bit_rep;
    cr->label = label;
    cr->start = start;
    cr->end = end;
}

void print_core(const struct core *cr) {
    if (cr->bit_rep & 0x8000000000000000) { // if printing 1-level cores
        uint64_t middle_count = (0x7FFFFFFFFFFFFFFF & cr->bit_rep) >> 6;
        uint64_t middle_val = (cr->bit_rep >> 2) & 3;
        printf("%" PRIu64 "\n", ((cr->bit_rep >> 5) & 1));
        printf("%" PRIu64 "\n", ((cr->bit_rep >> 4) & 1));
        for (uint64_t i=0; i<middle_count; i++) {
            printf("%" PRIu64 "\n", ((middle_val >> 1) & 1));
            printf("%" PRIu64 "\n", (middle_val & 1));           
        }
        printf("%" PRIu64 "\n", ((cr->bit_rep >> 1) & 1));
        printf("%" PRIu64 "\n", (cr->bit_rep & 1));
    } else {
        for (ubit_size index = cr->bit_size - 1; 0 < index; index--) {
            printf("%" PRIu64 "\n", ((cr->bit_rep >> index) & 1));
        }
        printf("%" PRIu64 "\n", (cr->bit_rep & 1));
    }
}

static inline int core_eq(const struct core *lhs, const struct core *rhs) {
    return lhs->bit_rep == rhs->bit_rep;
}

static inline int core_neq(const struct core *lhs, const struct core *rhs) {
    return lhs->bit_rep != rhs->bit_rep;
}

static inline int core_gt(const struct core *lhs, const struct core *rhs) {
    return lhs->bit_rep > rhs->bit_rep;
}

static inline int core_lt(const struct core *lhs, const struct core *rhs) {
    return lhs->bit_rep < rhs->bit_rep;
}

static inline int core_geq(const struct core *lhs, const struct core *rhs) {
    return lhs->bit_rep >= rhs->bit_rep;
}

static inline int core_leq(const struct core *lhs, const struct core *rhs) {
    return lhs->bit_rep <= rhs->bit_rep;
}

static inline ubit_size umin(ubit_size a, ubit_size b) { 
    return a < b ? a : b; 
}

/**
 * @brief Compute the bit-length of a 64-bit unsigned integer.
 *
 * Bit-length is defined as the position of the most significant set bit
 * (1-based). Returns 0 if x == 0.
 *
 * Uses __builtin_clzll for efficient leading-zero counting.
 */
static inline ubit_size bitlen_u64(uint64_t x) {
    return x ? (ubit_size)(64u - (ubit_size)__builtin_clzll(x)) : 0u;
}

/**
 * @brief Compute the bit-length of x with a minimum result of 2.
 *
 * Ensures the returned bit-length is at least 2, even if x has fewer
 * significant bits (or is zero).
 */
static inline ubit_size bitlen_min2(uint64_t x) {
    ubit_size bl = bitlen_u64(x);
    return bl < 2u ? 2u : bl;
}

/**
 * @brief Extract a 2-bit symbol from the low 6 bits of a packed value.
 *
 * The low 6 bits are interpreted as three 2-bit slots:
 *   k = 0 → rightmost 2 bits
 *   k = 1 → middle 2 bits
 *   k = 2 → leftmost 2 bits
 */
static inline uint64_t sym2(uint64_t rep, unsigned k) {
    return k ? (rep >> k) & 3ull : rep & 3ull;
}

/**
 * @brief Extract the middle repetition count stored above the low 6 bits.
 *
 * Bits 0–5 are reserved for 2-bit symbols. Bits above 6 contain the
 * repetition count. The most significant bit is masked off before shifting.
 */
static inline uint64_t mid_count(uint64_t rep) {
    return (rep & 0x7FFFFFFFFFFFFFFFull) >> 6;
}

/**
 * @brief Emit an encoded index-bit value based on symbol comparison.
 *
 * Computes:
 *   result = i + selected_bit or 2 + i + selected_bit, depending on 
 *  from where the bit is selected
 *
 * If the least significant bits of a2 and b2 differ, select (b2 & 1).
 * Otherwise, select ((b2 >> 1) & 1).
 *
 * @param a2 First 2-bit symbol.
 * @param b2 Second 2-bit symbol.
 * @param i  Base index multiplier.
 * @return Encoded value 2 + i plus chosen bit from b2.
 */
static inline uint64_t emit_idx_bit(uint64_t a2, uint64_t b2, uint64_t i) {
    if ((a2 & 1) != (b2 & 1)) {
        return i + (b2 & 1);
    }
    return 2 + i + ((b2 >> 1) & 1);
}

/**
 * @brief Perform level-1 compression of two adjacent cores.
 *
 * Implements the 3-symbol + middle-count encoding scheme.
 * The caller guarantees that both `left` and `right` are already
 * level-1 encoded.
 */
void core_compress_level1(const struct core *left, struct core *right) {
    uint64_t L = left->bit_rep;
    uint64_t R = right->bit_rep;

    uint64_t L3 = sym2(L, 0), L2 = sym2(L, 2), L1 = sym2(L, 4);
    uint64_t R3 = sym2(R, 0), R2 = sym2(R, 2), R1 = sym2(R, 4);

    uint64_t Lm = mid_count(L);
    uint64_t Rm = mid_count(R);

    uint64_t out;

    if (L3 != R3) {
        // base = 0 (2*i with i=0)
        out = emit_idx_bit(L3, R3, 0);
        right->bit_rep = out;
        right->bit_size = 2;        
    }
    else if (L2 != R2) {
        // Your original used base=4/6 -> i=2 (since 2*i = 4)
        out = emit_idx_bit(L2, R2, 4);
        right->bit_rep = out;
        right->bit_size = bitlen_min2(out);
    }
    else if (Lm != Rm) {
        // Preserve your “compare across boundary depending on which count is smaller”
        if (Lm < Rm) {
            // compare left L1 with right R2; index = Lm + 1
            out = emit_idx_bit(L1, R2, 4 * Lm + 4);
        } else {
            // compare left L2 with right R1; index = Rm + 1
            out = emit_idx_bit(L2, R1, 4 * Rm + 4);
        }
        right->bit_rep = out;
        right->bit_size = bitlen_min2(out);
    }
    else if (L1 != R1) {
        // left mismatch: index = Lm + 1
        out = emit_idx_bit(L1, R1, 4 * Lm + 4);
        right->bit_rep = out;
        right->bit_size = bitlen_min2(out);
    }
    else {
        // same
        out = 2ull * (uint64_t)right->bit_size;
        right->bit_rep = out;
        right->bit_size = bitlen_min2(out);
    }

    right->start = left->start;
}

/**
 * @brief Perform upper-level compression via bounded first-difference encoding.
 *
 * The caller guarantees that both `left` and `right` are already
 * upper-level encoded.
 */
void core_compress_upper(const struct core *left, struct core *right) {
    ubit_size bound = umin(right->bit_size, umin(left->bit_size, 64u));

    ubit_size idx;
    if (left->bit_rep == right->bit_rep) {
        idx = bound;
    } else {
        uint64_t x = left->bit_rep ^ right->bit_rep;   // nonzero here
        idx = (ubit_size)__builtin_ctzll(x);
        idx = umin(idx, bound);
    }

    uint64_t out = 2ull * (uint64_t)idx + ((right->bit_rep >> idx) & 1ull);
    right->bit_rep = out;
    right->bit_size = bitlen_min2(out);
    right->start = left->start;
}

// LCP part

int parse1(const char *begin, const char *end, struct core *cores, uint64_t offset) {

    const char *it1 = begin;
    const char *it2 = end;
    int core_index = 0;
    int last_invalid_char_index = -1;

    // find lcp cores
    for (; it1 + 2 < end; it1++) {

        // skip invalid character
        if (alphabet[(unsigned char)*it1] == 4) {
            last_invalid_char_index = it1 - begin;
            continue;
        }

        if (alphabet[(unsigned char)*it1] == alphabet[(unsigned char)*(it1+1)]) {
            continue;
        }

        // check for RINT core
        if (alphabet[(unsigned char)*(it1+1)] == alphabet[(unsigned char)*(it1+2)]) {

            // count middle characters
            uint32_t middle_count = 1;
            const char *temp = it1 + 2;
            while (temp < end && alphabet[(unsigned char)*(temp-1)] == alphabet[(unsigned char)*temp]) {
                temp++;
                middle_count++;
            }
            if (temp != end) {
                // check if there is any SSEQ cores left behind
                if (it2 < it1 && last_invalid_char_index < it2 - begin - 1) {
                    init_core1(&(cores[core_index]), it2-1, it1-it2+2, it2-begin-1+offset, it1-begin+1+offset);
                    core_index++;
                }

                // create RINT core
                it2 = it1 + 2 + middle_count;
                init_core1(&(cores[core_index]), it1, 2+middle_count, it1-begin+offset, it2-begin+offset);
                core_index++;

                continue;
            }
        }

        if (alphabet[(unsigned char)*it1] > alphabet[(unsigned char)*(it1+1)] &&
            alphabet[(unsigned char)*(it1+1)] < alphabet[(unsigned char)*(it1+2)]) {

            // check if there is any SSEQ cores left behind
            if (it2 < it1 && last_invalid_char_index < it2 - begin - 1) {
                init_core1(&(cores[core_index]), it2-1, it1-it2+2, it2-begin-1+offset, it1-begin+1+offset);
                core_index++;
            }

            // create LMIN core
            it2 = it1 + 3;
            init_core1(&(cores[core_index]), it1, 3, it1-begin+offset, it2-begin+offset);
            core_index++;

            continue;
        }

        if (begin == it1) {
            continue;
        }

        // check for LMAX
        if (it1+3 < end &&
            alphabet[(unsigned char)*it1] < alphabet[(unsigned char)*(it1+1)] &&
            alphabet[(unsigned char)*(it1+1)] > alphabet[(unsigned char)*(it1+2)] &&
            alphabet[(unsigned char)*(it1-1)] <= alphabet[(unsigned char)*(it1)] &&
            alphabet[(unsigned char)*(it1+2)] >= alphabet[(unsigned char)*(it1+3)]) {

            // check if there is any SSEQ cores left behind
            if (it2 < it1 && last_invalid_char_index < it2 - begin - 1) {
                init_core1(&(cores[core_index]), it2-1, it1-it2+2, it2-begin-1+offset, it1-begin+1+offset);
                core_index++;
            }

            // create LMAX core
            it2 = it1 + 3;
            init_core1(&(cores[core_index]), it1, 3, it1-begin+offset, it2-begin+offset);
            core_index++;

            continue;
        }
    }

    return core_index;
}

int parse2(const char *begin, const char *end, struct core *cores, uint64_t offset) {

    const char *it1 = end - 1;
    const char *it2 = begin - 1;
    int core_index = 0;

    // find lcp cores
    for (; begin <= it1 - 2; it1--) {

        // skip invalid character
        if (alphabet[(unsigned char)*it1] == alphabet[(unsigned char)*(it1-1)]) {
            continue;
        }

        // check for RINT core
        if (alphabet[(unsigned char)*(it1-1)] == alphabet[(unsigned char)*(it1-2)]) {

            // count middle characters
            uint32_t middle_count = 1;
            const char *temp = it1 - 2;
            while (begin <= temp && alphabet[(unsigned char)*(temp+1)] == alphabet[(unsigned char)*temp]) {
                temp--;
                middle_count++;
            }
            if (begin <= temp) {
                // check if there is any SSEQ cores left behind
                if (it1 < it2) {
                    init_core2(&(cores[core_index]), it2+1, it2-it1+2, end-it2-1+offset, end-it1-1+offset);
                    core_index++;
                }

                // create RINT core
                it2 = it1 - 2 - middle_count;
                init_core2(&(cores[core_index]), it1, 2+middle_count, end-it1-1+offset, end-it2-1+offset);
                core_index++;

                continue;
            }
        }

        if ((3ULL ^ alphabet[(unsigned char)*it1]) > (3ULL ^ alphabet[(unsigned char)*(it1-1)]) &&
            (3ULL ^ alphabet[(unsigned char)*(it1-1)]) < (3ULL ^ alphabet[(unsigned char)*(it1-2)])) {

            // check if there is any SSEQ cores left behind
            if (it1 < it2) {
                init_core2(&(cores[core_index]), it2+1, it2-it1+2, end-it2-1+offset, end-it1-1+offset);
                core_index++;
            }

            // create LMIN core
            it2 = it1 - 3;
            init_core2(&(cores[core_index]), it1, 3, end-it1-1+offset, end-it2-1+offset);
            core_index++;

            continue;
        }

        if (begin == it1) {
            continue;
        }

        // check for LMAX
        if (begin <= it1-3 &&
            (3ULL ^ alphabet[(unsigned char)*it1]) < (3ULL ^ alphabet[(unsigned char)*(it1-1)]) &&
            (3ULL ^ alphabet[(unsigned char)*(it1-1)]) > (3ULL ^ alphabet[(unsigned char)*(it1-2)]) &&
            (3ULL ^ alphabet[(unsigned char)*(it1+1)]) <= (3ULL ^ alphabet[(unsigned char)*(it1)]) &&
            (3ULL ^ alphabet[(unsigned char)*(it1-2)]) >= (3ULL ^ alphabet[(unsigned char)*(it1-3)])) {

            // check if there is any SSEQ cores left behind
            if (it1 < it2) {
                init_core2(&(cores[core_index]), it2+1, it2-it1+2, end-it2-1+offset, end-it1-1+offset);
                core_index++;
            }

            // create LMAX core
            it2 = it1 - 3;
            init_core2(&(cores[core_index]), it1, 3, end-it1-1+offset, end-it2-1+offset);
            core_index++;

            continue;
        }
    }

    return core_index;
}

int parse3(struct core *begin, struct core *end, struct core *cores) {

    struct core *it1 = begin;
    struct core *it2 = end;
    int core_index = 0;

    // find lcp cores
    for (; it1 + 2 < end; it1++) {

        // skip invalid character
        if (core_eq(it1, it1+1)) {
            continue;
        }

        // check for RINT core
        if (core_eq(it1+1, it1+2)) {

            // count middle characters
            uint32_t middle_count = 1;
            struct core *temp = it1 + 2;
            while (temp < end && core_eq(temp-1, temp)) {
                temp++;
                middle_count++;
            }
            if (temp != end) {
                // check if there is any SSEQ cores left behind
                if (it2 < it1) {
                    init_core3(&(cores[core_index]), it2-1, it1-it2+2);
                    core_index++;
                }

                // create RINT core
                it2 = it1 + 2 + middle_count;
                init_core3(&(cores[core_index]), it1, it2-it1);
                core_index++;

                continue;
            }
        }

        // check for LMIN
        if (core_gt(it1, it1+1) && core_lt(it1+1, it1+2)) {
            
            // check if there is any SSEQ cores left behind
            if (it2 < it1) {
                init_core3(&(cores[core_index]), it2-1, it1-it2+2);
                core_index++;
            }

            // create LMIN core
            it2 = it1 + 3;
            init_core3(&(cores[core_index]), it1, it2-it1);
            core_index++;

            continue;
        }

        if (begin == it1) {
            continue;
        }

        // check for LMAX
        if (it1+3 < end &&
            core_lt(it1, it1+1) &&
            core_gt(it1+1, it1+2) &&
            core_leq(it1-1, it1) &&
            core_geq(it1+2, it1+3)) {

            // check if there is any SSEQ cores left behind
            if (it2 < it1) {
                init_core3(&(cores[core_index]), it2-1, it1-it2+2);
                core_index++;
            }

            // create LMAX core
            it2 = it1 + 3;
            init_core3(&(cores[core_index]), it1, it2-it1);
            core_index++;

            continue;
        }
    }
    return core_index;
}

/**
 * @brief Performs Deterministic Coin Tossing (DCT) compression on binary sequences.
 *
 * This function is a central part of the LCP (Locally Consisted Parsing) algorithm. It identifies differences
 * between consecutive binary strings, compressing the information by focusing on the position and value of
 * the first divergent bit from the right-end of the strings. This difference is used to generate a compact
 * 'core' that encapsulates the unique elements of each sequence in a smaller binary form.
 *
 * This compression significantly reduces redundant information, making further analysis of the sequences
 * within the LCP framework more efficient and manageable.
 *
 * @return 0 if dct is performed, -1 if no enough cores are available for dct.
 */
int lcp_dct(struct lps *lps_ptr) {

    // at least 2 cores are needed for compression
    if (lps_ptr->size < DCT_ITERATION_COUNT + 1) {
        return -1;
    }

    if (lps_ptr->level == 1) {
        for (uint64_t dct_index = 0; dct_index < DCT_ITERATION_COUNT; dct_index++) {
            struct core *it_left = lps_ptr->cores + lps_ptr->size - 2, *it_right = lps_ptr->cores + lps_ptr->size - 1;

            for (; lps_ptr->cores + dct_index <= it_left; it_left--, it_right--) {
                core_compress_level1(it_left, it_right);
            }
        }
    } else {
        for (uint64_t dct_index = 0; dct_index < DCT_ITERATION_COUNT; dct_index++) {
            struct core *it_left = lps_ptr->cores + lps_ptr->size - 2, *it_right = lps_ptr->cores + lps_ptr->size - 1;

            for (; lps_ptr->cores + dct_index <= it_left; it_left--, it_right--) {
                core_compress_upper(it_left, it_right);
            }
        }
    }

    return 0;
}

int lps_deepen1_dct_iters(struct lps *lps_ptr, int dct_iteration_count) {
    // compress cores
    if (lcp_dct(lps_ptr) < 0) {
        lps_ptr->size = 0;
        lps_ptr->level++;
        return 0;
    }

    // find new cores
    int new_size = parse3(lps_ptr->cores + dct_iteration_count, lps_ptr->cores + lps_ptr->size, lps_ptr->cores);
    int temp = new_size;

    // remove old cores
    while(temp < lps_ptr->size) {
        temp++;
    }
    lps_ptr->size = new_size;

    lps_ptr->level++;

    if (lps_ptr->size)
        lps_ptr->cores = (struct core*)realloc(lps_ptr->cores, lps_ptr->size * sizeof(struct core));

    return 1;
}

int lps_deepen1(struct lps *lps_ptr) {
    return lps_deepen1_dct_iters(lps_ptr, DCT_ITERATION_COUNT);
}

int lps_deepen_dct_iters(struct lps *lps_ptr, int lcp_level, int dct_iteration_count) {

    if (lcp_level <= lps_ptr->level)
        return 0;

    while (lps_ptr->level < lcp_level && lps_deepen1_dct_iters(lps_ptr, dct_iteration_count))
        ;

    return 1;
}

int lps_deepen(struct lps *lps_ptr, int lcp_level) {
    return lps_deepen_dct_iters(lps_ptr, lcp_level, DCT_ITERATION_COUNT);
}

void init_lps(struct lps *lps_ptr, const char *str, int len) {   
    lps_ptr->level = 1;
    lps_ptr->size = 0;
    lps_ptr->cores = (struct core *)malloc((len/CONSTANT_FACTOR)*sizeof(struct core));
    lps_ptr->size = parse1(str, str+len, lps_ptr->cores, 0);
}

void init_lps_offset(struct lps *lps_ptr, const char *str, int len, uint64_t offset) {   
    lps_ptr->level = 1;
    lps_ptr->size = 0;
    lps_ptr->cores = (struct core *)malloc((len/CONSTANT_FACTOR)*sizeof(struct core));
    lps_ptr->size = parse1(str, str+len, lps_ptr->cores, offset);
}

void init_lps2(struct lps *lps_ptr, const char *str, int len) {   
    lps_ptr->level = 1;
    lps_ptr->size = 0;
    lps_ptr->cores = (struct core *)malloc((len/CONSTANT_FACTOR)*sizeof(struct core));
    lps_ptr->size = parse2(str, str+len, lps_ptr->cores, 0);
}

void init_lps3(struct lps *lps_ptr, FILE *in) {
    // read the level from the binary file
    if (fread(&(lps_ptr->level), sizeof(int), 1, in) != 1) {
        fprintf(stderr, "Error reading level from file\n");
        exit(EXIT_FAILURE);
    }

    // read the size (number of cores)
    if(fread(&(lps_ptr->size), sizeof(int), 1, in) != 1) {
        fprintf(stderr, "Error reading size from file\n");
        exit(EXIT_FAILURE);
    }

    lps_ptr->cores = NULL;

    if (lps_ptr->size) {
        // allocate memory for the cores array
        lps_ptr->cores = (struct core *)malloc(lps_ptr->size * sizeof(struct core));
        if (fread(lps_ptr->cores, lps_ptr->size * sizeof(struct core), 1, in) != 1) {
            fprintf(stderr, "Error reading cores from file\n");
            exit(EXIT_FAILURE);
        }
    }
}

void init_lps4(struct lps *lps_ptr, const char *str, int len, int lcp_level, int chunk_size) {

    if (lcp_level < 1)
        return;

    lps_ptr->level = 1;
    lps_ptr->size = 0; 
    int estimated_size = (int)(len / pow((double)CONSTANT_FACTOR, lcp_level));
    lps_ptr->cores = (struct core *)malloc(estimated_size*sizeof(struct core));

    int str_index = 0, core_index = 0;

    {
        int str_len = minimum(chunk_size, len);
        struct lps temp_lps;
        init_lps_offset(&temp_lps, str, str_len, 0);
        lps_deepen(&temp_lps, lcp_level);

        if (temp_lps.size) {
            memcpy(lps_ptr->cores, temp_lps.cores, (temp_lps.size)*sizeof(struct core));
            core_index = (temp_lps.size);
            lps_ptr->size = (temp_lps.size);
            if (temp_lps.size>1)
                str_index = lps_ptr->cores[core_index-2].start;
            else 
                str_index = lps_ptr->cores[core_index-1].start;
        }
        free(temp_lps.cores);
    }

    while (str_index < len) {
        int str_len = minimum(chunk_size, len-str_index);
        struct lps temp_lps;
        init_lps_offset(&temp_lps, str+str_index, str_len, str_index);
        lps_deepen(&temp_lps, lcp_level);

        if (1<temp_lps.size) {
            int overlap = 2;
            while (0<overlap) {
                if (lps_ptr->cores[core_index-overlap].start == temp_lps.cores[0].start)
                    break;
                overlap--;
            }
            memcpy(lps_ptr->cores+core_index, temp_lps.cores+overlap, (temp_lps.size-overlap)*sizeof(struct core));
            core_index += (temp_lps.size-overlap);
            lps_ptr->size += (temp_lps.size-overlap);

            if ((uint64_t)str_index < lps_ptr->cores[core_index-2].start) {
                str_index = lps_ptr->cores[core_index-2].start;
                free(temp_lps.cores);
                continue;
            } 
        }
        
        // find next start point
        for(int i=str_index+str_len-1; str_index <= i; i--) {
            if (alphabet[(unsigned char)*(str+i)] == 4) {
                str_index = i+1;
                break;
            }
        }
        if (alphabet[(unsigned char)*(str+str_index)] != 4) { // all of the characters are valid, so not valid cores found
            str_index += str_len;
        }
        
        free(temp_lps.cores);
    }

    if (lps_ptr->size)
        lps_ptr->cores = (struct core*)realloc(lps_ptr->cores, lps_ptr->size * sizeof(struct core));
}

void free_lps(struct lps *lps_ptr) {
    free(lps_ptr->cores);
    lps_ptr->size = 0;
}

int64_t lps_memsize(const struct lps *lps_ptr) {
    return sizeof(struct lps) + lps_ptr->size * sizeof(struct core);
}

void print_lps(const struct lps *lps_ptr) {
    printf("Level: %d \n", lps_ptr->level);
    for(int i=0; i<lps_ptr->size; i++) {
        print_core(&(lps_ptr->cores[i]));
        printf(" ");
    }
}

int lps_eq(const struct lps *lhs, const struct lps *rhs) {
    if (lhs->size != rhs->size) {
        return 0;
    }

    for(int i=0; i<lhs->size; i++) {
        if (core_neq(&(lhs->cores[i]), &(rhs->cores[i])) != 0) {
            return 0;
        }
    }

    return 1;
}

int lps_neq(const struct lps *lhs, const struct lps *rhs) {
    if (lhs->size != rhs->size) {
        return 1;
    }

    for(int i=0; i<lhs->size; i++) {
        if (core_neq(&(lhs->cores[i]), &(rhs->cores[i])) != 0) {
            return 1;
        }
    }

    return 0;
}