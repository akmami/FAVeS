#ifndef __UTILS_SKETCHES_HPP__
#define __UTİLS_SKETCHES_HPP__

#include <iomanip>
#include <fstream>
#include <iostream>
#include <sstream>
#include <array>
#include <vector>
#include <iostream>
#include <algorithm>
#include <math.h>
#include <stdint.h>
#include "../utils/bvec.h"

#define DISTANCE_LENGTH 65536

#define RELAXATION_SPAN 5

#ifndef __UINT128_T__
#define __UINT128_T__
typedef struct {
	uint64_t x; // kmer(50) + span(14)
	uint64_t y; // reference_id(32) + index(31) + strand(1)
} uint128_t;
#endif


typedef struct {
	uint64_t minimum;
	uint64_t maximum;
} gap;

typedef struct {
    uint64_t relaxed_unique_count;
    uint64_t total_relaxed_span;
	uint64_t total_relaxed_gap_span; 
	uint64_t total_relaxed_gap_count;
	gap      total_relaxed_gap;
} seed_stats_t;

typedef struct {
	std::string     name;
	uint64_t        len = 0;
	uint128_t*      seeds;
	uint8_t*        is_unique;
    uint64_t        seeds_len = 0;
} chr_info_t;

BVEC_INIT(chrs, chr_info_t*)

typedef struct {
	uint64_t     genome_len = 0;
	uint64_t     valid_genome_len = 0;
	chrs_bvec_t* chroms;
} genome_info_t;

#define bed9_write(out, chr, start, end, name, rgb) {	\
	out << chr << '\t' 									\
        << start << '\t'								\
        << end << '\t'									\
        << name << '\t'									\
        << 0 << '\t'									\
        << '+' << '\t'									\
        << start << '\t'								\
        << end << '\t'									\
        << rgb << '\n';									\
}

#define set_min_max(curr_gap, value) { 					        \
	(curr_gap).minimum = std::min((curr_gap).minimum, value);	\
	(curr_gap).maximum = std::max((curr_gap).maximum, value);	\
}

void analyze(uint128_t *seeds, uint64_t len,
             int &contiguous_counts,
             int (&lengths)[DISTANCE_LENGTH],
			 gap &length_gaps) {

	if (len > 0) {

		bool isOverlapped = false;

		lengths[__sketch_get_length(seeds[0])] += 1;

		for (uint64_t i = 1; i < len; i++) {

			uint64_t previous_start = __sketch_get_index(seeds[i-1]);
			uint64_t previous_end = previous_start + __sketch_get_length(seeds[i-1]);
			uint64_t current_start = __sketch_get_index(seeds[i]);
			uint64_t current_end = current_start + __sketch_get_length(seeds[i]);

			if (current_start <= previous_end) {
				contiguous_counts += 1;
				if (!isOverlapped) {
					contiguous_counts += 1;
					isOverlapped = true;
				}
			}

			// length
			uint64_t length = current_end - current_start;
			lengths[length] += 1;

			length_gaps.minimum = std::min(length_gaps.minimum, length);
			length_gaps.maximum = std::max(length_gaps.maximum, length);
		}
	}
};

static inline uint64_t i64_low(uint64_t index, int radius) {
	if (radius <= index) return index - radius;
	return 0;
}

static inline uint64_t i64_high(uint64_t index, int radius, uint64_t end) {
	if (index + radius < end) return index + radius;
	return end - 1;
}

void relaxed_uniqueness_paint(chr_info_t *info,
							  seed_stats_t (&seed_statistics)[RELAXATION_SPAN+1],
                              std::array<std::ofstream, RELAXATION_SPAN+1> &bed_out, int print_bed) {

	uint64_t chr_length = info->len;
	uint128_t *seeds = info->seeds;
	uint64_t seeds_len = info->seeds_len;
	const uint8_t *is_unique = info->is_unique;

    for (int radius = 0; radius <= RELAXATION_SPAN; radius++) {
        uint8_t *relaxed = (uint8_t*)calloc(seeds_len, 1);
        if (!relaxed) abort();

        for (uint64_t seed_index = 0; seed_index < seeds_len; seed_index++) {
            if (!is_unique[seed_index]) continue;
			
            uint64_t left_index = i64_low(seed_index, radius);
            uint64_t right_index = i64_high(seed_index, radius, seeds_len);

            for (uint64_t relaxed_seed_index = left_index; relaxed_seed_index <= right_index; relaxed_seed_index++) {
                relaxed[relaxed_seed_index] = 1;
            }
        }

		// compute relaxed seeds count
		uint64_t relaxed_seeds_len = 0;
        for (uint64_t seed_index = 0; seed_index < seeds_len; seed_index++) {
            relaxed_seeds_len += relaxed[seed_index];
        }

        // compute span of relaxed seeds
		std::vector<std::pair<uint64_t,uint64_t>> intervals;
		intervals.reserve(relaxed_seeds_len);

		for (uint64_t seed_index= 0; seed_index < seeds_len; seed_index++) {
			if (!relaxed[seed_index]) continue;
			uint64_t s = __sketch_get_index(seeds[seed_index]);
			uint64_t e = s + __sketch_get_length(seeds[seed_index]);
			intervals.emplace_back(s, e);
		}

		std::sort(intervals.begin(), intervals.end());

		uint64_t total_span = 0;

		if (!intervals.empty()) {
			std::vector<std::pair<uint64_t,uint64_t>> filled_blocks;
            filled_blocks.reserve(intervals.size());

			uint64_t span_start = intervals[0].first;
			uint64_t span_end   = intervals[0].second;

			seed_statistics[radius].total_relaxed_gap_span += span_start;
			set_min_max(seed_statistics[radius].total_relaxed_gap, span_start);
			if (span_start) seed_statistics[radius].total_relaxed_gap_count++;

			for (size_t k = 1; k < intervals.size(); k++) {
				std::pair<uint64_t,uint64_t> inter = intervals[k];
				if (inter.first <= span_end) {
					span_end = std::max(span_end, inter.second);
				} else {
					filled_blocks.emplace_back(span_start, span_end);
					total_span += span_end - span_start;
					
					uint64_t g = inter.first - span_end;
					seed_statistics[radius].total_relaxed_gap_span += g;
					seed_statistics[radius].total_relaxed_gap_count++;
					set_min_max(seed_statistics[radius].total_relaxed_gap, g);
					
					span_start = inter.first;
					span_end   = inter.second;
				}
			}

			filled_blocks.emplace_back(span_start, span_end);
			total_span += span_end - span_start;

			uint64_t tail_gap = (span_end < chr_length) ? (chr_length - span_end) : 0;
			seed_statistics[radius].total_relaxed_gap_span += tail_gap;
			if (span_end < chr_length) seed_statistics[radius].total_relaxed_gap_count++;
			set_min_max(seed_statistics[radius].total_relaxed_gap, tail_gap);

			if (print_bed) {
				std::ostream &out = bed_out[radius];

				uint64_t prev_end = 0;
				uint64_t filled_i = 0;
				uint64_t gap_i = 0;

				for (const auto &blk : filled_blocks) {
					const uint64_t s = blk.first;
					const uint64_t e = blk.second;

					// gap before block
					if (s > prev_end) {
						// "gap"
						bed9_write(out, info->name, prev_end, s, "gap", "255,0,0");
						gap_i++;
					}

					// "filled"
					bed9_write(out, info->name, s, e, "filled", "0,0,255");
					filled_i++;

					prev_end = e;
				}

				// final gap
				if (prev_end < chr_length) {
					bed9_write(out, info->name, prev_end, chr_length, "gap", "255,0,0");
				}
			}
		} else {
            // no relaxed hits => entire chromosome is a gap
            std::ostream &out = bed_out[radius];
			if (print_bed) {
            	bed9_write(out, info->name, 0, chr_length, "gap", "255,0,0");
			}

            // stats: all gap, no span
			seed_statistics[radius].total_relaxed_gap_span += chr_length;
			seed_statistics[radius].total_relaxed_gap_count++;
			set_min_max(seed_statistics[radius].total_relaxed_gap, chr_length);
        }

		// finalize stats
		seed_statistics[radius].relaxed_unique_count += relaxed_seeds_len;
		seed_statistics[radius].total_relaxed_span += total_span;
        
		// cleanup
		free(relaxed);
    }
}

// ---------------------------------------------------------------
// ---------------------------------------------------------------
// Other utils
// ---------------------------------------------------------------
// ---------------------------------------------------------------

std::string format_int(int value) {
    std::stringstream ss;
    ss.imbue(std::locale(""));
    ss << std::fixed << value;
    return ss.str();
};

std::string format_double(double value, size_t precision = 2) {
    std::stringstream ss;
    ss.imbue(std::locale(""));
    ss << std::fixed << std::setprecision(precision) << value;
    return ss.str();
};

double mean(int (&numbers)[DISTANCE_LENGTH], std::vector<uint64_t> numbersXL = {}) {
    double sum = 0;
    double count = 0;
    for (size_t i = 0; i < DISTANCE_LENGTH; i++) {
        sum += (i * numbers[i]);
        count += numbers[i];
    }
    for (size_t i = 0; i < numbersXL.size(); i++) {
        sum += numbersXL[i];
    }
    count += numbersXL.size();
    return sum / count;
};

double stdev(int (&numbers)[DISTANCE_LENGTH], std::vector<uint64_t> numbersXL = {}) {
    double mean_value = mean(numbers, numbersXL);
    double count = 0;
    for (size_t i = 0; i < DISTANCE_LENGTH; i++) {
        count += numbers[i];
    }
    count += numbersXL.size();
    double variance = 0;
    for (size_t i = 0; i < DISTANCE_LENGTH; i++) {
        variance += ((mean_value - i) * (mean_value - i) * numbers[i]);
    }
    for (size_t i = 0; i < numbersXL.size(); i++) {
        variance += ((mean_value - numbersXL[i]) * (mean_value - numbersXL[i]));
    }
    return sqrt(variance / count);
};


#endif