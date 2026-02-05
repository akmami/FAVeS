#include <chrono>
#include <fstream>
#include <iostream>
#include <vector>
#include <unordered_map>
#include <algorithm>
#include "../sketch/blend.h"
#include "../utils/bvec.h"
#include "utils.h"


#define RELAXATION_SPAN 5

BVEC_INIT(fseeds, uint128_t)

#define BLEND_NEIGHBOR_NUMBER 5

void analyze(uint128_t *fuzzy_seeds, uint64_t len,
             int &contiguous_counts,
             int (&lengths)[DISTANCE_LENGTH],
             std::vector<uint64_t> &lengthsXL,
			 gap &length_gaps) {

	if (len > 0) {

		bool isOverlapped = false;

		if (__blend_get_length(fuzzy_seeds[0]) < DISTANCE_LENGTH) {
			lengths[__blend_get_length(fuzzy_seeds[0])] += 1;
		} else {
			lengthsXL.push_back(__blend_get_length(fuzzy_seeds[0]));
		}

		for (uint64_t i = 1; i < len; i++) {

			uint64_t current_start = __blend_get_index(fuzzy_seeds[i]);
			uint64_t current_end = current_start + __blend_get_length(fuzzy_seeds[i]);
			uint64_t previous_start = __blend_get_index(fuzzy_seeds[i-1]);
			uint64_t previous_end = previous_start + __blend_get_length(fuzzy_seeds[i-1]);

			if (current_start <= previous_end) {
				contiguous_counts += 1;
			}

			// Length
			uint64_t length = current_end - current_start;
			if (length < DISTANCE_LENGTH) {
				lengths[length] += 1;
			} else {
				lengthsXL.push_back(length);
			}

			length_gaps.minimum = std::min(length_gaps.minimum, length);
			length_gaps.maximum = std::max(length_gaps.maximum, length);
		}

		if (isOverlapped) {
			contiguous_counts += 1;
		}
	}
};

void process(std::string &sequence,
             int window, int kmer_size, int blend_bits,
             int &core_counts,
             int &contiguous_counts,
             std::unordered_map<uint32_t, size_t> &distinct_cores,
             std::chrono::milliseconds &durations,
             int (&lengths)[DISTANCE_LENGTH],
             std::vector<uint64_t> &lengthsXL,
			 gap &length_gaps,
             double &sizes,
			 fseeds_bvec_t **seeds, int chr_idx) {

	auto start = std::chrono::high_resolution_clock::now();

    uint128_t *fuzzy_seeds;
    uint64_t fuzzy_seeds_len = 0;
    
    fuzzy_seeds_len = blend_sb_sketch(sequence.c_str(), sequence.length(), window, kmer_size, blend_bits, BLEND_NEIGHBOR_NUMBER, 0, &fuzzy_seeds);

	auto extraction_end = std::chrono::high_resolution_clock::now();
	durations += std::chrono::milliseconds(std::chrono::duration_cast<std::chrono::milliseconds>(extraction_end - start).count());
	core_counts += fuzzy_seeds_len;
    sizes += sizeof(uint128_t) * fuzzy_seeds_len;
	
	analyze(fuzzy_seeds, fuzzy_seeds_len, contiguous_counts, lengths, lengthsXL, length_gaps);

    for (uint64_t index = 0; index < fuzzy_seeds_len; index++) {
        distinct_cores[__blend_get_kmer(fuzzy_seeds[index])]++;
    }

	sequence.clear();

	if (fuzzy_seeds_len) {
		seeds[chr_idx] = fseeds_init(fuzzy_seeds_len);
		if (!seeds[chr_idx]) {
			printf("Coulnd't allocate array\n");
			exit(1);
		}
		memcpy(seeds[chr_idx]->array, fuzzy_seeds, sizeof(uint128_t) * fuzzy_seeds_len);
		seeds[chr_idx]->size = fuzzy_seeds_len;
		free(fuzzy_seeds);
	} else {
		seeds[chr_idx] = NULL;
	}
}

static inline int64_t clamp_i64(int64_t x, int64_t lo, int64_t hi) {
    if (x < lo) return lo;
    if (x > hi) return hi;
    return x;
}

void relaxed_uniqueness_paint(uint128_t *fuzzy_seeds, uint64_t fuzzy_seeds_len, std::unordered_map<uint32_t, size_t> &distinct_cores, uint64_t (&relaxed_unique_count)[RELAXATION_SPAN], uint64_t (&total_relaxed_span)[RELAXATION_SPAN]) {

    // Precompute which positions are truly unique
    uint8_t *is_unique = (uint8_t*)malloc(fuzzy_seeds_len);
    if (!is_unique) abort();
	
    for (uint64_t i = 0; i < fuzzy_seeds_len; i++) {
		auto it = distinct_cores.find(__blend_get_kmer(fuzzy_seeds[i]));
		is_unique[i] = (it != distinct_cores.end() && it->second == 1) ? 1 : 0;
    }

    for (int r = 1; r <= RELAXATION_SPAN; r++) {
        uint8_t *relaxed = (uint8_t*)calloc(fuzzy_seeds_len, 1);
        if (!relaxed) abort();

        for (uint64_t i = 0; i < fuzzy_seeds_len; i++) {
            if (!is_unique[i]) continue;

			
            int64_t left_index = clamp_i64((int64_t)i - r, 0, (int64_t)fuzzy_seeds_len - 1);
            int64_t right_index = clamp_i64((int64_t)i + r, 0, (int64_t)fuzzy_seeds_len - 1);

            for (int64_t j = left_index; j <= right_index; j++) {
                relaxed[j] = 1;
            }
        }

        // count ones and calculate span
        uint64_t ones = 0;
        for (uint64_t i = 0; i < fuzzy_seeds_len; i++) {
            ones += relaxed[i];
        }
		std::vector<std::pair<uint64_t,uint64_t>> intervals;
		intervals.reserve(ones);

		for (uint64_t i = 0; i < fuzzy_seeds_len; i++) {
			if (!relaxed[i]) continue;
			uint64_t s = __blend_get_index(fuzzy_seeds[i]);
			uint64_t e = s + __blend_get_length(fuzzy_seeds[i]);
			intervals.emplace_back(s, e);
		}

		std::sort(intervals.begin(), intervals.end());

		uint64_t total_span = 0;
		uint64_t span_start = intervals[0].first;
		uint64_t span_end   = intervals[0].second;

		for (size_t k = 1; k < intervals.size(); k++) {
			auto [s, e] = intervals[k];
			if (s <= span_end) {
				span_end = std::max(span_end, e);
			} else {
				total_span += span_end - span_start;
				span_start = s;
				span_end   = e;
			}
		}
		total_span += span_end - span_start;

        relaxed_unique_count[r-1] += ones;
		total_relaxed_span[r-1] += total_span;
        free(relaxed);
    }

    free(is_unique);
}

int main(int argc, char **argv) {

	if (argc < 5) {
		std::cerr << "Wrong format: " << argv[0] << " [infile] [kmer-size] [window] [blend-bits]" << std::endl;
		return -1;
	}

	std::ifstream input(argv[1]);
	if (!input.good()) {
		std::cerr << "Error opening: " << argv[1] << " . You have failed." << std::endl;
		return -1;
	}

    int kmer_size = atoi(argv[2]);
    int window = atoi(argv[3]);
	int blend_bits = atoi(argv[4]);

	// variables
	std::string line;

	std::fstream genome;
	genome.open(argv[1], std::ios::in);

	uint64_t total_genone_length = 0;

    // section 1
    int core_counts = 0;
	int contiguous_counts = 0;
    std::unordered_map<uint32_t, size_t> distinct_cores;
    std::chrono::milliseconds durations;
	// section 2
	int lengths[DISTANCE_LENGTH] = {0};
    std::vector<uint64_t> lengthsXL;
	gap length_gaps;
    length_gaps.minimum = UINT64_MAX;
    length_gaps.maximum = 0;
    // section 3
    double sizes = 0;

    // other
	size_t count_once = 0;

	// fuzzy seeds related variables
	fseeds_bvec_t **seeds;
	uint64_t seeds_len = 0, seeds_cap = 128;
	seeds = (fseeds_bvec_t **)malloc(sizeof(fseeds_bvec_t *) * seeds_cap);
	int chr_idx = 0;

	// read file
	if (genome.is_open()) {

		// process genome
		std::string sequence, id;
		sequence.reserve(250000000);

		while (getline(genome, line)) {

			if (line[0] == '>') {
				// process previous chromosome before moving into new one
				if (sequence.size() != 0) {
					if (seeds_len == seeds_cap) {
						seeds_cap *= 2;
						fseeds_bvec_t **temp = (fseeds_bvec_t **)realloc(seeds, sizeof(fseeds_bvec_t *) * seeds_cap);
						if (!temp) {
							printf("realloc failed\n");
							exit(1);
						}
						seeds = temp;
					}
					total_genone_length += sequence.size() - std::count(sequence.begin(), sequence.end(), 'N');
					process(sequence, window, kmer_size, blend_bits, core_counts, contiguous_counts, distinct_cores, durations, lengths, lengthsXL, length_gaps, sizes, seeds, chr_idx);
					chr_idx++;
				}

				id = line.substr(1);
				// std::cout << "Processing started for " << id << std::endl;
				continue;

			} else if (line[0] != '>') {
				sequence += line;
			}
		}

		if (sequence.size() != 0) {
			if (seeds_len == seeds_cap) {
				seeds_cap *= 2;
				fseeds_bvec_t **temp = (fseeds_bvec_t **)realloc(seeds, sizeof(fseeds_bvec_t *) * seeds_cap);
				if (!temp) {
					printf("realloc failed\n");
					exit(1);
				}
				seeds = temp;
			}
			total_genone_length += sequence.size() - std::count(sequence.begin(), sequence.end(), 'N');
			process(sequence, window, kmer_size, blend_bits, core_counts, contiguous_counts, distinct_cores, durations, lengths, lengthsXL, length_gaps, sizes, seeds, chr_idx);
			chr_idx++;
		}

		for (const auto& [value, count] : distinct_cores) {
			if (count == 1) {
				++count_once;
			}
		}

		genome.close();		
	}

	uint64_t relaxed_counts[RELAXATION_SPAN] = {0ULL, 0ULL, 0ULL, 0ULL, 0ULL};
	uint64_t total_relaxed_span[RELAXATION_SPAN] = {0ULL, 0ULL, 0ULL, 0ULL, 0ULL};
	uint64_t total_span = 0;
	// test for uniqueness
	for (int i = 0; i < chr_idx; i++) {
		if (seeds[i] != NULL) {
			std::vector<std::pair<uint64_t,uint64_t>> intervals;
			intervals.reserve(seeds[i]->size);

			for (uint64_t i = 0; i < seeds[i]->size; i++) {
				auto it = distinct_cores.find(__blend_get_kmer(seeds[i]->array[i]));
				if (it != distinct_cores.end() && it->second == 1) {
					uint64_t s = __blend_get_index(seeds[i]->array[i]);
					uint64_t e = s + __blend_get_length(seeds[i]->array[i]);
					intervals.emplace_back(s, e);
				}
			}

			std::sort(intervals.begin(), intervals.end());

			uint64_t span_start = intervals[0].first;
			uint64_t span_end = intervals[0].second;

			for (size_t k = 1; k < intervals.size(); k++) {
				auto [s, e] = intervals[k];
				if (s <= span_end) {
					span_end = std::max(span_end, e);
				} else {
					total_span += span_end - span_start;
					span_start = s;
					span_end   = e;
				}
			}
			total_span += span_end - span_start;

			relaxed_uniqueness_paint(seeds[i]->array, seeds[i]->size, distinct_cores, relaxed_counts, total_relaxed_span);
			fseeds_free(seeds[i]);
		}
	}

	free(seeds);

	std::cout << std::endl;

	std::setprecision(5);

	// Param Settings
	std::cout << "k: " << kmer_size << ", w: " << window << ", b: " << blend_bits << ", n: " << BLEND_NEIGHBOR_NUMBER << std::endl;

	// Total Cores
	std::cout << "Total \\# Cores: " << format_int(core_counts) << std::endl;

	// Contiguous Cores
	std::cout << "Contiguous Cores: " << format_int(contiguous_counts) << std::endl;

	// Total Span
	std::cout << "Total Span Ratio %: " << ((double)total_span) / total_genone_length << std::endl;

	// Distinct Cores
    std::cout << "Distinct Cores: " << format_int(distinct_cores.size()) << std::endl;

	// Unique Cores
    std::cout << "Unique Cores: " << format_int(count_once) << std::endl;

	// Uniqueness Ratio
    std::cout << "Uniqueness %: " << format_double(((double)count_once) / ((double)core_counts), 5) << std::endl;

	// Relaxed Uniqueness Ratio
	for (int i = 0; i < RELAXATION_SPAN; i++) {
		std::cout << "Relaxed-Uniqueness (-" << i + 1 << ",+" << i + 1 << ") %: " << format_double(((double)relaxed_counts[i]) / ((double)core_counts), 5) << ", span: " << ((double)total_relaxed_span[i]) / total_genone_length << std::endl;
	}
	
	// Execution Time
	std::cout << "Exec. Time (sec): " << format_double(((double)durations.count()) / 1000) << std::endl;

	// Mean Core Length
	std::cout << "Avg Length: " << format_double(mean(lengths, lengthsXL)) << std::endl;

	// Std Dev of Lengths
	std::cout << "StdDev Length: " << format_double(stdev(lengths, lengthsXL)) << std::endl;

	// Min/Max Core Length
	std::cout << "Min/Max Length: " << format_int(length_gaps.minimum) << "/" << format_int(length_gaps.maximum) << std::endl;

	// Total Sizes
	std::cout << "Total Size (GB): " << format_double(sizes / (1024.0 * 1024.0 * 1024.0)) << std::endl << std::endl;

	return 0;
};