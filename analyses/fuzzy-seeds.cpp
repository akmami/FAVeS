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
#define BLEND_NEIGHBOR_NUMBER 5

BVEC_INIT(fseeds, uint128_t)
BVEC_INIT(uint64, uint64_t)

#define validate_bvec(seeds, seeds_len, seeds_cap) { \
	if (seeds_len == seeds_cap) { \
		seeds_cap *= 2; \
		fseeds_bvec_t **temp = (fseeds_bvec_t **)realloc(seeds, sizeof(fseeds_bvec_t *) * seeds_cap); \
		if (!temp) { \
			printf("realloc failed\n"); \
			exit(1); \
		} \
		seeds = temp; \
	} \
}

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

			// length
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
			 fseeds_bvec_t **seeds, uint64_t chr_idx) {

	auto start = std::chrono::high_resolution_clock::now();

    uint128_t *fuzzy_seeds;
    uint64_t fuzzy_seeds_len = 0;
    
    fuzzy_seeds_len = blend_sb_sketch(sequence.c_str(), sequence.length(), window, kmer_size, blend_bits, BLEND_NEIGHBOR_NUMBER, chr_idx, &fuzzy_seeds);

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
		seeds[chr_idx] = fseeds_init(0);
		if (!seeds[chr_idx]) {
			printf("Coulnd't allocate array\n");
			exit(1);
		}
		seeds[chr_idx]->array = fuzzy_seeds;
		seeds[chr_idx]->size = fuzzy_seeds_len;
		seeds[chr_idx]->capacity = fuzzy_seeds_len;
	} else {
		seeds[chr_idx] = NULL;
		seeds[chr_idx]->size = 0;
	}
}

static inline int64_t clamp_i64(int64_t index, int64_t low, int64_t high) {
    if (index < low) return low;
    if (index > high) return high;
    return index;
}

void relaxed_uniqueness_paint(uint128_t *fuzzy_seeds,
						      uint64_t fuzzy_seeds_len,
    						  const uint8_t *is_unique,
    						  uint64_t (&relaxed_unique_count)[RELAXATION_SPAN],
    						  uint64_t (&total_relaxed_span)[RELAXATION_SPAN],
							  uint64_t (&total_relaxed_gap_span)[RELAXATION_SPAN], uint64_t chr_length,
							  uint64_t (&total_relaxed_gap_count)[RELAXATION_SPAN],
							  gap (&total_relaxed_gap)[RELAXATION_SPAN]) {

    for (int relaxation_span = 1; relaxation_span <= RELAXATION_SPAN; relaxation_span++) {
        uint8_t *relaxed = (uint8_t*)calloc(fuzzy_seeds_len, 1);
        if (!relaxed) abort();

        for (uint64_t seed_index = 0; seed_index < fuzzy_seeds_len; seed_index++) {
            if (!is_unique[seed_index]) continue;
			
            int64_t left_index = clamp_i64((int64_t)seed_index - relaxation_span, 0, (int64_t)fuzzy_seeds_len - 1);
            int64_t right_index = clamp_i64((int64_t)seed_index + relaxation_span, 0, (int64_t)fuzzy_seeds_len - 1);

            for (int64_t relaxed_seed_index = left_index; relaxed_seed_index <= right_index; relaxed_seed_index++) {
                relaxed[relaxed_seed_index] = 1;
            }
        }

		// compute relaxed seeds count
		uint64_t relaxed_seeds_len = 0;
        for (uint64_t seed_index = 0; seed_index < fuzzy_seeds_len; seed_index++) {
            relaxed_seeds_len += relaxed[seed_index];
        }

        // compute span of relaxed seeds
		std::vector<std::pair<uint64_t,uint64_t>> intervals;
		intervals.reserve(relaxed_seeds_len);

		for (uint64_t seed_index= 0; seed_index < fuzzy_seeds_len; seed_index++) {
			if (!relaxed[seed_index]) continue;
			uint64_t s = __blend_get_index(fuzzy_seeds[seed_index]);
			uint64_t e = s + __blend_get_length(fuzzy_seeds[seed_index]);
			intervals.emplace_back(s, e);
		}

		std::sort(intervals.begin(), intervals.end());

		uint64_t total_span = 0;

		if (!intervals.empty()) {
			uint64_t span_start = intervals[0].first;
			uint64_t span_end   = intervals[0].second;

			total_relaxed_gap_span[relaxation_span-1] += span_start;

			total_relaxed_gap[relaxation_span-1].minimum = std::min(total_relaxed_gap[relaxation_span-1].minimum, span_start);
			total_relaxed_gap[relaxation_span-1].maximum = std::max(total_relaxed_gap[relaxation_span-1].maximum, span_start);

			if (span_start) total_relaxed_gap_count[relaxation_span-1]++;

			for (size_t k = 1; k < intervals.size(); k++) {
				auto [s, e] = intervals[k];
				if (s <= span_end) {
					span_end = std::max(span_end, e);
				} else {
					total_span += span_end - span_start;
					total_relaxed_gap_span[relaxation_span-1] += s - span_end;
					total_relaxed_gap_count[relaxation_span-1]++;
					total_relaxed_gap[relaxation_span-1].minimum = std::min(total_relaxed_gap[relaxation_span-1].minimum, s - span_end);
					total_relaxed_gap[relaxation_span-1].maximum = std::max(total_relaxed_gap[relaxation_span-1].maximum, s - span_end);
					span_start = s;
					span_end = e;
				}
			}
			total_span += span_end - span_start;
			total_relaxed_gap_span[relaxation_span-1] += chr_length - span_end;

			if (span_end < chr_length) total_relaxed_gap_count[relaxation_span-1]++;
			total_relaxed_gap[relaxation_span-1].minimum = std::min(total_relaxed_gap[relaxation_span-1].minimum, chr_length - span_end);
			total_relaxed_gap[relaxation_span-1].maximum = std::max(total_relaxed_gap[relaxation_span-1].maximum, chr_length - span_end);
		}
		// fianlize stats
        relaxed_unique_count[relaxation_span-1] += relaxed_seeds_len;
		total_relaxed_span[relaxation_span-1] += total_span;
        
		// cleanup
		free(relaxed);
    }
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
	uint64_t seeds_cap = 128;
	seeds = (fseeds_bvec_t **)malloc(sizeof(fseeds_bvec_t *) * seeds_cap);
	uint64_t chr_idx = 0;
	uint64_bvec_t *chr_lengths = uint64_init(128);

	// read file
	if (genome.is_open()) {

		// process genome
		std::string sequence, id;
		sequence.reserve(250000000);

		while (getline(genome, line)) {

			if (line[0] == '>') {
				// process previous chromosome before moving into new one
				if (sequence.size() != 0) {
					validate_bvec(seeds, chr_idx, seeds_cap)
					total_genone_length += sequence.size() - std::count(sequence.begin(), sequence.end(), 'N');
					uint64_add(chr_lengths, sequence.size());
					process(sequence, window, kmer_size, blend_bits, core_counts, contiguous_counts, distinct_cores, durations, lengths, lengthsXL, length_gaps, sizes, seeds, chr_idx);
					chr_idx++;
				}

				id = line.substr(1);
				continue;

			} else if (line[0] != '>') {
				sequence += line;
			}
		}

		if (sequence.size() != 0) {
			validate_bvec(seeds, chr_idx, seeds_cap)
			total_genone_length += sequence.size() - std::count(sequence.begin(), sequence.end(), 'N');
			uint64_add(chr_lengths, sequence.size());
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

	// compute unique seeds
	uint8_t **is_unique = (uint8_t **)malloc(sizeof(uint8_t *) * chr_idx);
	if (!is_unique) abort();

	for (uint64_t chrom = 0; chrom < chr_idx; chrom++) {
		if (seeds[chrom] == NULL) {
			is_unique[chrom] = NULL;
			continue;
		}

		uint64_t seed_len = seeds[chrom]->size;
		is_unique[chrom] = (uint8_t *)malloc(seed_len);
		if (!is_unique[chrom]) abort();
	}
	for (uint64_t chrom = 0; chrom < chr_idx; chrom++) {
		if (seeds[chrom] == NULL) continue;

		for (uint64_t seed_index = 0; seed_index < seeds[chrom]->size; seed_index++) {
			auto it = distinct_cores.find(__blend_get_kmer(seeds[chrom]->array[seed_index]));
			is_unique[chrom][seed_index] = (it != distinct_cores.end() && it->second == 1) ? 1 : 0;
		}
	}

	// compute relaxations
	uint64_t relaxed_counts[RELAXATION_SPAN] = {0ULL, 0ULL, 0ULL, 0ULL, 0ULL};
	uint64_t total_relaxed_span[RELAXATION_SPAN] = {0ULL, 0ULL, 0ULL, 0ULL, 0ULL};
	uint64_t total_relaxed_gap_span[RELAXATION_SPAN] = {0ULL, 0ULL, 0ULL, 0ULL, 0ULL};
	uint64_t total_relaxed_gap_count[RELAXATION_SPAN] = {0ULL, 0ULL, 0ULL, 0ULL, 0ULL};
	gap total_relaxed_gap[RELAXATION_SPAN] = {{UINT64_MAX, 0}, {UINT64_MAX, 0}, {UINT64_MAX, 0}, {UINT64_MAX, 0}, {UINT64_MAX, 0}};
	uint64_t total_span = 0, total_gap_span = 0, total_gap_count = 0;

	for (uint64_t chrom = 0; chrom < chr_idx; chrom++) {
		if (seeds[chrom] != NULL) {
			std::vector<std::pair<uint64_t,uint64_t>> intervals;
			intervals.reserve(seeds[chrom]->size);

			for (uint64_t seed_index = 0; seed_index < seeds[chrom]->size; seed_index++) {
				if (is_unique[chrom][seed_index]) {
					uint64_t s = __blend_get_index(seeds[chrom]->array[seed_index]);
					uint64_t e = s + __blend_get_length(seeds[chrom]->array[seed_index]);
					intervals.emplace_back(s, e);
				}
			}

			std::sort(intervals.begin(), intervals.end());

			if (!intervals.empty()) {

				uint64_t span_start = intervals[0].first;
				uint64_t span_end   = intervals[0].second;

				total_gap_span += span_start;

				if (span_start) total_gap_count++;

				for (size_t interval_index = 1; interval_index < intervals.size(); interval_index++) {
					auto [s, e] = intervals[interval_index];
					if (s <= span_end) {
						span_end = std::max(span_end, e);
					} else {
						total_span += span_end - span_start;
						total_gap_span += s - span_end;
						total_gap_count++;
						span_start = s;
						span_end = e;
					}
				}
				total_span += span_end - span_start;
				total_gap_span += chr_lengths->array[chrom] - span_end;
				if (span_end < chr_lengths->array[chrom]) total_gap_count++;

			} else {
				std::cerr << "Chromosome " << chrom << " has no unique seeds\n";
			}

			relaxed_uniqueness_paint(seeds[chrom]->array, seeds[chrom]->size, is_unique[chrom], relaxed_counts, total_relaxed_span, total_relaxed_gap_span, chr_lengths->array[chrom], total_relaxed_gap_count, total_relaxed_gap);
			fseeds_free(seeds[chrom]);
		}
	}

	free(seeds);

	for (uint64_t chrom = 0; chrom < chr_idx; chrom++) {
		free(is_unique[chrom]);
	}
	free(is_unique);

	std::cout << std::endl;

	std::setprecision(5);

	// Param Settings
	std::cout << "k: " << kmer_size << ", w: " << window << ", b: " << blend_bits << ", n: " << BLEND_NEIGHBOR_NUMBER << std::endl;
	std::cout << "fasta: " << argv[1] << std::endl;

	// Total Cores
	std::cout << "Total \\# Cores: " << format_int(core_counts) << std::endl;

	// Contiguous Cores
	std::cout << "Contiguous Cores: " << format_int(contiguous_counts) << std::endl;

	// Distinct Cores
    std::cout << "Distinct Cores: " << format_int(distinct_cores.size()) << std::endl;

	// Unique Cores
    std::cout << "Unique Cores: " << format_int(count_once) << std::endl;

	// Uniqueness Ratio
    std::cout << "Uniqueness %: " << format_double(((double)count_once) / ((double)core_counts), 5)  << ", span: " << ((double)total_span) / total_genone_length << ", gap: " << ((double)total_gap_span) / total_genone_length << ", gap #: " << total_gap_count << ", gap mean: " << ((double)total_gap_span) / total_gap_count << std::endl;

	// Relaxed Uniqueness Ratio
	for (int i = 0; i < RELAXATION_SPAN; i++) {
		std::cout << "Relaxed-Uniqueness (-" << i + 1 << ",+" << i + 1 << ") %: " << format_double(((double)relaxed_counts[i]) / ((double)core_counts), 5) << ", span: " << ((double)total_relaxed_span[i]) / total_genone_length<< ", gap %: " << ((double)total_relaxed_gap_span[i]) / total_genone_length << ", gap #: " << total_relaxed_gap_count[i] << ", gap mean: " << ((double)total_relaxed_gap_span[i]) / total_relaxed_gap_count[i] << ", gap: " << total_relaxed_gap[i].minimum << "-" << total_relaxed_gap[i].maximum << std::endl;
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