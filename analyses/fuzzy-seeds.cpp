#include <chrono>
#include <fstream>
#include <iostream>
#include <array>
#include <vector>
#include <unordered_map>
#include <algorithm>
#include "../sketch/blend.h"
#include "../utils/bvec.h"
#include "utils.h"


#define RELAXATION_SPAN 5

typedef struct {
    uint64_t relaxed_unique_count;
    uint64_t total_relaxed_span;
	uint64_t total_relaxed_gap_span; 
	uint64_t total_relaxed_gap_count;
	gap total_relaxed_gap;
} seed_stats_t;

typedef struct {
	std::string name;
	uint64_t len = 0;
	uint128_t *fuzzy_seeds;
	uint8_t *is_unique;
    uint64_t fuzzy_seeds_len = 0;
} chr_info_t;

BVEC_INIT(chrs, chr_info_t*)

typedef struct {
	uint64_t genome_len = 0;
	uint64_t valid_genome_len = 0;
	chrs_bvec_t *chroms;
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

#define set_min_max(curr_gap, value) { 								\
	(curr_gap).minimum = std::min((curr_gap).minimum, value);	\
	(curr_gap).maximum = std::max((curr_gap).maximum, value);	\
}

void analyze(uint128_t *fuzzy_seeds, uint64_t len,
             int &contiguous_counts,
             int (&lengths)[DISTANCE_LENGTH],
			 gap &length_gaps) {

	if (len > 0) {

		bool isOverlapped = false;

		lengths[__blend_get_length(fuzzy_seeds[0])] += 1;

		for (uint64_t i = 1; i < len; i++) {

			uint64_t previous_start = __blend_get_index(fuzzy_seeds[i-1]);
			uint64_t previous_end = previous_start + __blend_get_length(fuzzy_seeds[i-1]);
			uint64_t current_start = __blend_get_index(fuzzy_seeds[i]);
			uint64_t current_end = current_start + __blend_get_length(fuzzy_seeds[i]);

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

void process(chr_info_t *info, std::string &sequence,
             int window, int kmer_size, int blend_bits, int n_number,
             int &core_counts,
             int &contiguous_counts,
             std::unordered_map<uint32_t, size_t> &distinct_cores,
             std::chrono::milliseconds &durations,
             int (&lengths)[DISTANCE_LENGTH],
			 gap &length_gaps,
             double &sizes,
			 uint64_t chr_idx) {

	auto start = std::chrono::high_resolution_clock::now();

    uint128_t *fuzzy_seeds;
    uint64_t fuzzy_seeds_len = 0;
    
    fuzzy_seeds_len = blend_sb_sketch(sequence.c_str(), sequence.length(), window, kmer_size, blend_bits, n_number, chr_idx, &fuzzy_seeds);

	auto extraction_end = std::chrono::high_resolution_clock::now();
	durations += std::chrono::milliseconds(std::chrono::duration_cast<std::chrono::milliseconds>(extraction_end - start).count());
	core_counts += fuzzy_seeds_len;
    sizes += sizeof(uint128_t) * fuzzy_seeds_len;
	
	analyze(fuzzy_seeds, fuzzy_seeds_len, contiguous_counts, lengths, length_gaps);

    for (uint64_t index = 0; index < fuzzy_seeds_len; index++) {
        distinct_cores[__blend_get_kmer(fuzzy_seeds[index])]++;
    }

	sequence.clear();

	if (fuzzy_seeds_len) {
		info->fuzzy_seeds = fuzzy_seeds;
		info->fuzzy_seeds_len = fuzzy_seeds_len;
		info->is_unique = (uint8_t *)calloc(fuzzy_seeds_len, sizeof(uint8_t));
	} else {
		info->fuzzy_seeds = NULL;
		info->fuzzy_seeds_len = 0;
		info->is_unique = NULL;
	}
}

static inline int64_t clamp_i64(int64_t index, int64_t low, int64_t high) {
    if (index < low) return low;
    if (index > high) return high;
    return index;
}

void relaxed_uniqueness_paint(chr_info_t *info,
							  seed_stats_t (&seed_statistics)[RELAXATION_SPAN+1],
                              std::array<std::ofstream, RELAXATION_SPAN+1> &bed_out, int print_bed) {

	uint64_t chr_length = info->len;
	uint128_t *fuzzy_seeds = info->fuzzy_seeds;
	uint64_t fuzzy_seeds_len = info->fuzzy_seeds_len;
	const uint8_t *is_unique = info->is_unique;

    for (int radius = 0; radius <= RELAXATION_SPAN; radius++) {
        uint8_t *relaxed = (uint8_t*)calloc(fuzzy_seeds_len, 1);
        if (!relaxed) abort();

        for (uint64_t seed_index = 0; seed_index < fuzzy_seeds_len; seed_index++) {
            if (!is_unique[seed_index]) continue;
			
            int64_t left_index = clamp_i64((int64_t)seed_index - radius, 0, (int64_t)fuzzy_seeds_len - 1);
            int64_t right_index = clamp_i64((int64_t)seed_index + radius, 0, (int64_t)fuzzy_seeds_len - 1);

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

int main(int argc, char **argv) {

	if (argc < 6) {
		std::cerr << "Wrong format: " << argv[0] << " [infile] [kmer-size] [window] [blend-bits] [n-number]" << std::endl;
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
	int n_number = atoi(argv[5]);
	int print_bed = 0;

	if (argc == 7) {
		print_bed = 1;
	}

	// variables
	std::string line;

	std::fstream genome;
	genome.open(argv[1], std::ios::in);

	genome_info_t g_info;
	g_info.chroms = chrs_init(128);

    // section 1
    int core_counts = 0;
	int contiguous_counts = 0;
    std::unordered_map<uint32_t, size_t> distinct_cores;
    std::chrono::milliseconds durations;
	// section 2
	int lengths[DISTANCE_LENGTH] = {0};
	gap length_gaps = {UINT64_MAX, 0};
    // section 3
    double sizes = 0;

    // other
	size_t count_once = 0;

	// bed files
	std::array<std::ofstream, RELAXATION_SPAN+1> bed_files;

	if (print_bed) {
		for (int i = 0; i < RELAXATION_SPAN+1; i++) {
			std::string fname = "../out/relaxed_span_" + std::to_string(i) + ".bed";
			bed_files[i].open(fname, std::ios::out);
			if (!bed_files[i].good()) {
				std::cerr << "Error opening BED output: " << fname << "\n";
				return -1;
			}

			// IGV track line
			bed_files[i] << "track name=\"relaxed_span_" << (i)
						<< "\" description=\"Relaxed uniqueness span " << (i)
						<< "\" itemRgb=\"On\"\n";
		}
	}

	// read file
	if (genome.is_open()) {

		// process genome
		std::string sequence, id;
		sequence.reserve(250000000);

		while (getline(genome, line)) {

			if (line[0] == '>') {
				// process previous chromosome before moving into new one
				if (sequence.size() != 0) {
					size_t len = sequence.size();

					chr_info_t *info = (chr_info_t *)malloc(sizeof(chr_info_t));
					info->name = id;
					info->len = len;
					
					g_info.valid_genome_len += len - std::count(sequence.begin(), sequence.end(), 'N');
					g_info.genome_len += len;
					
					process(info, sequence, window, kmer_size, blend_bits, n_number, 
							core_counts, contiguous_counts, distinct_cores, 
							durations, lengths, length_gaps, sizes, 
							g_info.chroms->size);

					chrs_add(g_info.chroms, info);
				}

				id = line.substr(1);
				continue;

			} else if (line[0] != '>') {
				sequence += line;
			}
		}

		if (sequence.size() != 0) {
			size_t len = sequence.size();

			chr_info_t *info = (chr_info_t *)malloc(sizeof(chr_info_t));
			info->name = id;
			info->len = len;

			g_info.valid_genome_len += len - std::count(sequence.begin(), sequence.end(), 'N');
			g_info.genome_len += len;

			process(info, sequence, window, kmer_size, blend_bits, n_number, 
				core_counts, contiguous_counts, distinct_cores, 
				durations, lengths, length_gaps, sizes, 
				g_info.chroms->size);

			chrs_add(g_info.chroms, info);
		}

		for (const auto& [value, count] : distinct_cores) {
			if (count == 1) {
				++count_once;
			}
		}

		genome.close();		
	}

	// compute unique seeds
	{
		uint64_t chr_size = g_info.chroms->size;
		for (uint64_t chrom = 0; chrom < chr_size; chrom++) {
			if (g_info.chroms->array[chrom]->fuzzy_seeds == NULL) continue;
			if (g_info.chroms->array[chrom]->is_unique == NULL) { fprintf(stderr, "should not happend\n"); continue; }

			// good for readability
			uint128_t *fuzzy_seeds = g_info.chroms->array[chrom]->fuzzy_seeds;
			uint64_t seed_len = g_info.chroms->array[chrom]->fuzzy_seeds_len;
			uint8_t *is_unique = g_info.chroms->array[chrom]->is_unique;
			
			for (uint64_t seed_index = 0; seed_index < seed_len; seed_index++) {
				auto it = distinct_cores.find(__blend_get_kmer(fuzzy_seeds[seed_index]));
				is_unique[seed_index] = (it != distinct_cores.end() && it->second == 1) ? 1 : 0;
			}
		}
	}

	// compute relaxations
	seed_stats_t seed_statistics[RELAXATION_SPAN+1];
	memset(seed_statistics, 0, sizeof(seed_stats_t) * (RELAXATION_SPAN + 1));
	for (int i = 0; i <= RELAXATION_SPAN; i++) {
		seed_statistics[i].total_relaxed_gap.minimum = UINT64_MAX;
	}

	{
		uint64_t chr_size = g_info.chroms->size;
		for (uint64_t chrom = 0; chrom < chr_size; chrom++) {
			relaxed_uniqueness_paint(g_info.chroms->array[chrom], seed_statistics, bed_files, print_bed);
		}
	}

	for (uint64_t chrom = 0; chrom < g_info.chroms->size; chrom++) {
		if (g_info.chroms->array[chrom]) {
			if (g_info.chroms->array[chrom]->fuzzy_seeds) free(g_info.chroms->array[chrom]->fuzzy_seeds);
			if (g_info.chroms->array[chrom]->is_unique) free(g_info.chroms->array[chrom]->is_unique);
			free(g_info.chroms->array[chrom]);
			g_info.chroms->array[chrom] = NULL;
		}
	}
	chrs_free(g_info.chroms);

	std::cout << std::endl;

	std::setprecision(5);

	// Param Settings
	std::cout << "k: " << kmer_size << ", w: " << window << ", b: " << blend_bits << ", n: " << n_number << std::endl;
	std::cout << "fasta: " << argv[1] << std::endl;

	// Total Cores
	std::cout << "Total \\# Cores: " << format_int(core_counts) << std::endl;

	// Contiguous Cores
	std::cout << "Contiguous Cores: " << format_int(contiguous_counts) << std::endl;

	// Distinct Cores
    std::cout << "Distinct Cores: " << format_int(distinct_cores.size()) << " - " << format_double((double)distinct_cores.size() / core_counts, 5) << std::endl;

	// Unique Cores
    std::cout << "Unique Cores: " << format_int(count_once) << " - " << format_double((double)count_once / core_counts, 5) << std::endl;

	// Genome Stats
	std::cout << "Genome len: " << g_info.genome_len << 
			", valid regions: " << g_info.valid_genome_len <<  
			", N: " << g_info.genome_len - g_info.valid_genome_len << 
			" - " << (double)(g_info.genome_len - g_info.valid_genome_len) / g_info.genome_len << " % " << std::endl;

	// Relaxed Uniqueness Ratio
	for (int i = 0; i <= RELAXATION_SPAN; i++) {
		std::cout << "Relaxed-Uniqueness (-" << i << ",+" << i << ") %: " << 
			format_double(((double)seed_statistics[i].relaxed_unique_count) / ((double)core_counts), 5) << 
			", span: " << ((double)seed_statistics[i].total_relaxed_span) / g_info.genome_len << 
			", gap %: " << ((double)g_info.valid_genome_len - seed_statistics[i].total_relaxed_span) / g_info.genome_len << 
			", gap #: " << seed_statistics[i].total_relaxed_gap_count << 
			", gap mean: " << ((double)seed_statistics[i].total_relaxed_gap_span) / seed_statistics[i].total_relaxed_gap_count << 
			", gap: " << seed_statistics[i].total_relaxed_gap.minimum << 
			"-" << seed_statistics[i].total_relaxed_gap.maximum << std::endl;
	}
	
	// Execution Time
	std::cout << "Exec. Time (sec): " << format_double(((double)durations.count()) / 1000) << std::endl;

	// Mean Core Length
	std::cout << "Avg Length: " << format_double(mean(lengths)) << std::endl;

	// Std Dev of Lengths
	std::cout << "StdDev Length: " << format_double(stdev(lengths)) << std::endl;

	// Min/Max Core Length
	std::cout << "Min/Max Length: " << format_int(length_gaps.minimum) << "/" << format_int(length_gaps.maximum) << std::endl;

	// Total Sizes
	std::cout << "Total Size (GB): " << format_double(sizes / (1024.0 * 1024.0 * 1024.0)) << std::endl << std::endl;

	return 0;
};