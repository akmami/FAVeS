#include <chrono>
#include <unordered_map>
#include <algorithm>
#include <string>
#include "../sketch/blend.h"


#define __sketch_get_kmer(kmer) __blend_get_kmer(kmer)
#define __sketch_get_length(kmer) __blend_get_length(kmer)
#define __sketch_get_reference_id(kmer) __blend_get_reference_id(kmer)
#define __sketch_get_index(kmer) __blend_get_index(kmer)
#define __sketch_get_strand(kmer) __blend_get_strand(kmer)


// this to be included after sketch macro definitions
#include "utils.hpp"


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
    
    fuzzy_seeds_len = sketch_blend(sequence.c_str(), sequence.length(), window, kmer_size, blend_bits, n_number, chr_idx, &fuzzy_seeds);

	auto extraction_end = std::chrono::high_resolution_clock::now();
	durations += std::chrono::milliseconds(std::chrono::duration_cast<std::chrono::milliseconds>(extraction_end - start).count());
	core_counts += fuzzy_seeds_len;
    sizes += sizeof(uint128_t) * fuzzy_seeds_len;
	
	analyze(fuzzy_seeds, fuzzy_seeds_len, contiguous_counts, lengths, length_gaps);

    for (uint64_t index = 0; index < fuzzy_seeds_len; index++) {
		auto key = __blend_get_kmer(fuzzy_seeds[index]);
		auto it = distinct_cores.find(key);

		if (it != distinct_cores.end()) {
			it->second++;
		} else {
			distinct_cores.emplace(key, 1);
		}
    }


	sequence.clear();

	if (fuzzy_seeds_len) {
		info->seeds = fuzzy_seeds;
		info->seeds_len = fuzzy_seeds_len;
		info->is_unique = (uint8_t *)calloc(fuzzy_seeds_len, sizeof(uint8_t));
	} else {
		info->seeds = NULL;
		info->seeds_len = 0;
		info->is_unique = NULL;
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
			if (g_info.chroms->array[chrom]->seeds == NULL) continue;
			if (g_info.chroms->array[chrom]->is_unique == NULL) { fprintf(stderr, "should not happend\n"); continue; }

			// good for readability
			uint128_t *fuzzy_seeds = g_info.chroms->array[chrom]->seeds;
			uint64_t seed_len = g_info.chroms->array[chrom]->seeds_len;
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

	{
		uint64_t chr_size = g_info.chroms->size;
		for (uint64_t chrom = 0; chrom < chr_size; chrom++) {
			relaxed_uniqueness_paint(g_info.chroms->array[chrom], seed_statistics, bed_files, print_bed);
		}
	}

	for (uint64_t chrom = 0; chrom < g_info.chroms->size; chrom++) {
		if (g_info.chroms->array[chrom]) {
			if (g_info.chroms->array[chrom]->seeds) free(g_info.chroms->array[chrom]->seeds);
			if (g_info.chroms->array[chrom]->is_unique) free(g_info.chroms->array[chrom]->is_unique);
			free(g_info.chroms->array[chrom]);
			g_info.chroms->array[chrom] = NULL;
		}
	}
	chrs_free(g_info.chroms);

	std::setprecision(5);

	std::cout << "{" << std::endl;

	// Param Settings
	std::cout << '\t' << "\"k\": " << kmer_size  << "," << std::endl;
	std::cout << '\t' << "\"w\": " << window 	 << "," << std::endl;
	std::cout << '\t' << "\"b\": " << blend_bits << "," << std::endl;
	std::cout << '\t' << "\"n\": " << n_number 	 << "," << std::endl;

	// Total Cores
	std::cout << '\t' << "\"total_seeds\": \"" << format_int(core_counts) << "\"," << std::endl;

	// Contiguous Cores
	std::cout << '\t' << "\"contig_seeds\": \"" << format_int(contiguous_counts) << "\"," << std::endl;

	// Distinct Cores
    std::cout << '\t' << "\"distinct_seeds\": \"" << format_int(distinct_cores.size()) << "\"," << std::endl;
	std::cout << '\t' << "\"distinct_seeds_ratio\": " << format_double((double)distinct_cores.size() / core_counts * 100, 5) << "," << std::endl;

	// Unique Cores
    std::cout << '\t' << "\"unique_seeds\": \"" << format_int(count_once) << "\"," << std::endl;
	std::cout << '\t' << "\"unique_seeds_ratio\": " << format_double((double)count_once / core_counts, 5) << "," << std::endl;

	// Genome Stats
	std::cout << '\t' << "\"genome_len\": " << g_info.genome_len << "," << std::endl;
	std::cout << '\t' << "\"valid_regions\": " << g_info.valid_genome_len << "," << std::endl;
	std::cout << '\t' << "\"genome_gap\": " << g_info.genome_len - g_info.valid_genome_len << "," << std::endl;
	std::cout << '\t' << "\"genome_gap_ratio\": " << (double)(g_info.genome_len - g_info.valid_genome_len) / g_info.genome_len * 100 << "," << std::endl;

	// Relaxed Uniqueness Ratio
	std::cout << '\t' << "\"details\": [" << std::endl;
	for (int i = 0; i <= RELAXATION_SPAN; i++) {
		std::cout << '\t' << '\t' << "{" << std::endl;
		std::cout << '\t' << '\t' << '\t' << "\"unique_seed_ratio\": " << format_double(((double)seed_statistics[i].relaxed_unique_count) / ((double)core_counts), 5) << "," << std::endl;
		std::cout << '\t' << '\t' << '\t' << "\"span_ratio\": " << format_double(((double)seed_statistics[i].total_relaxed_span) / g_info.genome_len, 5) << "," << std::endl;
		std::cout << '\t' << '\t' << '\t' << "\"gap_ratio\": " << format_double(((double)g_info.valid_genome_len - seed_statistics[i].total_relaxed_span) / g_info.genome_len, 5) << "," << std::endl;
		std::cout << '\t' << '\t' << '\t' << "\"gap_number\": " << seed_statistics[i].total_relaxed_gap_count << "," << std::endl;
		double total_relaxed_gap_span = g_info.genome_len - seed_statistics[i].total_relaxed_span;
		std::cout << '\t' << '\t' << '\t' << "\"gap_mean\": " << format_double(total_relaxed_gap_span / seed_statistics[i].total_relaxed_gap_count, 5) << std::endl;
		if (i == RELAXATION_SPAN)
			std::cout << '\t' << '\t' << "}" << std::endl;
		else 
			std::cout << '\t' << '\t' << "}," << std::endl;
	}
	std::cout << '\t' << "]," << std::endl;
	
	// Execution Time
	std::cout << '\t' << "\"time\": " << format_double(((double)durations.count()) / 1000) << "," << std::endl;

	// Mean Core Length
	std::cout << '\t' << "\"avg_len\": " << format_double(mean(lengths)) << "," << std::endl;

	// Std Dev of Lengths
	std::cout << '\t' << "\"std_dev_len\": " << format_double(stdev(lengths)) << "," << std::endl;

	// Min/Max Core Length
	std::cout << '\t' << "\"min_max_len\": \"" << format_int(length_gaps.minimum) << "/" << format_int(length_gaps.maximum) << "\"," << std::endl;

	// Total Sizes
	std::cout << '\t' << "\"size\": " << format_double(sizes / (1024.0 * 1024.0 * 1024.0)) << std::endl;

	std::cout << "}," << std::endl;

	return 0;
};