#include <chrono>
#include <unordered_map>
#include <algorithm>
#include <string>


#ifdef BLEND
#include "../sketch/blend.h"

#define __sketch_get_kmer(kmer) __blend_get_kmer(kmer)
#define __sketch_get_length(kmer) __blend_get_length(kmer)
#define __sketch_get_reference_id(kmer) __blend_get_reference_id(kmer)
#define __sketch_get_index(kmer) __blend_get_index(kmer)
#define __sketch_get_strand(kmer) __blend_get_strand(kmer)

#elif defined(MINIMIZER)
#include "../sketch/minimizer.h"

#define __sketch_get_kmer(kmer) __minimizer_get_kmer(kmer)
#define __sketch_get_length(kmer) __minimizer_get_length(kmer)
#define __sketch_get_reference_id(kmer) __minimizer_get_reference_id(kmer)
#define __sketch_get_index(kmer) __minimizer_get_index(kmer)
#define __sketch_get_strand(kmer) __minimizer_get_strand(kmer)

#elif defined(SYNCMER)
#include "../sketch/syncmer.h"

#define __sketch_get_kmer(kmer) __syncmer_get_kmer(kmer)
#define __sketch_get_length(kmer) __syncmer_get_length(kmer)
#define __sketch_get_reference_id(kmer) __syncmer_get_reference_id(kmer)
#define __sketch_get_index(kmer) __syncmer_get_index(kmer)
#define __sketch_get_strand(kmer) __syncmer_get_strand(kmer)

#elif defined(STROBEMER)
#include <tuple>
#include <vector>
#include "../sketch/strobemer/strobemer.hpp"

#define __sketch_get_kmer(kmer) ((kmer).x)
#define __sketch_get_length(kmer) (((kmer).y & 0xFFFFFFFF) + 22 - ((kmer).y >> 32))
#define __sketch_get_reference_id(kmer) 0
#define __sketch_get_index(kmer) ((kmer).y >> 32)
#define __sketch_get_strand(kmer) 0

using strobes_vector = std::vector<std::tuple<uint64_t, unsigned int, unsigned int, unsigned int, unsigned int>>;

#elif defined(STROBEMER2)
#include <vector>
#include <cstring>
#include "../sketch/strobemer/strobemer2.hpp"

#define __sketch_get_kmer(kmer) ((kmer).x)
#define __sketch_get_length(kmer) ((kmer).y & 0xFFFFFFFF)
#define __sketch_get_reference_id(kmer) 0
#define __sketch_get_index(kmer) ((kmer).y >> 32)
#define __sketch_get_strand(kmer) 0

using namespace strobemer2;

#elif defined(LCP)
#include <vector>
#include <cstring>
#include "../sketch/lps.h"

#define __sketch_get_kmer(kmer) ((kmer).x)
#define __sketch_get_length(kmer) ((kmer).y & 0xFFFFFFFF)
#define __sketch_get_reference_id(kmer) 0
#define __sketch_get_index(kmer) ((kmer).y >> 32)
#define __sketch_get_strand(kmer) 0


#else 
// dummy macros
#define __sketch_get_kmer(kmer) 0
#define __sketch_get_length(kmer) 0
#define __sketch_get_reference_id(kmer) 0
#define __sketch_get_index(kmer) 0
#define __sketch_get_strand(kmer) 0

#endif

// this to be included after sketch macro definitions
#include "utils.hpp"


void process(chr_info_t *info, std::string &sequence,
#ifdef BLEND 
             int window, int kmer_size, int blend_bits, int n_number,
#elif defined(MINIMIZER)
             int window, int kmer_size,
#elif defined(SYNCMER)
			 int smer_size, uint32_t kmer_size,
#elif defined(STROBEMER)
			 int num_strobe, int kmer_size, uint32_t w_min, uint32_t w_max,
#elif defined(STROBEMER2)
			 int read_len, int aux_len,
#elif defined(LCP)
			 int lcp_level, int dct_count,
#endif
			 uint64_t &core_counts,
             uint64_t &contiguous_counts,
             std::unordered_map<uint32_t, size_t> &distinct_cores,
             std::chrono::milliseconds &durations,
             uint64_t (&lengths)[DISTANCE_LENGTH],
			 gap &length_gaps,
             double &sizes,
			 uint64_t chr_idx) {
	
	auto start = std::chrono::high_resolution_clock::now();

	uint128_t *seeds;
    uint64_t seeds_len = 0;
#ifdef BLEND   
    seeds_len = sketch_blend(sequence.c_str(), sequence.length(), window, kmer_size, blend_bits, n_number, chr_idx, &seeds);
#elif defined(MINIMIZER)
    seeds_len = sketch_minimizers(sequence.c_str(), sequence.length(), window, kmer_size, 1, chr_idx, &seeds);
#elif defined(SYNCMER)
	seeds_len = sketch_syncmers(sequence.c_str(), sequence.length(), smer_size, kmer_size, chr_idx, &seeds);	
#elif defined(STROBEMER)
	strobes_vector minstrobes = seq_to_minstrobes2(num_strobe, kmer_size, w_min, w_max, sequence, chr_idx);
	seeds_len = minstrobes.size();
	seeds = (uint128_t *)malloc(sizeof(uint128_t) * seeds_len);
	uint64_t index = 0;
	for (const auto& t : minstrobes) {
		uint64_t strobemer_hash = std::get<0>(t);
        uint64_t strobe1_pos = std::get<2>(t);
        uint64_t strobe2_pos = std::get<3>(t);
		
		seeds[index].x = strobemer_hash;
		seeds[index++].y = (strobe1_pos << 32) | strobe2_pos;
	}
#elif defined(STROBEMER2)
	(void)chr_idx;
    int k = IndexParameters::DEFAULT, s = IndexParameters::DEFAULT;
    int l = IndexParameters::DEFAULT, u = IndexParameters::DEFAULT;
    int c = IndexParameters::DEFAULT, m = IndexParameters::DEFAULT;

	IndexParameters params = IndexParameters::from_read_length(read_len, k, s, l, u, c, m, aux_len);
	const int kk = params.syncmer.k;

	uint64_t seeds_cap = 1000000000;
	seeds = (uint128_t *)malloc(sizeof(uint128_t) * seeds_cap);

	RandstrobeGenerator gen(sequence, params.syncmer, params.randstrobe, RandstrobeHashMode::Classic);
	Randstrobe rs;
	while ((rs = gen.next()) != gen.end()) {
		if (seeds_len == seeds_cap) {
			seeds_cap *= 2;
			uint128_t *tmp = (uint128_t *)realloc(seeds, sizeof(uint128_t) * seeds_cap);
			if (!tmp) {
				std::cout << "Realloc failed" << std::endl;
				exit(EXIT_FAILURE);
			}
			seeds = tmp;
		}
		seeds[seeds_len].x = rs.hash;
		seeds[seeds_len++].y = ((uint64_t)rs.strobe1_pos << 32) | (rs.strobe2_pos + kk - rs.strobe1_pos);
	}
#elif defined(LCP)
	(void)chr_idx;
	struct lps lps_str;
	init_lps_segmented(&lps_str, sequence.c_str(), sequence.size(), lcp_level, dct_count);
	seeds = (uint128_t *)malloc(sizeof(uint128_t) * lps_str.size);
	for (int i = 0; i < lps_str.size; i++) {
		seeds[seeds_len].x = lps_str.cores[i].label;
		seeds[seeds_len++].y = ((uint64_t)lps_str.cores[i].start << 32) | ((uint64_t)lps_str.cores[i].end - (uint64_t)lps_str.cores[i].start);
	}
#else 
	printf("No method macro defined\n");
#endif
	auto extraction_end = std::chrono::high_resolution_clock::now();
	durations += std::chrono::milliseconds(std::chrono::duration_cast<std::chrono::milliseconds>(extraction_end - start).count());
	core_counts += seeds_len;
    sizes += sizeof(uint128_t) * seeds_len;
	
	analyze(seeds, seeds_len, contiguous_counts, lengths, length_gaps);

    for (uint64_t index = 0; index < seeds_len; index++) {
		auto key = __sketch_get_kmer(seeds[index]);
		auto it = distinct_cores.find(key);

		if (it != distinct_cores.end()) {
			it->second++;
		} else {
			distinct_cores.emplace(key, 1);
		}
    }

	sequence.clear();

	if (seeds_len) {
		info->seeds = seeds;
		info->seeds_len = seeds_len;
		info->is_unique = (uint8_t *)calloc(seeds_len, sizeof(uint8_t));
	} else {
		info->seeds = NULL;
		info->seeds_len = 0;
		info->is_unique = NULL;
	}
}

int main(int argc, char **argv) {

#ifdef BLEND 
	if (argc < 6) {
		std::cerr << "Wrong format: " << argv[0] << " [infile] [kmer-size] [window] [blend-bits] [n-number]" << std::endl;
		return -1;
	}
#elif defined(MINIMIZER)
	if (argc < 4) {
		std::cerr << "Wrong format: " << argv[0] << " [infile] [kmer-size] [window]" << std::endl;
		return -1;
	}
#elif defined(SYNCMER)
	if (argc < 4) {
		std::cerr << "Wrong format: " << argv[0] << " [infile] [kmer-size] [smer-size]" << std::endl;
		return -1;
	}
#elif defined(STROBEMER)
	if (argc < 6) {
		std::cerr << "Wrong format: " << argv[0] << " [infile] [num-strobe] [kmer-size] [w-min] [w-max]" << std::endl;
		return -1;
	}
#elif defined(STROBEME2)
	if (argc < 4) {
		std::cerr << "Wrong format: " << argv[0] << " [infile] [read-len] [aux-len]" << std::endl;
		return -1;
	}
#elif defined(LCP)
	if (argc < 4) {
		std::cerr << "Wrong format: " << argv[0] << " [infile] [lcp-level] [dct-count]" << std::endl;
		return -1;
	}
#endif

	std::ifstream input(argv[1]);
	if (!input.good()) {
		std::cerr << "Error opening: " << argv[1] << " . You have failed." << std::endl;
		return -1;
	}

	int print_bed = 0;

#ifdef BLEND 
    int kmer_size = atoi(argv[2]);
    int window = atoi(argv[3]);
	int blend_bits = atoi(argv[4]);
	int n_number = atoi(argv[5]);

	if (argc == 7) {
		print_bed = 1;
	}
#elif defined(MINIMIZER)
    int kmer_size = atoi(argv[2]);
    int window = atoi(argv[3]);

	if (argc == 5) {
		print_bed = 1;
	}
#elif defined(SYNCMER)
    int kmer_size = atoi(argv[2]);
    int smer_size = atoi(argv[3]);

	if (argc == 5) {
		print_bed = 1;
	}
#elif defined(STROBEMER)
	int num_strobe = atoi(argv[2]);
    int kmer_size = atoi(argv[3]);
    int w_min = atoi(argv[4]);
    int w_max = atoi(argv[5]);

	if (argc == 7) {
		print_bed = 1;
	}
#elif defined(STROBEMER2)
	int read_len = atoi(argv[2]);
	int aux_len = atoi(argv[3]);

	if (argc == 5) {
		print_bed = 1;
	}
#elif defined(LCP)
	int lcp_level = atoi(argv[2]);
	int dct_count = atoi(argv[3]);

	if (argc == 5) {
		print_bed = 1;
	}
#endif

	// variables
	std::string line;

	std::fstream genome;
	genome.open(argv[1], std::ios::in);

	genome_info_t g_info;
	g_info.chroms = chrs_init(128);

    // section 1
    uint64_t core_counts = 0;
	uint64_t contiguous_counts = 0;
    std::unordered_map<uint32_t, size_t> distinct_cores;
    std::chrono::milliseconds durations;
	// section 2
	uint64_t lengths[DISTANCE_LENGTH] = {0};
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

					chr_info_t *info = new chr_info_t;
					info->name = id;
					info->len = len;
					
					g_info.valid_genome_len += len - std::count(sequence.begin(), sequence.end(), 'N');
					g_info.genome_len += len;
					
					process(info, sequence, 
#ifdef BLEND 
							window, kmer_size, blend_bits, n_number, 
#elif defined(MINIMIZER)
             				window, kmer_size,
#elif defined(SYNCMER)
			  				smer_size, kmer_size,
#elif defined(STROBEMER)
							num_strobe, kmer_size, w_min, w_max,
#elif defined(STROBEMER2)
							read_len, aux_len,
#elif defined(LCP)
							lcp_level, dct_count,
#endif
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

			chr_info_t *info = new chr_info_t;
			info->name = id;
			info->len = len;

			g_info.valid_genome_len += len - std::count(sequence.begin(), sequence.end(), 'N');
			g_info.genome_len += len;

			process(info, sequence, 
#ifdef BLEND 
					window, kmer_size, blend_bits, n_number, 
#elif defined(MINIMIZER)
             		window, kmer_size,
#elif defined(SYNCMER)
			  		smer_size, kmer_size,
#elif defined(STROBEMER)
					num_strobe, kmer_size, w_min, w_max,
#elif defined(STROBEMER2)
					read_len, aux_len,
#elif defined(LCP)
					lcp_level, dct_count,
#endif
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
			uint128_t *seeds = g_info.chroms->array[chrom]->seeds;
			uint64_t seed_len = g_info.chroms->array[chrom]->seeds_len;
			uint8_t *is_unique = g_info.chroms->array[chrom]->is_unique;
			
			for (uint64_t seed_index = 0; seed_index < seed_len; seed_index++) {
				auto it = distinct_cores.find(__sketch_get_kmer(seeds[seed_index]));
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
#ifdef BLEND 
	std::cout << '\t' << "\"k\": " << kmer_size  << "," << std::endl;
	std::cout << '\t' << "\"w\": " << window 	 << "," << std::endl;
	std::cout << '\t' << "\"b\": " << blend_bits << "," << std::endl;
	std::cout << '\t' << "\"n\": " << n_number 	 << "," << std::endl;
#elif defined(MINIMIZER)
	std::cout << '\t' << "\"k\": " << kmer_size  << "," << std::endl;
	std::cout << '\t' << "\"w\": " << window 	 << "," << std::endl;
#elif defined(SYNCMER)
	std::cout << '\t' << "\"k\": " << kmer_size  << "," << std::endl;
	std::cout << '\t' << "\"s\": " << smer_size  << "," << std::endl;
#elif defined(STROBEMER)
	std::cout << '\t' << "\"n\": " 	   << num_strobe << "," << std::endl;
	std::cout << '\t' << "\"k\": " 	   << kmer_size  << "," << std::endl;
	std::cout << '\t' << "\"w_min\": " << w_min  	 << "," << std::endl;
	std::cout << '\t' << "\"w_max\": " << w_max  	 << "," << std::endl;
#elif defined(STROBEME2)
	std::cout << '\t' << "\"r\": " << read_len << "," << std::endl;
	std::cout << '\t' << "\"a\": " << aux_len  << "," << std::endl;
#elif defined(LCP)
	std::cout << '\t' << "\"l\": " << lcp_level << "," << std::endl;
	std::cout << '\t' << "\"d\": " << dct_count  << "," << std::endl;
#endif

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
		std::cout << '\t' << '\t' << '\t' << "\"gap_mean\": " << format_double_classic(total_relaxed_gap_span / seed_statistics[i].total_relaxed_gap_count, 5) << std::endl;
		if (i == RELAXATION_SPAN)
			std::cout << '\t' << '\t' << "}" << std::endl;
		else 
			std::cout << '\t' << '\t' << "}," << std::endl;
	}
	std::cout << '\t' << "]," << std::endl;
	
	// Execution Time
	std::cout << '\t' << "\"time\": " << format_double_classic(((double)durations.count()) / 1000) << "," << std::endl;

	// Mean Core Length
	std::cout << '\t' << "\"avg_len\": " << format_double_classic(mean(lengths)) << "," << std::endl;

	// Std Dev of Lengths
	std::cout << '\t' << "\"std_dev_len\": " << format_double_classic(stdev(lengths)) << "," << std::endl;

	// Min/Max Core Length
	std::cout << '\t' << "\"min_max_len\": \"" << format_int(length_gaps.minimum) << "/" << format_int(length_gaps.maximum) << "\"," << std::endl;

	// Total Sizes
	std::cout << '\t' << "\"size\": " << format_double(sizes / (1024.0 * 1024.0 * 1024.0)) << std::endl;

	std::cout << "}," << std::endl;

	return 0;
};