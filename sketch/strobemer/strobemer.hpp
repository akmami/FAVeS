//
//  strobemer.hpp
//  cpp_course
//
//  Created by Kristoffer Sahlin on 4/21/21.
//

#ifndef strobemer_hpp
#define strobemer_hpp

#include <stdio.h>
#include <string>
#include <vector>
#include <deque>
#include <tuple>
#include <iostream>
#include <math.h>       /* pow */
#include <bitset>
#include "robin_hood.h"


typedef std::vector< std::tuple<uint64_t, unsigned int, unsigned int, unsigned int, unsigned int>> mers_vector;

mers_vector seq_to_kmers(int k, std::string &seq, unsigned int ref_index);
mers_vector seq_to_randstrobes2(int n, int k, int w_min, int w_max, std::string &seq, unsigned int ref_index);
mers_vector seq_to_minstrobes2(int n, int k, int w_min, int w_max, std::string &seq, unsigned int ref_index);
mers_vector seq_to_hybridstrobes2(int n, int k, int w_min, int w_max, std::string &seq, unsigned int ref_index);
mers_vector seq_to_randstrobes3(int n, int k, int w_min, int w_max, std::string &seq, unsigned int ref_index);
mers_vector seq_to_hybridstrobes3(int n, int k, int w_min, int w_max, std::string &seq, unsigned int ref_index);

typedef robin_hood::unordered_map< unsigned int, std::vector< std::tuple<uint64_t, unsigned int, unsigned int, unsigned int, unsigned int>>> three_pos_index;
mers_vector construct_flat_vector_three_pos(three_pos_index &tmp_index, uint64_t &unique_elements);
typedef robin_hood::unordered_map< uint64_t, std::tuple<uint64_t, unsigned int >> kmer_lookup;
void index_vector_three_pos(mers_vector  &mers_vector, kmer_lookup &mers_index);


struct hit {
    unsigned int query_s;
    unsigned int query_e;
    unsigned int ref_s;
    unsigned int ref_e;
};

struct nam {
    unsigned int ref_id;
    unsigned int query_s;
    unsigned int query_e;
    unsigned int ref_s;
    unsigned int ref_e;
    unsigned int n_hits = 0;
    unsigned int previous_query_start;
    unsigned int previous_ref_start;
};
#endif /* strobemer_hpp */