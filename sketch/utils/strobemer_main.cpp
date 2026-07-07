//
//  strobemer_main.cpp
//  
//  The implementation refined from https://github.com/ksahlin/strobemers 
//

#include <iostream>
#include <string>
#include <vector>
#include <tuple>

#include "strobemer.hpp"

using strobes_vector = std::vector<std::tuple<uint64_t, unsigned int, unsigned int, unsigned int, unsigned int>>;

int main(int argc, char* argv[]) {
    if (argc < 2) {
        std::cerr << "Usage:\n";
        std::cerr << "  " << argv[0] << " rs2\n";
        std::cerr << "  " << argv[0] << " rs3\n";
        std::cerr << "  " << argv[0] << " ms2\n";
        return 1;
    }

    std::string method = argv[1];

    std::string seq = "ACGCGTACGAATCACGCCGGGTGTGTGTGATCGGGGCTATCAGCTACGTACTATGCTAGCTACGGACGGCGATTTTTTTTCATATCGTACGCTAGCTAGCTAGCTGCGATCGATTCG";

    unsigned int n = 3;
    unsigned int k = 15;
    unsigned int w_min = 16;
    unsigned int w_max = 30;
    unsigned int seq_id = 0;

    strobes_vector strobes;

    if (method == "rs2") {
        n = 2;
        strobes = seq_to_randstrobes2(n, k, w_min, w_max, seq, seq_id);
    } else if (method == "rs3") {
        n = 3;
        strobes = seq_to_randstrobes3(n, k, w_min, w_max, seq, seq_id);
    } else if (method == "ms2") {
        n = 2;
        strobes = seq_to_minstrobes2(n, k, w_min, w_max, seq, seq_id);
    } else {
        std::cerr << "Unknown method: " << method << "\n";
        return 1;
    }

    std::cout << "Method: " << method << "\n";
    std::cout << "Number of strobes: " << strobes.size() << "\n\n";

    for (const auto& t : strobes) {
        uint64_t strobemer_hash = std::get<0>(t);
        unsigned int strobe1_pos = std::get<2>(t);
        unsigned int strobe2_pos = std::get<3>(t);
        unsigned int strobe3_pos = std::get<4>(t);

        std::string strobemer;

        if (n == 2) {
            strobemer =
                seq.substr(strobe1_pos, k) +
                seq.substr(strobe2_pos, k);
        } else {
            strobemer =
                seq.substr(strobe1_pos, k) +
                seq.substr(strobe2_pos, k) +
                seq.substr(strobe3_pos, k);
        }

        std::cout << "Hash: " << strobemer_hash << "\n";
        std::cout << "Positions: ";

        if (n == 2) {
            std::cout << strobe1_pos << ", " << strobe2_pos << "\n";
        } else {
            std::cout << strobe1_pos << ", "
                      << strobe2_pos << ", "
                      << strobe3_pos << "\n";
        }

        std::cout << "Sequence: " << strobemer << "\n\n";
    }

    return 0;
}
