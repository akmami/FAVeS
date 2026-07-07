//
//  strobemer2_main.cpp
//  
//  The implementation refined from https://github.com/ksahlin/strobealign
//
//  For a randstrobe start = strobe1_pos, end = strobe2_pos + k
//  For a syncmer start = position, end = position + k.
//
//  Usage:
//      strobemer2_main [options] <ref.fasta>
//        --syncmers            dump syncmers instead of randstrobes
//        --mode classic|strand randstrobe hash scheme (default: classic,
//                              which matches strobealign's C++ dumpstrobes)
//        -r INT                read-length preset (default 150)
//        -k INT -s INT         override k / s
//        -l INT -u INT         override w_min / w_max
//        -c INT                override syncmer bitcount (q = 2^c - 1)
//        -m INT                override max seed length
//        -o FILE               write to FILE instead of stdout
//
//  Build:
//      g++ -std=c++17 -O2 strobemer2.cpp strobemer2_main.cpp -o strobemer2_main
//

#include <cstdint>
#include <cstring>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "strobemer2.hpp"

using namespace strobemer2;

// ---- FASTA reading (replicates strobealign cpp/refs.cpp semantics) ----------
// - concatenates sequence lines per record
// - cuts the contig name at the first whitespace
// - uppercases the sequence with `c & ~32` (as strobealign does)
struct Contig {
    std::string name;
    std::string seq;
};

static void to_uppercase(std::string& s) {
    for (char& c : s) {
        c = static_cast<char>(static_cast<unsigned char>(c) & ~32);
    }
}

static std::vector<Contig> read_fasta(std::istream& in) {
    std::vector<Contig> contigs;
    std::string line, name, seq;
    auto flush = [&]() {
        if (!seq.empty()) {
            to_uppercase(seq);
            contigs.push_back(Contig{name, seq});
            seq.clear();
        }
    };
    bool first = true;
    while (std::getline(in, line)) {
        if (!line.empty() && line.back() == '\r') {
            line.pop_back();  // tolerate CRLF
        }
        if (!line.empty() && line[0] == '>') {
            flush();
            std::string::size_type space = line.find_first_of(" \t\f\v\n\r");
            name = space < line.length() ? line.substr(1, space - 1) : line.substr(1);
            first = false;
        } else if (!first) {
            seq += line;
        }
    }
    flush();
    return contigs;
}

static int parse_int(const char* a) { return static_cast<int>(std::strtol(a, nullptr, 10)); }

int main(int argc, char** argv) {
    bool syncmers = false;
    RandstrobeHashMode mode = RandstrobeHashMode::Classic;
    int r = 150;
    int k = IndexParameters::DEFAULT, s = IndexParameters::DEFAULT;
    int l = IndexParameters::DEFAULT, u = IndexParameters::DEFAULT;
    int c = IndexParameters::DEFAULT, m = IndexParameters::DEFAULT;
    int aux_len = 17;
    std::string ref_path, out_path;

    for (int i = 1; i < argc; i++) {
        std::string arg = argv[i];
        auto need = [&](const char* what) -> const char* {
            if (i + 1 >= argc) {
                std::cerr << "Missing value for " << what << "\n";
                std::exit(2);
            }
            return argv[++i];
        };
        if (arg == "--syncmers") syncmers = true;
        else if (arg == "--mode") {
            std::string mv = need("--mode");
            if (mv == "classic") mode = RandstrobeHashMode::Classic;
            else if (mv == "strand" || mv == "strand-aware") mode = RandstrobeHashMode::StrandAware;
            else { std::cerr << "Unknown mode: " << mv << "\n"; return 2; }
        }
        else if (arg == "-r") r = parse_int(need("-r"));
        else if (arg == "-k") k = parse_int(need("-k"));
        else if (arg == "-s") s = parse_int(need("-s"));
        else if (arg == "-l") l = parse_int(need("-l"));
        else if (arg == "-u") u = parse_int(need("-u"));
        else if (arg == "-c") c = parse_int(need("-c"));
        else if (arg == "-m") m = parse_int(need("-m"));
        else if (arg == "-o") out_path = need("-o");
        else if (!arg.empty() && arg[0] == '-') {
            std::cerr << "Unknown option: " << arg << "\n";
            return 2;
        }
        else ref_path = arg;
    }

    if (ref_path.empty()) {
        std::cerr << "Usage: strobemer2_dump [--syncmers] [--mode classic|strand] "
                     "[-r N] [-k N -s N -l N -u N -c N -m N] [-o out] <ref.fasta>\n";
        return 2;
    }

    IndexParameters params = IndexParameters::from_read_length(r, k, s, l, u, c, m, aux_len);

    std::ifstream in(ref_path);
    if (!in) {
        std::cerr << "Cannot open " << ref_path << "\n";
        return 1;
    }
    auto contigs = read_fasta(in);

    std::ofstream fout;
    std::ostream* os = &std::cout;
    if (!out_path.empty()) {
        fout.open(out_path);
        if (!fout) { std::cerr << "Cannot write " << out_path << "\n"; return 1; }
        os = &fout;
    }

    const int kk = params.syncmer.k;
    for (const auto& contig : contigs) {
        if (syncmers) {
            SyncmerIterator it(contig.seq, params.syncmer);
            Syncmer syncmer;
            while (!(syncmer = it.next()).is_end()) {
                *os << contig.name << '\t' << syncmer.position << '\t'
                    << (syncmer.position + kk) << '\t' << syncmer.hash << '\n';
            }
        } else {
            // Use RandstrobeGenerator to mirror dumpstrobes' default path exactly.
            RandstrobeGenerator gen(contig.seq, params.syncmer, params.randstrobe, mode);
            Randstrobe rs;
            while ((rs = gen.next()) != gen.end()) {
                *os << contig.name << '\t' << rs.strobe1_pos << '\t'
                    << (rs.strobe2_pos + kk) << '\t' << rs.hash << '\n';
            }
        }
    }
    return 0;
}
