//
//  strobemer2.hpp
//
//  Self-contained extraction of the strobemer (randstrobe) sketching core
//  from the strobealign project (https://github.com/ksahlin/strobealign).
//
//  This file is *not* the ksahlin/strobemers reference implementation (that is
//  strobemer.hpp/.cpp in this directory). It is the seeding core that strobealign
//  actually uses: canonical open-syncmers + randstrobes ("method 3'").
//
//  Provenance / consistency note
//  -----------------------------
//  strobealign currently ships two copies of the sketcher:
//    * cpp/randstrobes.cpp  (C++, now frozen/legacy)
//    * src/seeding/*.rs      (Rust, actively developed)
//
//  Their SYNCMER generation and their STROBE-PAIRING (window + method-3' popcount
//  minimizer) are byte-for-byte equivalent. They differ in exactly ONE place:
//  how the two chosen syncmer hashes are combined into the final randstrobe hash.
//
//    * Classic (C++ cpp/randstrobes.cpp):
//        hash = ((h1 & main_mask) | (h2 & ~main_mask)) & RANDSTROBE_HASH_MASK
//      No strand information is stored in the hash.
//
//    * StrandAware (Rust src/seeding/strobes.rs):
//        The layout becomes
//          <strobe1 hash><strobe1 orientation><strobe2 hash><strobe2 orientation><offset>
//        i.e. one orientation bit for each strobe is folded into the hash.
//
//  Because that is the only divergence, this file implements ONE shared core and
//  exposes both hash schemes through a `RandstrobeHashMode` selector, so you can
//  reproduce either strobealign build exactly.
//
//  API (see bottom of file):
//    * randstrobes(seq, params, mode)        -> forward-strand randstrobes
//    * randstrobes_query(seq, params, mode)  -> {forward, reverse-complement}
//    * canonical_syncmers(seq, syncmer_params)
//    * make_index_parameters(read_length, ...) / IndexParameters::from_read_length
//

#ifndef __STROBEMER2_HPP__
#define __STROBEMER2_HPP__

#include <array>
#include <cstdint>
#include <deque>
#include <iostream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

namespace strobemer2 {

using syncmer_hash_t = uint64_t;
using randstrobe_hash_t = uint64_t;

// Top 56 bits of the hash; the low 8 bits are reserved (strobe2 offset region).
static constexpr uint64_t RANDSTROBE_HASH_MASK = 0xFFFFFFFFFFFFFF00ULL;
// Number of low bits reserved for the strobe2 offset (matches Rust STROBE2_OFFSET_BITS).
static constexpr uint32_t STROBE2_OFFSET_BITS = 8;

// Selects how the two syncmer hashes are combined into a randstrobe hash.
enum class RandstrobeHashMode {
    Classic,      // strobealign C++ (cpp/randstrobes.cpp) - no orientation bits
    StrandAware,  // strobealign Rust (src/seeding/strobes.rs) - orientation bits folded in
};

// --------------------------------------------------------------------------
// Parameters
// --------------------------------------------------------------------------

struct SyncmerParameters {
    int k;
    int s;
    int t_syncmer;  // = (k - s) / 2 + 1

    SyncmerParameters(int k, int s)
        : k(k), s(s), t_syncmer((k - s) / 2 + 1) {
        verify();
    }

    void verify() const {
        if (k <= 7 || k > 32) {
            throw std::invalid_argument("k not in [8,32]");
        }
        if (s > k) {
            throw std::invalid_argument("s is larger than k");
        }
        if ((k - s) % 2 != 0) {
            throw std::invalid_argument(
                "(k - s) must be even to create canonical syncmers "
                "(set s to k-2, k-4, k-6, ...)");
        }
    }

    bool operator==(const SyncmerParameters& o) const {
        return k == o.k && s == o.s && t_syncmer == o.t_syncmer;
    }
};

struct RandstrobeParameters {
    uint64_t q;
    int max_dist;
    unsigned w_min;
    unsigned w_max;

    // Mask for the strobe1 ("main") hash bits.
    uint64_t main_hash_mask;
    // StrandAware only: main hash mask including the strobe1 orientation bit.
    uint64_t forward_main_hash_mask;
    // StrandAware only: bit position of the strobe1 orientation bit.
    uint32_t partial_orientation_pos;

    RandstrobeParameters(uint64_t q, int max_dist, unsigned w_min, unsigned w_max,
                         uint64_t main_hash_mask)
        : q(q), max_dist(max_dist), w_min(w_min), w_max(w_max),
          main_hash_mask(main_hash_mask) {
        // trailing_zeros(main_hash_mask) - 1  (main_hash_mask = ~0 << (9 + aux_len))
        uint32_t tz = 0;
        while (tz < 64 && ((main_hash_mask >> tz) & 1ULL) == 0) {
            ++tz;
        }
        partial_orientation_pos = tz - 1;
        forward_main_hash_mask =
            main_hash_mask | (1ULL << partial_orientation_pos);
        if (w_min > w_max) {
            throw std::invalid_argument("w_min is greater than w_max");
        }
    }

    bool operator==(const RandstrobeParameters& o) const {
        return q == o.q && max_dist == o.max_dist && w_min == o.w_min &&
               w_max == o.w_max && main_hash_mask == o.main_hash_mask;
    }
};

// Bundles syncmer + randstrobe parameters, mirroring strobealign's IndexParameters.
struct IndexParameters {
    SyncmerParameters syncmer;
    RandstrobeParameters randstrobe;

    IndexParameters(SyncmerParameters syncmer, RandstrobeParameters randstrobe)
        : syncmer(syncmer), randstrobe(randstrobe) {}

    static constexpr int DEFAULT = INT32_MIN;

    // Reproduces strobealign's IndexParameters::from_read_length default profiles.
    // Any of k, s, w_min, w_max, c, max_seed_len, aux_len may be overridden;
    // leave them as DEFAULT to use the profile value.
    static IndexParameters from_read_length(
        int read_length, int k = DEFAULT, int s = DEFAULT, int w_min = DEFAULT,
        int w_max = DEFAULT, int c = DEFAULT, int max_seed_len = DEFAULT,
        int aux_len = DEFAULT);
};

// --------------------------------------------------------------------------
// Output structures (as in strobealign)
// --------------------------------------------------------------------------

struct Syncmer {
    syncmer_hash_t hash;
    size_t position;   // 0-based start position of the k-mer in the sequence
    bool is_forward;   // canonical k-mer taken from the forward strand?

    bool is_end() const { return hash == 0 && position == 0; }
    bool operator==(const Syncmer& o) const {
        return hash == o.hash && position == o.position;
    }
};

struct Randstrobe {
    randstrobe_hash_t hash;          // hash of (strobe1, strobe2)
    randstrobe_hash_t hash_revcomp;  // hash of (strobe2, strobe1) for the rev-comp strand
    unsigned int strobe1_pos;        // start of first strobe (k-mer)
    unsigned int strobe2_pos;        // start of second strobe (k-mer)

    bool operator==(const Randstrobe& o) const {
        return hash == o.hash && strobe1_pos == o.strobe1_pos &&
               strobe2_pos == o.strobe2_pos;
    }
    bool operator!=(const Randstrobe& o) const { return !(*this == o); }
};

// Query-side randstrobe: like Randstrobe but with [start, end) span coordinates
// (end = strobe2_pos + k), matching strobealign's QueryRandstrobe.
struct QueryRandstrobe {
    randstrobe_hash_t hash;
    randstrobe_hash_t hash_revcomp;
    unsigned int start;
    unsigned int end;
};

std::ostream& operator<<(std::ostream& os, const Syncmer& syncmer);
std::ostream& operator<<(std::ostream& os, const Randstrobe& randstrobe);
std::ostream& operator<<(std::ostream& os, const QueryRandstrobe& randstrobe);

// --------------------------------------------------------------------------
// Iterators
// --------------------------------------------------------------------------

// Streams canonical open-syncmers out of a sequence. `next()` returns the
// end marker (is_end() == true) once the sequence is exhausted.
class SyncmerIterator {
public:
    SyncmerIterator(std::string_view seq, SyncmerParameters parameters)
        : seq(seq), parameters(parameters) {}

    Syncmer next();

private:
    const std::string_view seq;
    const SyncmerParameters parameters;

    const uint64_t kmask = (1ULL << (2 * parameters.k)) - 1;
    const uint64_t smask = (1ULL << (2 * parameters.s)) - 1;
    const uint64_t kshift = (parameters.k - 1) * 2;
    const uint64_t sshift = (parameters.s - 1) * 2;
    std::deque<uint64_t> qs;  // s-mer hashes
    uint64_t qs_min_val = UINT64_MAX;
    size_t l = 0;
    uint64_t xk[2] = {0, 0};
    uint64_t xs[2] = {0, 0};
    size_t i = 0;
};

std::vector<Syncmer> canonical_syncmers(std::string_view seq,
                                        SyncmerParameters parameters);

// Iterates over randstrobes given a precomputed vector of syncmers.
class RandstrobeIterator {
public:
    RandstrobeIterator(const std::vector<Syncmer>& syncmers,
                       RandstrobeParameters parameters,
                       RandstrobeHashMode mode)
        : syncmers(syncmers), parameters(parameters), mode(mode) {
        if (parameters.w_min > parameters.w_max) {
            throw std::invalid_argument("w_min is greater than w_max");
        }
    }

    Randstrobe next() { return get(strobe1_index++); }
    bool has_next() const { return strobe1_index < syncmers.size(); }

private:
    Randstrobe get(unsigned int strobe1_index) const;
    const std::vector<Syncmer>& syncmers;
    const RandstrobeParameters parameters;
    const RandstrobeHashMode mode;
    unsigned strobe1_index = 0;
};

// Generates randstrobes while producing syncmers on the fly (low memory).
class RandstrobeGenerator {
public:
    RandstrobeGenerator(std::string_view seq,
                        SyncmerParameters syncmer_parameters,
                        RandstrobeParameters randstrobe_parameters,
                        RandstrobeHashMode mode)
        : syncmer_iterator(seq, syncmer_parameters),
          parameters(randstrobe_parameters),
          mode(mode) {}

    Randstrobe next();
    Randstrobe end() const { return Randstrobe{0, 0, 0, 0}; }

private:
    SyncmerIterator syncmer_iterator;
    const RandstrobeParameters parameters;
    const RandstrobeHashMode mode;
    std::deque<Syncmer> syncmers;
};

// --------------------------------------------------------------------------
// Convenience API  --  "give a string, get strobemers"
// --------------------------------------------------------------------------

// Forward-strand randstrobes for `seq`.
std::vector<Randstrobe> randstrobes(std::string_view seq,
                                    const IndexParameters& parameters,
                                    RandstrobeHashMode mode);

// Forward-strand randstrobes using a read-length preset (single call helper).
std::vector<Randstrobe> randstrobes(
    std::string_view seq, int read_length,
    RandstrobeHashMode mode = RandstrobeHashMode::StrandAware);

// Randstrobes for the forward sequence ([0]) and its reverse complement ([1]),
// with [start, end) coordinates. Mirrors strobealign::randstrobes_query.
std::array<std::vector<QueryRandstrobe>, 2> randstrobes_query(
    std::string_view seq, const IndexParameters& parameters,
    RandstrobeHashMode mode);

}  // namespace strobemer2

#endif  // STROBEMER2_HPP
