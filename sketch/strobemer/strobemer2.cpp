//
//  strobemer2.cpp
//
//  Self-contained extraction of strobealign's syncmer + randstrobe sketching.
//  See strobemer2.hpp for the provenance / consistency notes.
//
//  Ported verbatim from:
//    * cpp/randstrobes.cpp        (syncmer generation, strobe pairing, Classic hash)
//    * cpp/hash.hpp               (single-value xxh64)
//    * src/seeding/strobes.rs     (StrandAware hash + Syncmer.is_forward)
//    * cpp/indexparameters.cpp    (from_read_length default profiles)
//

#include "strobemer2.hpp"

#include <algorithm>
#include <bitset>
#include <cassert>
#include <cmath>
#include <limits>

namespace strobemer2 {

// --------------------------------------------------------------------------
// xxh64 (single 64-bit value) -- identical in cpp/hash.hpp and src/seeding/hash.rs
// --------------------------------------------------------------------------

#if !defined(NO_CLANG_BUILTIN) && defined(__has_builtin)
#  if __has_builtin(__builtin_rotateleft64)
#    define STROBEMER2_ROTL64 __builtin_rotateleft64
#  endif
#endif
#ifndef STROBEMER2_ROTL64
#  define STROBEMER2_ROTL64(x, r) (((x) << (r)) | ((x) >> (64 - (r))))
#endif

static constexpr uint64_t XXH_PRIME64_1 = 0x9E3779B185EBCA87ULL;
static constexpr uint64_t XXH_PRIME64_2 = 0xC2B2AE3D27D4EB4FULL;
static constexpr uint64_t XXH_PRIME64_3 = 0x165667B19E3779F9ULL;
static constexpr uint64_t XXH_PRIME64_4 = 0x85EBCA77C2B2AE63ULL;
static constexpr uint64_t XXH_PRIME64_5 = 0x27D4EB2F165667C5ULL;

static inline uint64_t xxh64(uint64_t input) {
    uint64_t result = XXH_PRIME64_5 + 8;
    input *= XXH_PRIME64_2;
    input = STROBEMER2_ROTL64(input, 31);
    result ^= input * XXH_PRIME64_1;
    result = STROBEMER2_ROTL64(result, 27);
    result = result * XXH_PRIME64_1 + XXH_PRIME64_4;
    result ^= result >> 33;
    result = result * XXH_PRIME64_2;
    result ^= result >> 29;
    result = result * XXH_PRIME64_3;
    result ^= result >> 32;
    return result;
}

static inline syncmer_hash_t syncmer_kmer_hash(uint64_t packed) {
    return xxh64(packed);
}
static inline syncmer_hash_t syncmer_smer_hash(uint64_t packed) {
    return xxh64(packed);
}

// a, A -> 0 ; c, C -> 1 ; g, G -> 2 ; t, T, u, U -> 3 ; anything else -> 4 (N)
static const unsigned char seq_nt4_table[256] = {
    0, 1, 2, 3,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  3, 3, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  3, 3, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4
};

// --------------------------------------------------------------------------
// randstrobe hash combinations
// --------------------------------------------------------------------------

// Classic strobealign (cpp/randstrobes.cpp): top bits from the main hash,
// bottom bits from the auxiliary hash, low 8 bits cleared. No orientation info.
static inline randstrobe_hash_t classic_randstrobe_hash(
    syncmer_hash_t hash1, syncmer_hash_t hash2, uint64_t main_hash_mask) {
    return ((hash1 & main_hash_mask) | (hash2 & ~main_hash_mask)) &
           RANDSTROBE_HASH_MASK;
}

// Strand-aware strobealign (src/seeding/strobes.rs). Layout:
//   <strobe1 hash><strobe1 orientation><strobe2 hash><strobe2 orientation><offset>
// The final mask (REF_RANDSTROBE_HASH_MASK << 1) clears the low 9 bits, so on the
// query side the strobe2 orientation bit and the offset region are zeroed -- this
// is faithful to the Rust Randstrobe::hash().
static inline randstrobe_hash_t strand_aware_randstrobe_hash(
    syncmer_hash_t hash1, syncmer_hash_t hash2, uint64_t is_forward1,
    uint64_t is_forward2, const RandstrobeParameters& p) {
    return ((hash1 & p.main_hash_mask) | (hash2 & ~p.forward_main_hash_mask) |
            (is_forward1 << p.partial_orientation_pos) |
            (is_forward2 << STROBE2_OFFSET_BITS)) &
           (RANDSTROBE_HASH_MASK << 1);
}

static inline Randstrobe make_randstrobe(const Syncmer& strobe1,
                                         const Syncmer& strobe2,
                                         const RandstrobeParameters& p,
                                         RandstrobeHashMode mode) {
    randstrobe_hash_t hash;
    randstrobe_hash_t hash_revcomp;
    if (mode == RandstrobeHashMode::Classic) {
        hash = classic_randstrobe_hash(strobe1.hash, strobe2.hash, p.main_hash_mask);
        hash_revcomp =
            classic_randstrobe_hash(strobe2.hash, strobe1.hash, p.main_hash_mask);
    } else {
        uint64_t f1 = strobe1.is_forward ? 1 : 0;
        uint64_t f2 = strobe2.is_forward ? 1 : 0;
        hash = strand_aware_randstrobe_hash(strobe1.hash, strobe2.hash, f1, f2, p);
        hash_revcomp = strand_aware_randstrobe_hash(strobe2.hash, strobe1.hash,
                                                    f1 ^ 1, f2 ^ 1, p);
    }
    return Randstrobe{hash, hash_revcomp,
                      static_cast<uint32_t>(strobe1.position),
                      static_cast<uint32_t>(strobe2.position)};
}

// --------------------------------------------------------------------------
// syncmer generation (verbatim from cpp/randstrobes.cpp, with is_forward from Rust)
// --------------------------------------------------------------------------

Syncmer SyncmerIterator::next() {
    for (; i < seq.length(); ++i) {
        int c = seq_nt4_table[(uint8_t)seq[i]];
        if (c < 4) {  // not an "N" base
            xk[0] = (xk[0] << 2 | c) & kmask;                  // forward strand
            xk[1] = xk[1] >> 2 | (uint64_t)(3 - c) << kshift;  // reverse strand
            xs[0] = (xs[0] << 2 | c) & smask;                  // forward strand
            xs[1] = xs[1] >> 2 | (uint64_t)(3 - c) << sshift;  // reverse strand
            if (++l < static_cast<size_t>(parameters.s)) {
                continue;
            }
            // we find an s-mer
            uint64_t ys = std::min(xs[0], xs[1]);
            uint64_t hash_s = syncmer_smer_hash(ys);
            qs.push_back(hash_s);
            // not enough hashes in the queue, yet
            if (qs.size() < static_cast<size_t>(parameters.k - parameters.s + 1)) {
                continue;
            }
            if (qs.size() == static_cast<size_t>(parameters.k - parameters.s + 1)) {
                // last s-mer within the first k-mer, decide if we add it
                for (size_t j = 0; j < qs.size(); j++) {
                    if (qs[j] <= qs_min_val) {
                        qs_min_val = qs[j];
                    }
                }
            } else {
                // update queue and current minimum and position
                uint64_t front = qs.front();
                qs.pop_front();
                if (front == qs_min_val) {
                    // we popped a minimum, find new brute force
                    qs_min_val = UINT64_MAX;
                    for (size_t j = 0; j < qs.size(); j++) {
                        if (qs[j] <= qs_min_val) {
                            qs_min_val = qs[j];
                        }
                    }
                } else if (hash_s < qs_min_val) {
                    // the new value added to queue is the new minimum
                    qs_min_val = hash_s;
                }
            }
            if (qs[parameters.t_syncmer - 1] == qs_min_val) {
                // occurs at t:th position in k-mer -> this is a syncmer
                uint64_t yk = std::min(xk[0], xk[1]);
                auto syncmer = Syncmer{syncmer_kmer_hash(yk),
                                       i - parameters.k + 1, xk[0] <= xk[1]};
                i++;
                return syncmer;
            }
        } else {
            // if there is an "N", restart
            qs_min_val = UINT64_MAX;
            l = xs[0] = xs[1] = xk[0] = xk[1] = 0;
            qs.clear();
        }
    }
    return Syncmer{0, 0, true};  // end marker
}

std::vector<Syncmer> canonical_syncmers(std::string_view seq,
                                        SyncmerParameters parameters) {
    std::vector<Syncmer> syncmers;
    SyncmerIterator syncmer_iterator{seq, parameters};
    Syncmer syncmer;
    while (!(syncmer = syncmer_iterator.next()).is_end()) {
        syncmers.push_back(syncmer);
    }
    return syncmers;
}

// --------------------------------------------------------------------------
// randstrobe pairing (method 3': popcount minimizer over the window)
// --------------------------------------------------------------------------

Randstrobe RandstrobeIterator::get(unsigned int strobe1_index) const {
    unsigned int w_end = std::min(
        static_cast<size_t>(strobe1_index + parameters.w_max), syncmers.size() - 1);

    auto strobe1 = syncmers[strobe1_index];
    auto max_position = strobe1.position + parameters.max_dist;
    unsigned int w_start = strobe1_index + parameters.w_min;
    uint64_t min_val = std::numeric_limits<uint64_t>::max();
    Syncmer strobe2 = strobe1;  // default if no nearby syncmer

    for (auto i = w_start; i <= w_end && syncmers[i].position <= max_position; i++) {
        assert(i < syncmers.size());
        // Method 3' skew sample more for prob exact matching
        std::bitset<64> b = (strobe1.hash ^ syncmers[i].hash) & parameters.q;
        uint64_t res = b.count();
        if (res < min_val) {
            min_val = res;
            strobe2 = syncmers[i];
        }
    }
    return make_randstrobe(strobe1, strobe2, parameters, mode);
}

Randstrobe RandstrobeGenerator::next() {
    while (syncmers.size() <= parameters.w_max) {
        Syncmer syncmer = syncmer_iterator.next();
        if (syncmer.is_end()) {
            break;
        }
        syncmers.push_back(syncmer);
    }
    if (syncmers.empty()) {
        return RandstrobeGenerator::end();
    }
    auto strobe1 = syncmers[0];
    auto max_position = strobe1.position + parameters.max_dist;
    uint64_t min_val = std::numeric_limits<uint64_t>::max();
    Syncmer strobe2 = strobe1;  // default if no nearby syncmer

    for (auto i = parameters.w_min;
         i < syncmers.size() && syncmers[i].position <= max_position; i++) {
        assert(i <= parameters.w_max);
        std::bitset<64> b = (strobe1.hash ^ syncmers[i].hash) & parameters.q;
        uint64_t res = b.count();
        if (res < min_val) {
            min_val = res;
            strobe2 = syncmers[i];
        }
    }
    syncmers.pop_front();
    return make_randstrobe(strobe1, strobe2, parameters, mode);
}

// --------------------------------------------------------------------------
// API
// --------------------------------------------------------------------------

std::vector<Randstrobe> randstrobes(std::string_view seq,
                                    const IndexParameters& parameters,
                                    RandstrobeHashMode mode) {
    std::vector<Randstrobe> result;
    if (seq.length() < parameters.randstrobe.w_max) {
        return result;
    }
    auto syncmers = canonical_syncmers(seq, parameters.syncmer);
    if (syncmers.empty()) {
        return result;
    }
    RandstrobeIterator it{syncmers, parameters.randstrobe, mode};
    while (it.has_next()) {
        result.push_back(it.next());
    }
    return result;
}

std::vector<Randstrobe> randstrobes(std::string_view seq, int read_length,
                                    RandstrobeHashMode mode) {
    return randstrobes(seq, IndexParameters::from_read_length(read_length), mode);
}

std::array<std::vector<QueryRandstrobe>, 2> randstrobes_query(
    std::string_view seq, const IndexParameters& parameters,
    RandstrobeHashMode mode) {
    std::array<std::vector<QueryRandstrobe>, 2> randstrobes;
    if (seq.length() < parameters.randstrobe.w_max) {
        return randstrobes;
    }

    // Syncmers for the forward sequence.
    auto syncmers = canonical_syncmers(seq, parameters.syncmer);
    if (syncmers.empty()) {
        return randstrobes;
    }

    // Forward randstrobes.
    RandstrobeIterator fwd{syncmers, parameters.randstrobe, mode};
    while (fwd.has_next()) {
        auto r = fwd.next();
        randstrobes[0].push_back(QueryRandstrobe{
            r.hash, r.hash_revcomp, r.strobe1_pos,
            r.strobe2_pos + static_cast<unsigned>(parameters.syncmer.k)});
    }

    // Canonical syncmers are invariant under reverse complement; reuse them and
    // only adjust coordinates and orientation.
    std::reverse(syncmers.begin(), syncmers.end());
    for (size_t i = 0; i < syncmers.size(); i++) {
        syncmers[i].position =
            seq.length() - syncmers[i].position - parameters.syncmer.k;
        syncmers[i].is_forward = !syncmers[i].is_forward;
    }

    // Randstrobes themselves cannot be reused for the reverse complement.
    RandstrobeIterator rc{syncmers, parameters.randstrobe, mode};
    while (rc.has_next()) {
        auto r = rc.next();
        randstrobes[1].push_back(QueryRandstrobe{
            r.hash, r.hash_revcomp, r.strobe1_pos,
            r.strobe2_pos + static_cast<unsigned>(parameters.syncmer.k)});
    }
    return randstrobes;
}

// --------------------------------------------------------------------------
// default parameter profiles (from cpp/indexparameters.cpp)
// --------------------------------------------------------------------------

namespace {
struct Profile {
    int canonical_read_length;
    int r_threshold;
    int k;
    int s_offset;
    int w_min;
    int w_max;
};
const Profile kProfiles[] = {
    {50, 70, 18, -4, 1, 4},   {75, 90, 20, -4, 1, 6},
    {100, 110, 20, -4, 2, 6}, {125, 135, 20, -4, 3, 8},
    {150, 175, 20, -4, 5, 11}, {250, 375, 22, -4, 6, 16},
    {400, std::numeric_limits<int>::max(), 23, -6, 5, 15},
};
}  // namespace

IndexParameters IndexParameters::from_read_length(int read_length, int k, int s,
                                                  int w_min, int w_max, int c,
                                                  int max_seed_len, int aux_len) {
    const int default_c = 8;
    int canonical_read_length = 50;
    for (const auto& p : kProfiles) {
        if (read_length <= p.r_threshold) {
            if (k == DEFAULT) k = p.k;
            if (s == DEFAULT) s = k + p.s_offset;
            if (w_min == DEFAULT) w_min = p.w_min;
            if (w_max == DEFAULT) w_max = p.w_max;
            canonical_read_length = p.canonical_read_length;
            break;
        }
    }

    int max_dist;
    if (max_seed_len == DEFAULT) {
        max_dist = std::clamp(canonical_read_length - 70, k, 255);
    } else {
        max_dist = max_seed_len - k;  // distance in start positions
    }
    uint64_t q = static_cast<uint64_t>(std::pow(2, c == DEFAULT ? default_c : c)) - 1;
    if (aux_len == DEFAULT) {
        aux_len = 17;
    }
    uint64_t main_hash_mask = ~0ULL << (9 + aux_len);

    return IndexParameters(SyncmerParameters(k, s),
                           RandstrobeParameters(q, max_dist, w_min, w_max,
                                                main_hash_mask));
}

// --------------------------------------------------------------------------
// Stream operators
// --------------------------------------------------------------------------

std::ostream& operator<<(std::ostream& os, const Syncmer& syncmer) {
    os << "Syncmer(hash=" << syncmer.hash << ", position=" << syncmer.position
       << ", is_forward=" << syncmer.is_forward << ")";
    return os;
}

std::ostream& operator<<(std::ostream& os, const Randstrobe& randstrobe) {
    os << "Randstrobe(hash=" << randstrobe.hash
       << ", strobe1_pos=" << randstrobe.strobe1_pos
       << ", strobe2_pos=" << randstrobe.strobe2_pos << ")";
    return os;
}

std::ostream& operator<<(std::ostream& os, const QueryRandstrobe& randstrobe) {
    os << "QueryRandstrobe(hash=" << randstrobe.hash << ", start=" << randstrobe.start
       << ", end=" << randstrobe.end << ")";
    return os;
}

}
