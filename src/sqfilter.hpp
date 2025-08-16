#pragma once
#include <cstdint>
#include <vector>

struct SqFilter {
    std::uint32_t mod{0};
    std::vector<std::uint64_t> bits;	// size = ceil(mod/64)

    static inline unsigned v2(std::uint32_t x){
        unsigned c=0; while ((x&1u)==0u) { x>>=1u; ++c; } return c;
    }

    // Initialize for modulus M; if odd_only==true, drop factors of 2.
    // Returns the effective modulus used (possibly M with powers of two removed).
    std::uint32_t init(std::uint32_t M, bool odd_only){
        if (M==0u) { mod=0; bits.clear(); return 0; }
        if (odd_only){
            unsigned t = v2(M);
            if (t>0) M >>= t;
        }
        mod = M;
        bits.assign((mod + 63u) / 64u, 0ull);
        if (mod==1u) { bits[0]=1ull; return mod; }

        for (std::uint32_t x=0; x<mod; ++x){
            std::uint32_t r = (std::uint32_t)((1ull*x*x) % mod);
            bits[r >> 6] |= (1ull << (r & 63u));
        }
        return mod;
    }

    inline bool enabled() const { return mod>0 && !bits.empty(); }

    inline bool is_qr(std::uint32_t r) const {
        if (!enabled()) return true;
        r %= mod;
        return ((bits[r >> 6] >> (r & 63u)) & 1ull) != 0ull;
    }
};
