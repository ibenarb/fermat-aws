// sqfilter.hpp
#pragma once
#include <vector>
#include <cstdint>

namespace sqf {

// Stell hier die Modulgröße ein.
// 0 -> Filter aus; 10080 ist unser Standard (2^5 * 3^2 * 5 * 7).
    constexpr unsigned M = 10080;

// Gibt true zurück, wenn r ein Quadratresiduum mod M ist.
// Implementierung: einmalig eine Maske der Quadrate mod M aufbauen.
    inline bool ok(unsigned r) {
        if (M == 0) return true;                  // Filter abgeschaltet
        static std::vector<std::uint8_t> mask;    // lazy init, threadsicher seit C++11
        if (mask.empty()) {
            mask.assign(M, 0);
            // Alle x (0..M-1) quadrieren und markieren: (x^2 mod M)
            for (unsigned x = 0; x < M; ++x) {
                unsigned s = static_cast<unsigned>((1ull * x * x) % M);
                mask[s] = 1;
            }
        }
        return mask[r % M] != 0;
    }

} // namespace sqf
