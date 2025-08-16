#include "progress.h"
#include <chrono>
#include <cmath>
#include <cstdio>
#include <thread>

void start_thread0_progress(const ProgressConfig& cfg) {
    std::thread([=]() {
        // Good enough for progress lines (approximate double)
        double sqrtN = std::sqrt(mpz_get_d(cfg.N));

        // Each surviving candidate advances ~ (2*T)/alpha in 'a'
        double scale = (2.0 * static_cast<double>(cfg.T)) / cfg.alpha;

        using clock = std::chrono::steady_clock;
        auto t0 = clock::now();
        auto period = std::chrono::duration<double>(cfg.period_minutes * 60.0);
        auto next = t0 + std::chrono::duration_cast<clock::duration>(period);

        for (;;) {
            std::this_thread::sleep_until(next);
            if (cfg.stop_flag && cfg.stop_flag->load(std::memory_order_relaxed)) break;

            uint64_t k = cfg.counter_thread0->load(std::memory_order_relaxed);	// post-sieve tests (thread 0)
            auto now = clock::now();
            double elapsed = std::chrono::duration<double>(now - t0).count();
            double rate = (elapsed > 0.0) ? (static_cast<double>(k) / elapsed) : 0.0;

            // Coverage (density-based approximation, thread 0 only)
            double delta_a = scale * static_cast<double>(k);						// Δa ≈ (2T/α) * k
            double delta_max = std::sqrt(8.0 * sqrtN * delta_a);					// δ_max ≈ sqrt(8 sqrtN Δa)

            std::fprintf(stderr,
                         "[progress] t=%.0fs k0=%llu rate=%.3ge/s Δa≈%.3g δ_max≈%.3g\n",
                         elapsed,
                         static_cast<unsigned long long>(k),
                         rate,
                         delta_a,
                         delta_max
            );

            next += std::chrono::duration_cast<clock::duration>(period);
        }
    }).detach();
}
