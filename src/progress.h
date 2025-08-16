#pragma once
#include <atomic>
#include <gmp.h>

// NOTE:
// mpz_srcptr is GMP's "const mpz_t" pointer type (typedef const __mpz_struct*).
// We pass N as mpz_srcptr (not as const mpz_t*), so call sites can pass N.get_mpz_t() directly.

struct ProgressConfig {
    std::atomic<uint64_t>*	counter_thread0;	// post-sieve counter from thread 0
    double					period_minutes;		// e.g., 60.0
    uint64_t				T;					// number of threads; stride is 2*T
    double					alpha;				// sieve survival fraction (H/L)
    mpz_srcptr				N;					// the big integer N (const)
    std::atomic<bool>*		stop_flag;			// set true to stop reporter
};

// Starts a detached reporter thread; zero work in worker hot paths.
void start_thread0_progress(const ProgressConfig& cfg);
