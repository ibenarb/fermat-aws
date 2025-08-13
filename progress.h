#pragma once
#include <gmp.h>
#include <atomic>
#include <cstdint>

struct ProgressConfig {
	std::atomic<std::uint64_t>* counter_thread0;	// increments only on thread 0 full tests
	double report_every_sec;						// period to print progress
	unsigned threads;								// total worker threads (T)
	double sieve_density;							// alpha ≈ surviving fraction after sieve (1.0 if off)
	mpz_srcptr Nptr;								// N for context
	std::atomic<bool>* stop_flag;					// set true to stop reporter
};

void start_thread0_progress(const ProgressConfig& cfg);
