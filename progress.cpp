#include "progress.h"
#include <atomic>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <gmp.h>
#include <iostream>
#include <thread>

static inline int digits10(mpz_srcptr n){
	return mpz_sizeinbase(n, 10);
}

void start_thread0_progress(const ProgressConfig& cfg){
	std::thread([cfg](){
		using clock = std::chrono::steady_clock;
		const double period = (cfg.report_every_sec > 0.1 ? cfg.report_every_sec : 10.0);

		std::uint64_t prev = cfg.counter_thread0 ? cfg.counter_thread0->load(std::memory_order_relaxed) : 0;
		auto t0 = clock::now();
		auto last = t0;

		const int Nd = digits10(cfg.Nptr);
		std::cerr << "[progress] N has " << Nd << " digits; T=" << cfg.threads
				  << "; alpha≈" << (cfg.sieve_density > 0 ? cfg.sieve_density : 1.0) << "\n";

		while(cfg.stop_flag && !cfg.stop_flag->load(std::memory_order_relaxed)){
			std::this_thread::sleep_for(std::chrono::duration<double>(period));

			auto now = clock::now();
			double dt = std::chrono::duration<double>(now - last).count();
			if (dt <= 0.0) continue;

			std::uint64_t cur = cfg.counter_thread0 ? cfg.counter_thread0->load(std::memory_order_relaxed) : 0;
			std::uint64_t d = (cur >= prev ? (cur - prev) : 0);

			// thread0 full-tests per second
			double t0_full_per_s = d / dt;

			// crude estimate of overall arithmetic-steps ("a" positions) per second:
			// scale to all threads, and undo sieve density to approximate pre-sieve rate
			double alpha = (cfg.sieve_density > 0.0 ? cfg.sieve_density : 1.0);
			double est_a_per_s = (t0_full_per_s * static_cast<double>(cfg.threads)) / alpha;

			double t_total = std::chrono::duration<double>(now - t0).count();

			std::cerr << "[progress] t=" << t_total
					  << " s; thread0_full=" << cur
					  << "; thread0_full/s=" << t0_full_per_s
					  << "; est_a_per_s≈" << est_a_per_s
					  << "\n";

			prev = cur;
			last = now;
		}
		std::cerr << "[progress] reporter stopped.\n";
	}).detach();
}
