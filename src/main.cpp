#include <gmpxx.h>
#include <gmp.h>
#include <atomic>
#include <chrono>
#include <cstdint>
#include <iostream>
#include <limits>
#include <mutex>
#include <string>
#include <thread>
#include <vector>
#include <algorithm>
#include <sstream>
#include <iomanip>

#include "progress.h"
#include "sqfilter.hpp"

// thread-0 progress counter and stop flag
std::atomic<uint64_t> g_thread0_tests{0};
static std::atomic<bool> g_stop_progress{false};

#if defined(__linux__)
#include <pthread.h>
#include <sched.h>
static void pin_to_core(unsigned cpu){
    cpu_set_t set; CPU_ZERO(&set); CPU_SET(cpu, &set);
    pthread_setaffinity_np(pthread_self(), sizeof(set), &set);
}
#endif

struct Result {
    bool found = false;
    bool limit_reached = false;		// hit --max-tests-per-thread
    bool time_limit = false;		// hit --seconds
    unsigned winner_thread = ~0u;
    std::uint64_t tries_to_find = 0;
    mpz_class p, q;
};

static inline bool is_square(const mpz_class& x) {
    return x >= 0 && mpz_perfect_square_p(x.get_mpz_t()) != 0;
}

static inline std::uint64_t gcd64(std::uint64_t a, std::uint64_t b){
    while (b) { std::uint64_t t = a % b; a = b; b = t; }
    return a;
}

// extended Euclid mod inverse for small p (int)
static long long modinv_ll(long long a, long long m){
    long long b = m, u = 1, v = 0;
    while (b) {
        long long t = a / b;
        a -= t * b; std::swap(a, b);
        u -= t * v; std::swap(u, v);
    }
    if (a != 1) return -1;
    if (u < 0) u += m;
    return u;
}

static std::vector<int> parse_prime_list(const std::string& s){
    std::vector<int> primes;
    std::stringstream ss(s);
    std::string tok;
    while (std::getline(ss, tok, ',')) {
        if (tok.empty()) continue;
        int p = std::stoi(tok);
        if (p >= 3) primes.push_back(p);
    }
    return primes;
}

static std::vector<int> odd_primes_up_to(int P){
    std::vector<int> ps;
    if (P < 3) return ps;
    std::vector<bool> sieve(static_cast<std::size_t>(P) + 1u, true);
    sieve[0] = sieve[1] = false;
    for (int i = 2; 1LL * i * i <= P; ++i) if (sieve[static_cast<std::size_t>(i)])
            for (long long j = 1LL * i * i; j <= P; j += i) sieve[static_cast<std::size_t>(j)] = false;
    for (int p = 3; p <= P; ++p) if (sieve[static_cast<std::size_t>(p)] && (p & 1)) ps.push_back(p);
    return ps;
}

// residues a mod p with (a^2 - N) quadratic residue mod p
static std::vector<int> allowed_residues_mod_p(const mpz_class& N, int p){
    std::vector<char> is_sq(static_cast<std::size_t>(p), 0);
    for (int x = 0; x < p; ++x) is_sq[static_cast<std::size_t>((1LL * x * x) % p)] = 1;
    unsigned long Np = mpz_fdiv_ui(N.get_mpz_t(), static_cast<unsigned long>(p));
    std::vector<int> allow;
    allow.reserve(static_cast<std::size_t>(p) / 2u + 2u);
    for (int a = 0; a < p; ++a) {
        int r = static_cast<int>((1LL * a * a - static_cast<long long>(Np)) % p);
        if (r < 0) r += p;
        if (is_sq[static_cast<std::size_t>(r)]) allow.push_back(a);
    }
    return allow;
}

// CRT combine residues modulo m with residues modulo p -> residues modulo m*p
static std::vector<std::uint32_t> crt_combine(const std::vector<std::uint32_t>& R, std::uint64_t m,
                                              const std::vector<int>& A, int p){
    std::vector<std::uint32_t> out;
    if (R.empty() || A.empty()) return out;
    auto Minv_mod_p = modinv_ll(static_cast<long long>(m % p), p);
    if (Minv_mod_p < 0) return {};
    out.reserve(R.size() * A.size());
    for (auto r : R) {
        int r_mod_p = static_cast<int>(r % static_cast<std::uint32_t>(p));
        for (int a : A) {
            int t = a - r_mod_p;
            if (t < 0) t += p;
            auto k = (1LL * t * Minv_mod_p) % p;
            auto x = static_cast<std::uint64_t>(r) + m * static_cast<std::uint64_t>(k);
            out.push_back(static_cast<std::uint32_t>(x));
        }
    }
    return out;
}

// density α = ∏_{p∈S} (p + χ_p(N)) / (2p), with χ_p the Legendre symbol (handle χ=0 -> α*=1)
static double sieve_alpha_legendre(const mpz_class& N, const std::vector<int>& primes){
    double a = 1.0;
    for (int p : primes) {
        if (p < 3 || (p & 1) == 0) continue;
        mpz_class P(p);
        int chi = mpz_legendre(N.get_mpz_t(), P.get_mpz_t()); // -1,0,1
        if (chi == 0) {
            // p | N: every a gives a square r=a^2 so all residues pass -> factor 1
            continue;
        }
        a *= (static_cast<double>(p + chi)) / (2.0 * static_cast<double>(p));
    }
    return a;
}

int main(int argc, char** argv) {
    if (argc < 2) {
        std::cerr << "usage: fermat N [--threads T] [--max-tests-per-thread K] [--seconds S] [--affinity] "
                     "[--mod2k K] [--sieve-up-to P | --sieve-primes a,b,c] [--sieve-cap M] [--force-large-sieve] [--sqfilter M]\n";
        return 2;
    }

    // -------- parse inputs --------
    mpz_class N;
    if (N.set_str(argv[1], 10) != 0) {
        std::cerr << "N: invalid decimal\n";
        return 2;
    }

    unsigned T = std::thread::hardware_concurrency();
    std::uint64_t max_tests_per_thread = 0;
    double seconds_limit = 0.0;
    bool use_affinity = false;

    unsigned mod2k = 0;					// OFF by default
    std::vector<int> sieve_primes;		// OFF by default
    std::uint64_t sieve_cap = 4000000000ull;	// raised default cap (~3.73 GiB address span in entries)
    bool force_large_sieve = false;		// require explicit opt-in for huge tables
    unsigned sqfilter_M = 0;			// 0=OFF; e.g. 315 or 10080

    for (int i = 2; i < argc; i++) {
        std::string s(argv[i]);
        if (s == "--threads" && i + 1 < argc) {
            T = std::stoul(argv[++i]);
        } else if (s == "--max-tests-per-thread" && i + 1 < argc) {
            max_tests_per_thread = std::stoull(argv[++i]);
        } else if (s == "--seconds" && i + 1 < argc) {
            seconds_limit = std::stod(argv[++i]);
        } else if (s == "--affinity") {
            use_affinity = true;
        } else if (s == "--mod2k" && i + 1 < argc) {
            mod2k = std::stoul(argv[++i]);
            if (mod2k != 0 && mod2k < 6) mod2k = 6;
            if (mod2k > 16) mod2k = 16;
        } else if (s == "--sieve-up-to" && i + 1 < argc) {
            int P = std::stoi(argv[++i]);
            if (P >= 3) sieve_primes = odd_primes_up_to(P);
        } else if (s == "--sieve-primes" && i + 1 < argc) {
            sieve_primes = parse_prime_list(argv[++i]);
        } else if (s == "--sieve-cap" && i + 1 < argc) {
            sieve_cap = std::stoull(argv[++i]);
            if (sieve_cap < 10000) sieve_cap = 10000;
            if (sieve_cap > 50000000000ull) sieve_cap = 50000000000ull; // hard stop ~50e9 residues
        } else if (s == "--force-large-sieve") {
            force_large_sieve = true;
        } else if (s == "--sqfilter" && i + 1 < argc) {
            sqfilter_M = static_cast<unsigned>(std::stoul(argv[++i]));
        }
    }
    if (T == 0) T = 1;

    // -------- precompute mod 2^k (only if enabled) --------
    const bool prefilter_on = (mod2k >= 6 && mod2k <= 16);
    const unsigned K = prefilter_on ? mod2k : 0;
    std::vector<std::uint64_t> sqmask;
    std::uint32_t M2 = 0;
    if (prefilter_on) {
        M2 = (K >= 31 ? 0x7FFFFFFFu : ((1u << K) - 1u));
        const std::size_t bits = 1u << K;
        sqmask.assign((bits + 63u) / 64u, 0ull);
        for (std::uint32_t x = 0; x < (1u << K); ++x) {
            std::uint32_t r = (x * x) & M2;
            sqmask[r >> 6] |= (1ull << (r & 63u));
        }
    }

    // -------- a0 = ceil(sqrt(N)); parity constraint --------
    mpz_class a0;
    mpz_sqrt(a0.get_mpz_t(), N.get_mpz_t());
    if (a0 * a0 < N) ++a0;

    unsigned long n_mod4 = mpz_fdiv_ui(N.get_mpz_t(), 4);
    bool need_a_odd  = (n_mod4 == 1);
    bool need_a_even = (n_mod4 == 3);

    // -------- build odd-prime sieve (optional) --------
    std::uint64_t sieve_mod = 1;
    std::vector<std::uint32_t> allowed_residues_modM; // residues a mod M that survive
    bool sieve_on = false;

    if (!sieve_primes.empty()) {
        std::sort(sieve_primes.begin(), sieve_primes.end());
        sieve_primes.erase(std::unique(sieve_primes.begin(), sieve_primes.end()), sieve_primes.end());

        allowed_residues_modM = {0};
        for (int p : sieve_primes) {
            if (p < 3 || (p % 2 == 0)) continue;
            if (!force_large_sieve && sieve_mod * static_cast<std::uint64_t>(p) > 50000000ull) {
                std::cerr << "[sieve] cap 50000000 reached without --force-large-sieve; dropping primes >= " << p << "\n";
                break;
            }
            if (sieve_mod * static_cast<std::uint64_t>(p) > sieve_cap) {
                std::cerr << "[sieve] cap " << sieve_cap << " reached; dropping primes >= " << p << "\n";
                break;
            }
            auto allow_p = allowed_residues_mod_p(N, p);
            if (allow_p.empty()) { allowed_residues_modM.clear(); break; }
            allowed_residues_modM = crt_combine(allowed_residues_modM, sieve_mod, allow_p, p);
            sieve_mod *= static_cast<std::uint64_t>(p);
            if (allowed_residues_modM.empty()) break;
        }
        if (!allowed_residues_modM.empty() && sieve_mod > 1) sieve_on = true;
    }

    // -------- compute density α for progress (approximate) --------
    double alpha = 1.0;
    if (sieve_on) alpha = sieve_alpha_legendre(N, sieve_primes);
    if (alpha <= 0.0) alpha = 1.0; // safety

    // -------- threading state --------
    std::atomic<bool> stop_all{false};
    Result res;
    std::mutex res_mx;

    std::vector<std::thread> workers;
    workers.reserve(T);
    std::vector<std::uint64_t> tries_by_thread(T, 0);

    // fast constants
    const unsigned long stride_ui = 2ul * T;
    const unsigned long s2_ui = stride_ui * stride_ui;
    const unsigned long delta_inc_ui = 2ul * s2_ui;

    // build jump table if sieve_on
    std::vector<unsigned> jump;			// jump[r] = minimal m>=1 to reach next allowed; 0 = dead cycle
    std::vector<char> allowed_mask;
    std::uint64_t stride_modM = 0;

    // Pre-print estimate if we know M
    if (sieve_on) {
        const long double M = static_cast<long double>(sieve_mod);
        long double persistent_est_bytes = M * (sizeof(unsigned) + sizeof(char)); // 5 bytes per residue
        long double temp_peak_est_bytes = M * 13.0L; // empirical upper bound of our implementation
        auto toGiB = [](long double b){ return static_cast<double>(b / (1024.0L*1024.0L*1024.0L)); };
        std::cerr << std::fixed << std::setprecision(6)
                  << "[sieve] M=" << sieve_mod
                  << " persistent≈" << toGiB(persistent_est_bytes) << " GiB,"
                  << " peak≈" << toGiB(temp_peak_est_bytes) << " GiB\n";
    }

    // -------- SQ filter (declare & init outside worker) --------
    SqFilter SQ;
    unsigned sqM_used = 0;
    if (sqfilter_M > 0) {
        const bool odd_only = prefilter_on; // if mod2k active, strip 2-adic factors
        sqM_used = SQ.init(sqfilter_M, odd_only);
    }

    // -------- start progress reporter (thread-0 only) --------
    {
        ProgressConfig pcfg{
                &g_thread0_tests,				// counter_thread0
                1.0,							// report every 10 minutes
                T,								// here T is both threads and AP stride factor (stride_ui = 2*T)
                alpha,							// sieve density
                N.get_mpz_t(),					// mpz_srcptr to N (see progress.h)
                &g_stop_progress
        };
        start_thread0_progress(pcfg);
    }

    auto t0 = std::chrono::steady_clock::now();

    // Build jump table with timing + actual memory accounting
    double jump_build_secs = 0.0;
    std::uint64_t persistent_actual_bytes = 0;
    std::uint64_t temp_peak_actual_bytes = 0;

    if (sieve_on) {
        auto jb_start = std::chrono::steady_clock::now();

        allowed_mask.assign(static_cast<std::size_t>(sieve_mod), 0);
        for (auto r : allowed_residues_modM) allowed_mask[static_cast<std::size_t>(r)] = 1;

        jump.assign(static_cast<std::size_t>(sieve_mod), 0);
        persistent_actual_bytes =
                static_cast<std::uint64_t>(allowed_mask.size() * sizeof(char)) +
                static_cast<std::uint64_t>(jump.size() * sizeof(unsigned));

        stride_modM = (sieve_mod ? static_cast<std::uint64_t>(stride_ui) % sieve_mod : 0);
        auto g = gcd64(sieve_mod, (stride_modM ? stride_modM : sieve_mod));	// g>=1 if sieve_mod>0
        if (g > 1) {
            std::cerr << "[sieve] note: gcd(2*T, M) = " << g
                      << " (M=" << sieve_mod << "). Choose T with gcd=1 for best speed.\n";
        }
        auto cycle_len = sieve_mod / g;

        std::vector<std::uint32_t> cycle; cycle.reserve(static_cast<std::size_t>(cycle_len));
        std::vector<unsigned> dist;	   dist.reserve(static_cast<std::size_t>(cycle_len));
        std::vector<std::size_t> pos;	 pos.reserve(static_cast<std::size_t>(cycle_len / 2 + 1));

        auto account_peak = [&](){
            std::uint64_t bytes =
                    static_cast<std::uint64_t>(cycle.capacity() * sizeof(std::uint32_t)) +
                    static_cast<std::uint64_t>(dist.capacity()  * sizeof(unsigned)) +
                    static_cast<std::uint64_t>(pos.capacity()   * sizeof(std::size_t));
            if (bytes > temp_peak_actual_bytes) temp_peak_actual_bytes = bytes;
        };
        account_peak();

        for (std::uint64_t c = 0; c < g; ++c) {
            cycle.clear();
            auto r = c % sieve_mod;
            do {
                cycle.push_back(static_cast<std::uint32_t>(r));
                r += stride_modM; if (r >= sieve_mod) r -= sieve_mod;
            } while (r != (c % sieve_mod));

            const std::size_t L = cycle.size();
            dist.assign(L, 0u);
            pos.clear();

            // grow capacities predictably, update temp peak
            if (cycle.capacity() < L) cycle.reserve(L);
            if (dist.capacity()  < L) dist.reserve(L);
            if (pos.capacity()   < L/2 + 1) pos.reserve(L/2 + 1);
            account_peak();

            for (std::size_t i = 0; i < L; ++i)
                if (allowed_mask[static_cast<std::size_t>(cycle[i])]) pos.push_back(i);

            if (pos.empty()) {
                for (std::size_t i = 0; i < L; ++i) jump[static_cast<std::size_t>(cycle[i])] = 0;
                continue;
            }

            for (std::size_t k = 0; k < pos.size(); ++k) {
                std::size_t s = pos[k];
                std::size_t t = pos[(k + 1) % pos.size()];
                std::size_t gap = (t >= s ? (t - s) : (L - (s - t)));
                for (std::size_t j = 0; j < gap; ++j) {
                    std::size_t idx = (s + j) % L;
                    std::size_t d = gap - j;		// minimal m>=1 from idx to next allowed
                    dist[idx] = static_cast<unsigned>(d);
                }
            }
            for (std::size_t i = 0; i < L; ++i) jump[static_cast<std::size_t>(cycle[i])] = dist[i];
        }

        auto jb_end = std::chrono::steady_clock::now();
        jump_build_secs = std::chrono::duration<double>(jb_end - jb_start).count();

        auto toGiB64 = [](std::uint64_t b){ return static_cast<double>(static_cast<long double>(b) / (1024.0L*1024.0L*1024.0L)); };
        std::cerr << std::fixed << std::setprecision(6)
                  << "[sieve] build_time_sec=" << jump_build_secs
                  << " persistent_actual≈" << toGiB64(persistent_actual_bytes) << " GiB,"
                  << " temp_peak_actual≈" << toGiB64(temp_peak_actual_bytes) << " GiB\n";
    }

    for (unsigned i = 0; i < T; i++) {
        workers.emplace_back([&, i]() {
#if defined(__linux__)
            if (use_affinity) pin_to_core(i % std::thread::hardware_concurrency());
#endif
            mpz_class a = a0;
            mpz_add_ui(a.get_mpz_t(), a.get_mpz_t(), 2ul * i);
            bool a_is_odd = mpz_odd_p(a.get_mpz_t());
            if (need_a_odd && !a_is_odd) { mpz_add_ui(a.get_mpz_t(), a.get_mpz_t(), 1); }
            if (need_a_even &&  a_is_odd) { mpz_add_ui(a.get_mpz_t(), a.get_mpz_t(), 1); }

            std::uint64_t r_modM = 0;
            if (sieve_on) r_modM = mpz_fdiv_ui(a.get_mpz_t(), static_cast<unsigned long>(sieve_mod));

            // SQ per-thread state
            std::uint32_t a_modSq = 0, N_modSq = 0, strideSq = 0;
            if (SQ.enabled()) {
                a_modSq = static_cast<std::uint32_t>(mpz_fdiv_ui(a.get_mpz_t(), static_cast<unsigned long>(SQ.mod)));
                N_modSq = static_cast<std::uint32_t>(mpz_fdiv_ui(N.get_mpz_t(), static_cast<unsigned long>(SQ.mod)));
                strideSq = static_cast<std::uint32_t>(stride_ui % SQ.mod);
            }

            mpz_class d = a * a - N;

            mpz_class delta;
            mpz_mul_ui(delta.get_mpz_t(), a.get_mpz_t(), stride_ui);
            mpz_mul_2exp(delta.get_mpz_t(), delta.get_mpz_t(), 1);
            mpz_add_ui(delta.get_mpz_t(), delta.get_mpz_t(), s2_ui);

            std::uint64_t local = 0;

            // advance by mm>=1 arithmetic steps (counts all skipped positions)
            auto advance_mm = [&](unsigned mm){
                if (mm == 0) return;
                local += mm;
                mpz_add_ui(a.get_mpz_t(), a.get_mpz_t(), static_cast<unsigned long>(mm * stride_ui));
                mpz_addmul_ui(d.get_mpz_t(), delta.get_mpz_t(), static_cast<unsigned long>(mm));
                auto t = static_cast<unsigned long>(mm);
                auto add2 = static_cast<unsigned long>(s2_ui * static_cast<unsigned long>(t * (t - 1)));
                mpz_add_ui(d.get_mpz_t(), d.get_mpz_t(), add2);
                mpz_add_ui(delta.get_mpz_t(), delta.get_mpz_t(), static_cast<unsigned long>(mm * delta_inc_ui));
                if (sieve_on) {
                    auto adv = (stride_modM * static_cast<std::uint64_t>(mm)) % sieve_mod;
                    r_modM += adv; if (r_modM >= sieve_mod) r_modM -= sieve_mod;
                }
                if (SQ.enabled()) {
                    std::uint64_t advSq = (static_cast<std::uint64_t>(strideSq) * mm) % SQ.mod;
                    a_modSq = static_cast<std::uint32_t>((a_modSq + advSq) % SQ.mod);
                }
            };

            while (!stop_all.load(std::memory_order_relaxed)) {
                // periodic time check
                if (seconds_limit > 0.0 && ((local & 0xFFFFu) == 0)) {
                    double secs_now = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
                    if (secs_now >= seconds_limit) {
                        tries_by_thread[i] = local;
                        {
                            std::scoped_lock lk(res_mx);
                            if (!res.found) res.time_limit = true;
                        }
                        stop_all.store(true, std::memory_order_relaxed);
                        return;
                    }
                }

                if (sieve_on) {
                    // ensure we are at an allowed residue; if not, jump to the next allowed
                    if (!allowed_mask[static_cast<std::size_t>(r_modM)]) {
                        unsigned m = jump[static_cast<std::size_t>(r_modM)];
                        if (m == 0) { // dead cycle for this thread
                            tries_by_thread[i] = local;
                            return;
                        }
                        advance_mm(m); // land on allowed residue
                    }

                    // optional SQ prefilter modulo sqfilter_M (before mod2k)
                    if (SQ.enabled() && mpz_sgn(d.get_mpz_t()) >= 0) {
                        std::uint32_t rM = static_cast<std::uint32_t>(( (1ull*a_modSq*a_modSq) + SQ.mod - N_modSq ) % SQ.mod);
                        if (!SQ.is_qr(rM)) {
                            unsigned m_next = jump[static_cast<std::size_t>(r_modM)];
                            advance_mm(m_next); // skip to next allowed
                            if (max_tests_per_thread && local >= max_tests_per_thread) {
                                { std::scoped_lock lk(res_mx); if (!res.found) res.limit_reached = true; }
                                tries_by_thread[i] = local; stop_all.store(true, std::memory_order_relaxed); return;
                            }
                            continue;
                        }
                    }

                    // optional mod 2^K prefilter (for the allowed candidate)
                    if (K > 0 && mpz_sgn(d.get_mpz_t()) >= 0) {
                        unsigned long low = mpz_get_ui(d.get_mpz_t());
                        auto rbits = static_cast<std::uint32_t>(low & M2);
                        if (((sqmask[rbits >> 6] >> (rbits & 63u)) & 1ull) == 0ull) {
                            // skip this allowed candidate, jump to next allowed
                            unsigned m_next = jump[static_cast<std::size_t>(r_modM)];
                            advance_mm(m_next); // counts current + skipped
                            if (max_tests_per_thread && local >= max_tests_per_thread) {
                                {
                                    std::scoped_lock lk(res_mx);
                                    if (!res.found) res.limit_reached = true;
                                }
                                tries_by_thread[i] = local;
                                stop_all.store(true, std::memory_order_relaxed);
                                return;
                            }
                            continue;
                        }
                    }

                    // full test at current allowed residue (post-sieve)
                    if (i == 0) g_thread0_tests.fetch_add(1, std::memory_order_relaxed);
                    if (d >= 0 && is_square(d)) {
                        mpz_class b;
                        mpz_sqrt(b.get_mpz_t(), d.get_mpz_t());
                        mpz_class p = a - b, q = a + b;
                        if (p > 0 && q > 0 && p * q == N) {
                            {
                                std::scoped_lock lk(res_mx);
                                if (!res.found) {
                                    res.found = true;
                                    res.p = p; res.q = q;
                                    res.winner_thread = i;
                                    res.tries_to_find = local + 1;
                                }
                            }
                            tries_by_thread[i] = local + 1;
                            stop_all.store(true, std::memory_order_relaxed);
                            return;
                        }
                    }

                    // move to next allowed residue
                    {
                        unsigned m_next = jump[static_cast<std::size_t>(r_modM)];
                        advance_mm(m_next);
                    }

                    // per-thread cap?
                    if (max_tests_per_thread && local >= max_tests_per_thread) {
                        {
                            std::scoped_lock lk(res_mx);
                            if (!res.found) res.limit_reached = true;
                        }
                        tries_by_thread[i] = local;
                        stop_all.store(true, std::memory_order_relaxed);
                        return;
                    }
                    continue; // back to top
                }

                // -------- no sieve path: 1 candidate per loop --------
                if (SQ.enabled() && mpz_sgn(d.get_mpz_t()) >= 0) {
                    std::uint32_t rM = static_cast<std::uint32_t>(( (1ull*a_modSq*a_modSq) + SQ.mod - N_modSq ) % SQ.mod);
                    if (!SQ.is_qr(rM)) {
                        advance_mm(1u);
                        if (max_tests_per_thread && local >= max_tests_per_thread) {
                            { std::scoped_lock lk(res_mx); if (!res.found) res.limit_reached = true; }
                            tries_by_thread[i] = local; stop_all.store(true, std::memory_order_relaxed); return;
                        }
                        continue;
                    }
                }

                if (K > 0 && mpz_sgn(d.get_mpz_t()) >= 0) {
                    unsigned long low = mpz_get_ui(d.get_mpz_t());
                    auto rbits = static_cast<std::uint32_t>(low & M2);
                    if (((sqmask[rbits >> 6] >> (rbits & 63u)) & 1ull) == 0ull) {
                        advance_mm(1u);
                        if (max_tests_per_thread && local >= max_tests_per_thread) {
                            {
                                std::scoped_lock lk(res_mx);
                                if (!res.found) res.limit_reached = true;
                            }
                            tries_by_thread[i] = local;
                            stop_all.store(true, std::memory_order_relaxed);
                            return;
                        }
                        continue;
                    }
                }

                // full check (no sieve)
                if (i == 0) g_thread0_tests.fetch_add(1, std::memory_order_relaxed);
                if (d >= 0 && is_square(d)) {
                    mpz_class b;
                    mpz_sqrt(b.get_mpz_t(), d.get_mpz_t());
                    mpz_class p = a - b, q = a + b;
                    if (p > 0 && q > 0 && p * q == N) {
                        {
                            std::scoped_lock lk(res_mx);
                            if (!res.found) {
                                res.found = true;
                                res.p = p; res.q = q;
                                res.winner_thread = i;
                                res.tries_to_find = local + 1;
                            }
                        }
                        tries_by_thread[i] = local + 1;
                        stop_all.store(true, std::memory_order_relaxed);
                        return;
                    }
                }

                advance_mm(1u);
                if (max_tests_per_thread && local >= max_tests_per_thread) {
                    {
                        std::scoped_lock lk(res_mx);
                        if (!res.found) res.limit_reached = true;
                    }
                    tries_by_thread[i] = local;
                    stop_all.store(true, std::memory_order_relaxed);
                    return;
                }
            }

            tries_by_thread[i] = local;
        });
    }

    for (auto& th : workers) th.join();

    // stop reporter
    g_stop_progress.store(true, std::memory_order_relaxed);

    auto t1 = std::chrono::steady_clock::now();
    double secs = std::chrono::duration<double>(t1 - t0).count();

    std::uint64_t sum = 0;
    for (auto v : tries_by_thread) sum += v;

    if (res.found) {
        std::cout << "FOUND\n";
        std::cout << "p=" << res.p << "\n";
        std::cout << "q=" << res.q << "\n";
        std::cout << "winner_thread=" << res.winner_thread << "\n";
        std::cout << "tries_to_find=" << res.tries_to_find << "\n";
    } else if (res.limit_reached) {
        std::cout << "NOT FOUND (per-thread test limit reached)\n";
    } else if (seconds_limit > 0.0) {
        std::cout << "TIME LIMIT\n";
    } else {
        std::cout << "ABORTED\n";
    }

    std::cout << "tests_sum=" << sum << "\n";
    std::cout << "tests_per_sec=" << (secs > 0 ? (static_cast<double>(sum) / secs) : 0.0) << "\n";
    std::cout << "time_sec=" << secs << "\n";
    std::cout << "mod2k=" << (prefilter_on ? static_cast<int>(K) : 0) << "\n";
    std::cout << "sqfilter_M=" << (SQ.enabled() ? static_cast<unsigned>(sqM_used) : 0) << "\n";
    if (!sieve_primes.empty()) {
        std::cout << "sieve_primes=";
        for (std::size_t i = 0; i < sieve_primes.size(); ++i) {
            if (i) std::cout << ",";
            std::cout << sieve_primes[i];
        }
        std::cout << "\n";
    }
    std::cout << "sieve_mod=" << (sieve_on ? sieve_mod : 0) << "\n";
    return res.found ? 0 : 1;
}
