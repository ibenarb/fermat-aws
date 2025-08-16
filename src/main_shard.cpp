#include <gmpxx.h>
#include <atomic>
#include <chrono>
#include <cstdint>
#include <iostream>
#include <limits>
#include <mutex>
#include <string>
#include <thread>
#include <vector>

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
    while(b){ std::uint64_t t=a%b; a=b; b=t; }
    return a;
}
static inline std::uint64_t lcm64(std::uint64_t a, std::uint64_t b){
    if(a==0 || b==0) return 0;
    return (a / gcd64(a,b)) * b;
}

int main(int argc, char** argv) {
    if (argc < 2) {
        std::cerr << "usage: fermat N [--threads T] [--max-tests-per-thread K] [--seconds S] [--affinity] [--stride S] [--offset O]\n";
        return 2;
    }

    mpz_class N;
    if (N.set_str(argv[1], 10) != 0) {
        std::cerr << "N: invalid decimal\n";
        return 2;
    }

    unsigned T = std::thread::hardware_concurrency();
    std::uint64_t max_tests_per_thread = 0;			// 0 = unlimited
    double seconds_limit = 0.0;						// 0 = unlimited
    bool use_affinity = false;
    std::uint64_t user_stride = 0;					// 0 => disabled
    std::uint64_t user_offset = 0;

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
        } else if (s == "--stride" && i + 1 < argc) {
            user_stride = std::stoull(argv[++i]);
        } else if (s == "--offset" && i + 1 < argc) {
            user_offset = std::stoull(argv[++i]);
        }
    }
    if (T == 0) T = 1;

    // a0 = ceil(sqrt(N))
    mpz_class a0;
    mpz_sqrt(a0.get_mpz_t(), N.get_mpz_t());
    if (a0 * a0 < N) ++a0;

    // Parity filter: ensure a^2 ≡ N (mod 4)
    unsigned long n_mod4 = mpz_fdiv_ui(N.get_mpz_t(), 4);
    bool need_a_odd  = (n_mod4 == 1);
    bool need_a_even = (n_mod4 == 3);

    // Compose local parity stride (2*T) with cluster stride using LCM.
    std::uint64_t base_stride = 2ull * T;
    std::uint64_t shard_stride = (user_stride ? user_stride : 1ull);
    std::uint64_t final_stride = lcm64(base_stride, shard_stride);

    // Normalize offset into [0, final_stride)
    std::uint64_t norm_offset = (final_stride ? (user_offset % final_stride) : 0ull);

    // Adjust offset once to satisfy required parity of 'a'
    {
        mpz_class a_test = a0;
        if (norm_offset) {
            mpz_class off(norm_offset);
            mpz_add(a_test.get_mpz_t(), a_test.get_mpz_t(), off.get_mpz_t());
        }
        bool odd = mpz_odd_p(a_test.get_mpz_t());
        if ((need_a_odd && !odd) || (need_a_even && odd)) {
            ++norm_offset;
            if (final_stride) norm_offset %= final_stride;
        }
    }

    std::atomic<bool> stop_all{false};
    Result res;
    std::mutex res_mx;

    std::vector<std::thread> workers;
    workers.reserve(T);
    std::vector<std::uint64_t> tries_by_thread(T, 0);

    // Small-const fast path if stride fits in unsigned long
    const bool ui_ok =
            final_stride > 0 &&
            final_stride <= static_cast<std::uint64_t>(std::numeric_limits<unsigned long>::max());

    const unsigned long stride_ui = ui_ok ? static_cast<unsigned long>(final_stride) : 0ul;
    const unsigned long s2_ui = ui_ok ? static_cast<unsigned long>(stride_ui * stride_ui) : 0ul;
    const unsigned long delta_inc_ui = ui_ok ? static_cast<unsigned long>(2ul * s2_ui) : 0ul;

    auto t0 = std::chrono::steady_clock::now();

    for (unsigned i = 0; i < T; i++) {
        workers.emplace_back([&, i]() {
#if defined(__linux__)
            if (use_affinity) pin_to_core(i % std::thread::hardware_concurrency());
#endif
            mpz_class a = a0;

            // start offset: cluster norm_offset + per-thread spread (2*i)
            if (norm_offset) {
                mpz_class off(norm_offset);
                mpz_add(a.get_mpz_t(), a.get_mpz_t(), off.get_mpz_t());
            }
            mpz_add_ui(a.get_mpz_t(), a.get_mpz_t(), 2ul * i);

            // Ensure parity (should already match)
            bool a_is_odd = mpz_odd_p(a.get_mpz_t());
            if (need_a_odd && !a_is_odd) { mpz_add_ui(a.get_mpz_t(), a.get_mpz_t(), 1); }
            if (need_a_even &&  a_is_odd) { mpz_add_ui(a.get_mpz_t(), a.get_mpz_t(), 1); }

            mpz_class d = a * a - N;

            mpz_class delta;
            mpz_class s2_big, stride_big, two_s2_big;

            if (ui_ok) {
                // delta = 2*a*stride + s2   (all *_ui ops)
                mpz_mul_ui(delta.get_mpz_t(), a.get_mpz_t(), stride_ui);
                mpz_mul_2exp(delta.get_mpz_t(), delta.get_mpz_t(), 1);
                mpz_add_ui(delta.get_mpz_t(), delta.get_mpz_t(), s2_ui);
            } else {
                // Big-int fallback
                stride_big = mpz_class(final_stride);
                s2_big = stride_big * stride_big;
                two_s2_big = 2 * s2_big;

                mpz_mul(delta.get_mpz_t(), a.get_mpz_t(), stride_big.get_mpz_t());
                mpz_mul_2exp(delta.get_mpz_t(), delta.get_mpz_t(), 1);
                mpz_add(delta.get_mpz_t(), delta.get_mpz_t(), s2_big.get_mpz_t());
            }

            std::uint64_t local = 0;

            while (!stop_all.load(std::memory_order_relaxed)) {
                // Time limit check every 65536 iterations
                if (seconds_limit > 0.0 && ((local & 0xFFFFu) == 0)) {
                    double secs = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
                    if (secs >= seconds_limit) {
                        tries_by_thread[i] = local;
                        {
                            std::scoped_lock lk(res_mx);
                            if (!res.found) res.time_limit = true;
                        }
                        stop_all.store(true, std::memory_order_relaxed);
                        return;
                    }
                }

                // Test current candidate
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

                // Optional per-thread limit
                ++local;
                if (max_tests_per_thread && local >= max_tests_per_thread) {
                    {
                        std::scoped_lock lk(res_mx);
                        if (!res.found) res.limit_reached = true;
                    }
                    tries_by_thread[i] = local;
                    stop_all.store(true, std::memory_order_relaxed);
                    return;
                }

                // Advance: d += delta; delta += 2*s2; a += stride
                if (ui_ok) {
                    mpz_add(d.get_mpz_t(), d.get_mpz_t(), delta.get_mpz_t());
                    mpz_add_ui(delta.get_mpz_t(), delta.get_mpz_t(), delta_inc_ui);
                    mpz_add_ui(a.get_mpz_t(), a.get_mpz_t(), stride_ui);
                } else {
                    mpz_add(d.get_mpz_t(), d.get_mpz_t(), delta.get_mpz_t());
                    mpz_add(delta.get_mpz_t(), delta.get_mpz_t(), two_s2_big.get_mpz_t());
                    mpz_add(a.get_mpz_t(), a.get_mpz_t(), stride_big.get_mpz_t());
                }
            }

            tries_by_thread[i] = local;
        });
    }

    for (auto& th : workers) th.join();

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
    } else if (res.time_limit) {
        std::cout << "TIME LIMIT\n";
    } else {
        std::cout << "ABORTED\n";
    }
    std::cout << "stride=" << (final_stride ? final_stride : base_stride) << "\n";
    std::cout << "offset=" << norm_offset << "\n";
    std::cout << "tests_sum=" << sum << "\n";
    std::cout << "tests_per_sec=" << (secs > 0 ? (sum / secs) : 0.0) << "\n";
    std::cout << "time_sec=" << secs << "\n";
    return res.found ? 0 : 1;
}
