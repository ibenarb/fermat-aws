cd ~/code/fermat-aws

cat > README.md <<'MD'
# fermat-aws

Great—here’s a drop-in README.md you can paste into the repo. I kept commands as single-line one-liners and documented the flags you’ve actually used in this project.

fermat-aws

High-performance Fermat-based integer factorization with optional modular sieving and square-filter prefilters.
Targets both local/WSL and AWS x86-64-v3 builds.

Features

Parallel Fermat search (a := ⌈√N⌉ + offset) with lock-free work distribution.

Fast residue sieving by small primes; optional explicit prime set.

Square-filter to prune a whose a² − N cannot be a square modulo a product of small primes.

Periodic progress telemetry (time, thread-0 full tests/s, estimated a rate).

Memory-aware “large sieve” mode for very big caps.

Requirements

CMake ≥ 3.22

GCC ≥ 12 (GCC 13 tested)

GMP (headers & libs)

Ninja (recommended)

Ubuntu (WSL) deps (one-liner):
sudo apt-get update && sudo apt-get install -y build-essential cmake ninja-build g++ libgmp-dev

Build (CLI)

Native / WSL-optimized:
cmake -S . -B build/native -G Ninja -DCMAKE_BUILD_TYPE=Release -DFERMAT_ARCH=native -DFERMAT_LTO=ON && cmake --build build/native -j

AWS x86-64-v3 optimized:
cmake -S . -B build/aws-v3 -G Ninja -DCMAKE_BUILD_TYPE=Release -DFERMAT_ARCH=x86-64-v3 -DFERMAT_LTO=ON && cmake --build build/aws-v3 -j

Notes:

FERMAT_ARCH selects -march= (e.g., native, x86-64-v3).

FERMAT_LTO=ON enables link-time optimization.

Release builds are strongly recommended.

## Reporting

Each run appends to a text log with a timestamped header and an end-of-run block.

### Where the log goes
Resolution precedence:
1. `--report-file <PATH>` (CLI)
2. `FERMAT_REPORT_FILE` (environment)
3. Repo root detected by walking up from the executable’s directory looking for `.git` or `CMakeLists.txt`
4. Executable directory (fallback)

Example:
```bash
./build/clion-Release/fermat 10000000000360999999998869 --threads 2 --sieve-up-to 19
# or override:
FERMAT_REPORT_FILE=/tmp/fermat_runs.txt ./build/clion-Release/fermat ...
# or:
./build/clion-Release/fermat ... --report-file /tmp/fermat_runs.txt


Quick start

Small demo:
./fermat 10000000036999999769 --threads 24 --max-tests-per-thread 1000000000000 --sieve-up-to 17 --sieve-cap 300000

Heavier run (as used):
./fermat 1000000000003429999999999001 --threads 23 --max-tests-per-thread 18446744073709551613 --sieve-up-to 19 --sieve-cap 4000000000

Large-sieve stress (prints memory):
./fermat <N> --threads 61 --sieve-up-to 29 --sieve-cap 4000000000 --force-large-sieve

Command-line options
Option	Type	Meaning	Typical values / remarks
--threads	int	Worker threads	Use physical cores for best throughput; avoid SMT saturation if memory-bound
--max-tests-per-thread	uint64	Per-thread cap on a-tests	Use very large for “run until found” workloads
--sieve-up-to	int	Generate sieve primes up to this prime	Examples used: 17, 19, 23, 29
--sieve-primes	CSV of ints	Explicit small-prime set (overrides --sieve-up-to)	Example: 17,19,23,29
--sieve-cap	uint64	Upper bound influencing sieve period/range	Larger ⇒ more memory and fewer modulus misses
--force-large-sieve	flag	Force high-memory sieve and print memory telemetry	Use on big iron only
--sqfilter	int	Square-filter modulus (often product of small primes)	Examples: 15015, 20677; tunes skip-rate
--progress-interval	seconds	Telemetry cadence	If implemented; default prints every ~10 s

Implementation detail: choose --threads = T such that gcd(2T, M) = 1 where M is the sieve period; otherwise some offsets repeat residue classes and waste work. The program emits a hint like gcd(2*T, M) = 3 if you pick a poor T.

Algorithm sketch

We search x = ⌈√N⌉ + a for minimal a ≥ 0 with x² − N a perfect square.

Residue sieve: for primes p ≤ P, precompute admissible a mod p such that x² − N can be a quadratic residue; skip the rest.

Square-filter: pick a squareful modulus m (often a product of small primes) and remove a where (x² − N) mod m cannot be a square.

Parallelization: each thread advances its own a subsequence with minimal synchronization; optional CPU pinning on Linux.

Performance tips

Prefer --threads equal to the number of physical cores; verify the sieve note about gcd(2T, M).

--sieve-cap and --sieve-up-to trade memory for skip-rate; increases can help until the sieve becomes memory-bound.

--sqfilter k can help or hurt depending on k; benchmark a few co-prime products (e.g., 3·5·7·11·13).

On AWS, use instances with high single-core perf and ample RAM bandwidth; compile with -DFERMAT_ARCH=x86-64-v3.

CLion (WSL) quick setup

Toolchain: WSL, Ninja, GCC.

CMake options (Release): -DFERMAT_ARCH=native -DFERMAT_LTO=ON

Build target: fermat binary under your chosen build/... directory.

Examples we actually ran

./fermat 1000000000003429999999999001 --threads 23 --max-tests-per-thread 18446744073709551613 --sieve-up-to 19 --sieve-cap 4000000000

./fermat 1000000000003429999999999001 --threads 23 --max-tests-per-thread 18446744073709551613 --sieve-up-to 19 --sieve-cap 4000000000 --sqfilter 20677

./fermat 1000000000003429999999999001 --threads 23 --max-tests-per-thread 18446744073709551613 --sieve-primes 17,19,23,29 --sieve-cap 4000000000 --sqfilter 15015

Output anatomy

[sieve] M=... persistent≈... GiB, peak≈... GiB

[progress] N has d digits; T=threads; alpha≈...

[progress] t=... s; thread0_full=...; thread0_full/s=...; est_a_per_s≈...

Final status prints whether a factor was found and (optionally) the last a per core if enabled.

Development

Language: C++17

Core deps: GMP

Build system: CMake + Ninja

Style: prefer tabs for code indentation; avoid long line continuations in shell examples.

Roadmap

Optional OpenMP backend.

Better auto-tuning of sieve parameters.

Progress hooks for per-thread last-a reporting (lightweight atomics).

License

Add a LICENSE file (e.g., MIT or Apache-2.0).

Citation

If this tool helps your research, please cite the repository URL and commit hash