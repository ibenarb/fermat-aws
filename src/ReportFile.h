#pragma once
#include <string>
#include <mutex>
#include <fstream>
#include <chrono>
#include <iomanip>
#include <sstream>
#include <cstdlib>
#include <filesystem>
#include <ctime>        // ← needed for gmtime_r / gmtime_s

#if defined(__linux__)
#include <unistd.h>   // readlink
#elif defined(__APPLE__)
#include <mach-o/dyld.h> // _NSGetExecutablePath
#elif defined(_WIN32)
    #include <windows.h>  // GetModuleFileNameA
#endif

class ReportFile {
public:
    // Optional: explicit init (safe to call multiple times).
    static void init_from_exe(const char* /*unused*/) {
        std::lock_guard<std::mutex> g(mu());
        ensure_initialized_unlocked();
    }

    static void set_path(const std::string& p) {
        std::lock_guard<std::mutex> g(mu());
        path() = p;
    }

    static std::string get_path() {
        std::lock_guard<std::mutex> g(mu());
        ensure_initialized_unlocked();
        return path();
    }

    // Thread-safe append (creates parent dir if missing).
    static void append_line(const std::string& line) {
        std::lock_guard<std::mutex> g(mu());
        ensure_initialized_unlocked();
        ensure_parent_dir_unlocked();
        std::ofstream out(path(), std::ios::out | std::ios::app);
        if (!out) return;
        out << line << '\n';
    }

    // Thread-safe append block (no extra newline).
    static void append_block(const std::string& block) {
        std::lock_guard<std::mutex> g(mu());
        ensure_initialized_unlocked();
        ensure_parent_dir_unlocked();
        std::ofstream out(path(), std::ios::out | std::ios::app);
        if (!out) return;
        out << block;
    }

    // UTC timestamp like 2025-08-16 14:32:05Z
    static std::string iso_utc_now() {
        using clock = std::chrono::system_clock;
        const auto now = clock::now();
        const std::time_t t = clock::to_time_t(now);
        std::tm tm{};
#if defined(_WIN32)
        gmtime_s(&tm, &t);
#else
        gmtime_r(&t, &tm);
#endif
        std::ostringstream os;
        os << std::put_time(&tm, "%Y-%m-%d %H:%M:%SZ");
        return os.str();
    }

private:
    static std::mutex& mu() {
        static std::mutex m;
        return m;
    }
    static std::string& path() {
        static std::string p; // lazily filled
        return p;
    }

    // Ensure path() is set exactly once, with precedence:
    // 1) env FERMAT_REPORT_FILE
    // 2) repo root (walk up from real exe dir searching .git or CMakeLists.txt)
    // 3) exe directory fallback
    static void ensure_initialized_unlocked() {
        if (!path().empty()) return;

        namespace fs = std::filesystem;
        std::error_code ec;

        // 1) ENV override
        if (const char* envp = std::getenv("FERMAT_REPORT_FILE")) {
            fs::path envpath = envp;
            if (!envpath.is_absolute()) {
                auto cwd = fs::current_path(ec);
                envpath = (cwd / envpath);
            }
            path() = envpath.lexically_normal().string();
            return;
        }

        // 2) Determine real executable directory
        fs::path exedir = get_exe_dir_unlocked();
        if (exedir.empty()) {
            exedir = fs::current_path(ec);
        }

        // Search for repo root markers upward
        fs::path best_root;
        for (fs::path p = exedir; !p.empty(); p = p.parent_path()) {
            std::error_code e2;
            bool has_git   = fs::exists(p / ".git", e2);
            bool has_cmake = fs::exists(p / "CMakeLists.txt", e2);
            if (has_git || has_cmake) { best_root = p; break; }
            if (p == p.root_path()) break;
        }

        if (!best_root.empty()) {
            path() = (best_root / "reports.txt").string();
        } else {
            // 3) Fallback next to the executable
            path() = (exedir / "reports.txt").string();
        }
    }

    // Real executable directory (platform-specific; ignores argv[0] & CWD).
    static std::filesystem::path get_exe_dir_unlocked() {
        namespace fs = std::filesystem;
#if defined(__linux__)
        char buf[4096];
        ssize_t n = readlink("/proc/self/exe", buf, sizeof(buf) - 1);
        if (n > 0) {
            buf[n] = '\0';
            fs::path exe(buf);
            return exe.has_parent_path() ? exe.parent_path() : fs::path();
        }
        return fs::path();
#elif defined(__APPLE__)
        uint32_t size = 0;
        _NSGetExecutablePath(nullptr, &size);
        std::string tmp(size, '\0');
        if (_NSGetExecutablePath(tmp.data(), &size) == 0) {
            fs::path exe(tmp.c_str());
            return exe.has_parent_path() ? exe.parent_path() : fs::path();
        }
        return fs::path();
    #elif defined(_WIN32)
        char buf[MAX_PATH];
        DWORD n = GetModuleFileNameA(nullptr, buf, MAX_PATH);
        if (n > 0) {
            fs::path exe(std::string(buf, buf + n));
            return exe.has_parent_path() ? exe.parent_path() : fs::path();
        }
        return fs::path();
    #else
        return fs::path();
#endif
    }

    static void ensure_parent_dir_unlocked() {
        namespace fs = std::filesystem;
        std::error_code ec;

        const std::string cur = path();  // snapshot
        fs::path pth{ cur };             // braces avoid most-vexing-parse
        fs::path parent = pth.has_parent_path() ? pth.parent_path() : fs::path(".");
        if (!fs::exists(parent, ec)) {
            fs::create_directories(parent, ec);
        }
    }
};
