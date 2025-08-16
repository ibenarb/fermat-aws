#pragma once
#include <string>
#include <sstream>
#include <cstring>

// Join argv into a single, reproducible command line for logging.
inline std::string join_argv_for_log(int argc, char** argv) {
    std::string s;
    s.reserve(256);
    for (int i = 0; i < argc; ++i) {
        if (i) s.push_back(' ');
        const char* a = argv[i];
        bool need_quotes = false;
        for (const char* p = a; *p; ++p) {
            if (*p == ' ' || *p == '\t' || *p == '"') { need_quotes = true; break; }
        }
        if (!need_quotes) {
            s += a;
        } else {
            s.push_back('"');
            for (const char* p = a; *p; ++p) {
                if (*p == '"') s += "\\\"";
                else s.push_back(*p);
            }
            s.push_back('"');
        }
    }
    return s;
}

// Parse a CLI option value with --key=value or --key value.
inline const char* find_cli_value(int argc, char** argv, const char* key) {
    const size_t klen = std::strlen(key);
    for (int i = 1; i < argc; ++i) {
        const char* a = argv[i];
        if (std::strncmp(a, "--", 2) != 0) continue;
        a += 2;
        if (std::strncmp(a, key, klen) == 0) {
            a += klen;
            if (*a == '=') return a + 1;
            if (*a == '\0' && i + 1 < argc) return argv[i + 1];
        }
    }
    return nullptr;
}
