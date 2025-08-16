cmake_minimum_required(VERSION 3.22)
project(Fermat LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 20)

find_package(PkgConfig REQUIRED)
pkg_check_modules(GMP REQUIRED gmp)
find_package(Threads REQUIRED)

# Tunable optimizations (no code duplication)
option(FERMAT_LTO "Enable LTO/IPO in Release" ON)
set(FERMAT_ARCH "native" CACHE STRING "CPU tune: native|znver2|x86-64|x86-64-v2|x86-64-v3")
set_property(CACHE FERMAT_ARCH PROPERTY STRINGS native znver2 x86-64 x86-64-v2 x86-64-v3)

# Build a LIST of arch flags (so they expand to separate options)
set(_arch_opts "")
if(FERMAT_ARCH STREQUAL "native")
  list(APPEND _arch_opts -march=native -mtune=native)
elseif(FERMAT_ARCH STREQUAL "znver2")
  list(APPEND _arch_opts -march=znver2 -mtune=znver2)
elseif(FERMAT_ARCH STREQUAL "x86-64")
  list(APPEND _arch_opts -march=x86-64)
elseif(FERMAT_ARCH STREQUAL "x86-64-v2")
  list(APPEND _arch_opts -march=x86-64-v2)
elseif(FERMAT_ARCH STREQUAL "x86-64-v3")
  list(APPEND _arch_opts -march=x86-64-v3)
endif()

# --- Targets ---
add_executable(fermat
  src/main.cpp
  src/progress.cpp
  src/progress.h
  src/sqfilter.hpp
)
target_include_directories(fermat PRIVATE ${GMP_INCLUDE_DIRS})
target_link_libraries(fermat PRIVATE gmpxx ${GMP_LIBRARIES} Threads::Threads)
target_compile_options(fermat PRIVATE
  -Wall -Wextra -Wpedantic
  $<$<CONFIG:Release>:-O3>
  $<$<CONFIG:Release>:-pipe>
  $<$<CONFIG:Release>:${_arch_opts}>
)

add_executable(fermat_shard src/main_shard.cpp)
target_include_directories(fermat_shard PRIVATE ${GMP_INCLUDE_DIRS})
target_link_libraries(fermat_shard PRIVATE gmpxx ${GMP_LIBRARIES} Threads::Threads)
target_compile_options(fermat_shard PRIVATE
  -Wall -Wextra -Wpedantic
  $<$<CONFIG:Release>:-O3>
  $<$<CONFIG:Release>:-pipe>
  $<$<CONFIG:Release>:${_arch_opts}>
)

# LTO/IPO after targets exist
include(CheckIPOSupported)
check_ipo_supported(RESULT ipo_ok OUTPUT ipo_msg)
if(FERMAT_LTO AND ipo_ok)
  set_property(TARGET fermat PROPERTY INTERPROCEDURAL_OPTIMIZATION_RELEASE TRUE)
  set_property(TARGET fermat_shard PROPERTY INTERPROCEDURAL_OPTIMIZATION_RELEASE TRUE)
elseif(FERMAT_LTO AND NOT ipo_ok)
  message(STATUS "IPO/LTO unavailable: ${ipo_msg}")
endif()
