cmake_minimum_required(VERSION 3.22)
project(fermat LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_BUILD_TYPE Release)

add_executable(fermat
		main.cpp
		progress.cpp
)

# Compile & link flags (GCC/Clang)
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU" OR CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
	add_compile_options(-O3 -pipe -flto -march=native -mtune=native)
	target_link_options(fermat PRIVATE -flto -Wl,-O1 -Wl,--as-needed)
endif()

find_library(GMP_LIB gmp REQUIRED)
find_library(GMPXX_LIB gmpxx REQUIRED)
target_link_libraries(fermat PRIVATE ${GMPXX_LIB} ${GMP_LIB} pthread)
