cmake_minimum_required(VERSION 3.10)

# set the project name and version
project(transform-SIMD VERSION 1.0)

# C++ setting
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)
add_compile_options(-Wall -Wextra -pedantic -Werror  -O2 -mavx2 -g)

list(APPEND CMAKE_PREFIX_PATH $ENV{HOME}/.local/)
find_package(benchmark REQUIRED)
find_package (Threads)


add_executable(RkJacobian_test 
	             RkJacobian_test.cxx)

add_executable(RkJacobian_bench
	             RkJacobian_bench.cxx)

target_link_libraries(RkJacobian_test)
target_link_libraries(RkJacobian_bench benchmark::benchmark)

