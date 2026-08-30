# Example build config for OpenLB using Intel C++ on HoreKa

CXX             := icpx
CC              := icx

CXXFLAGS        := -O3 -Wall -mavx2 -mavx512f
CXXFLAGS        += -std=c++20 --gcc-toolchain=/opt/gcc/13/

LDFLAGS += --gcc-toolchain=/opt/gcc/13/ 

PARALLEL_MODE   := NONE

PLATFORMS       := CPU_SISD CPU_SIMD

USE_EMBEDDED_DEPENDENCIES := ON
