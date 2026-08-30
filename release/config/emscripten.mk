# Compiler to use for C++ files, we use emscripten's wrapper
CXX             := em++
# Compiler to use for C files (used for emebedded dependencies)
CC              := emcc

# Suggested optimized build flags for Emscripten Compiler
CXXFLAGS        := -O3 -Wall -std=c++20
# We ignore this warning as it is not relevant for Emscripten, but we have code written in a way
# that triggers it
CXXFLAGS        += -Wno-missing-template-arg-list-after-template-kw

# optional linker flags
LDFLAGS         := -fPIC -g0
LDFLAGS         += -s ALLOW_MEMORY_GROWTH -s MAXIMUM_MEMORY=4294967296
LDFLAGS         += -s ENVIRONMENT="web"
LDFLAGS         += -s EXPORTED_RUNTIME_METHODS=ccall,cwrap,UTF8ToString
LDFLAGS         += -fwasm-exceptions
LDFLAGS         += -lembind
LDFLAGS         += -s WASM=1
# For further documentation consult the following:
# https://github.com/emscripten-core/emscripten/blob/main/src/settings.js
# https://emscripten.org/docs/tools_reference/emcc.html

#CXXFLAGS        += -s PTHREAD_POOL_SIZE=4

# Parallelization mode, must be: OFF
# See e.g. `config.git.mk` for further documentation.
PARALLEL_MODE   := NONE

# optional MPI and OpenMP flags
# No MPI flags, as it is not supported in Emscripten
# No OpenMP unless you are targeting threads within WASM
MPIFLAGS        :=
OMPFLAGS        :=

# Option: CPU_SISD
# See e.g. `config.git.mk` for further documentation.
# CPU_SISD must always be present.
PLATFORMS       := CPU_SISD

# Fundamental arithmetic data type
# Common options are float or double
FLOATING_POINT_TYPE := double

# Any entries are passed to the compiler as `-DFEATURE_*` declarations.
# Used to enable some alternative code paths and dependencies.
# In this case, we are using Emscripten-specific code paths.
FEATURES        := EMSCRIPTEN

# Set to OFF if libz and tinyxml2 are provided by the system (optional)
USE_EMBEDDED_DEPENDENCIES := ON
