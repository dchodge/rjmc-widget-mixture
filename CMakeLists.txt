cmake_minimum_required(VERSION 3.15)
project(rjmc_wasm)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# WebAssembly specific settings
set(CMAKE_SYSTEM_NAME WASM)
set(CMAKE_CXX_COMPILER "em++")
set(CMAKE_C_COMPILER "emcc")

# Emscripten specific flags
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -s WASM=1 -s EXPORTED_RUNTIME_METHODS='[\"ccall\",\"cwrap\"]' -s ALLOW_MEMORY_GROWTH=1 -s EXPORTED_FUNCTIONS='[\"_main\",\"_malloc\",\"_free\"]'")

# Include directories
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/src)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../rjmc/src/headers)

# Source files
set(SOURCES
    src/rjmc_wasm.cpp
    ../rjmc/src/rjmc_func.cpp
    ../rjmc/src/headers/rjmc.hpp
    ../rjmc/src/headers/mvn.hpp
)

# Create the WebAssembly module
add_executable(rjmc_wasm ${SOURCES})

# Link with Emscripten
target_link_libraries(rjmc_wasm)
