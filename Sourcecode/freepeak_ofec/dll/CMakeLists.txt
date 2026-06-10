cmake_minimum_required(VERSION 3.20)

project(ofec_mpmmo_dll LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

get_filename_component(OFEC_ROOT "${CMAKE_CURRENT_LIST_DIR}/.." ABSOLUTE)
get_filename_component(REPO_ROOT "${OFEC_ROOT}/../.." ABSOLUTE)

file(GLOB_RECURSE OFEC_CORE_SOURCES CONFIGURE_DEPENDS
    "${OFEC_ROOT}/core/*.cpp"
)

file(GLOB_RECURSE OFEC_INSTANCE_SOURCES CONFIGURE_DEPENDS
    "${OFEC_ROOT}/instance/*.cpp"
)
list(FILTER OFEC_INSTANCE_SOURCES EXCLUDE REGEX "parameter_function\\.cpp$")
list(FILTER OFEC_INSTANCE_SOURCES EXCLUDE REGEX "one_peak_Weierstrass_rugeness\\.cpp$")

file(GLOB_RECURSE OFEC_UTILITY_SOURCES CONFIGURE_DEPENDS
    "${OFEC_ROOT}/utility/*.cpp"
)

set(OFEC_RUN_SOURCES
    "${OFEC_ROOT}/run/interface.cpp"
    "${OFEC_ROOT}/run/register_instance.cpp"
)

add_library(ofec_mpmmo SHARED
    "${CMAKE_CURRENT_LIST_DIR}/ofec_mpmmo_api.cpp"
    ${OFEC_CORE_SOURCES}
    ${OFEC_INSTANCE_SOURCES}
    ${OFEC_UTILITY_SOURCES}
    ${OFEC_RUN_SOURCES}
)

target_include_directories(ofec_mpmmo PRIVATE
    "${OFEC_ROOT}"
)

if(MSVC)
    target_compile_options(ofec_mpmmo PRIVATE /permissive- /FS /wd4819)
endif()

add_executable(FreePeaks_Benchmark
    "${OFEC_ROOT}/run/main_solver.cpp"
    ${OFEC_CORE_SOURCES}
    ${OFEC_INSTANCE_SOURCES}
    ${OFEC_UTILITY_SOURCES}
    ${OFEC_RUN_SOURCES}
)

target_include_directories(FreePeaks_Benchmark PRIVATE
    "${OFEC_ROOT}"
    "${OFEC_ROOT}/run"
)

target_compile_definitions(FreePeaks_Benchmark PRIVATE
    _CRT_SECURE_NO_WARNINGS
    OFEC_DIR="."
)

if(MSVC)
    target_compile_options(FreePeaks_Benchmark PRIVATE /permissive- /FS /wd4819)
endif()

set_target_properties(ofec_mpmmo PROPERTIES
    OUTPUT_NAME "ofec_mpmmo"
    RUNTIME_OUTPUT_DIRECTORY "${REPO_ROOT}/build/ofec_mpmmo/bin"
    LIBRARY_OUTPUT_DIRECTORY "${REPO_ROOT}/build/ofec_mpmmo/bin"
    ARCHIVE_OUTPUT_DIRECTORY "${REPO_ROOT}/build/ofec_mpmmo/lib"
    RUNTIME_OUTPUT_DIRECTORY_DEBUG "${REPO_ROOT}/build/ofec_mpmmo/bin"
    LIBRARY_OUTPUT_DIRECTORY_DEBUG "${REPO_ROOT}/build/ofec_mpmmo/bin"
    ARCHIVE_OUTPUT_DIRECTORY_DEBUG "${REPO_ROOT}/build/ofec_mpmmo/lib"
    RUNTIME_OUTPUT_DIRECTORY_RELEASE "${REPO_ROOT}/build/ofec_mpmmo/bin"
    LIBRARY_OUTPUT_DIRECTORY_RELEASE "${REPO_ROOT}/build/ofec_mpmmo/bin"
    ARCHIVE_OUTPUT_DIRECTORY_RELEASE "${REPO_ROOT}/build/ofec_mpmmo/lib"
)

set_target_properties(FreePeaks_Benchmark PROPERTIES
    OUTPUT_NAME "FreePeaks_Benchmark"
    RUNTIME_OUTPUT_DIRECTORY "${REPO_ROOT}/build/freepeaks_benchmark/bin"
    RUNTIME_OUTPUT_DIRECTORY_DEBUG "${REPO_ROOT}/build/freepeaks_benchmark/bin"
    RUNTIME_OUTPUT_DIRECTORY_RELEASE "${REPO_ROOT}/build/freepeaks_benchmark/bin"
)
