cmake_minimum_required(VERSION 3.21)
project(fast-iir)

set(CMAKE_CXX_STANDARD 20)

# Option to select platform
option(BUILD_FOR_ARM "Build for ARM platform using CMSIS-DSP" ON)

# Include common headers
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include/filter_design)
#include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include/chebyshev_type_1)

if(BUILD_FOR_ARM)
    # Build for ARM Cortex-M
    add_subdirectory(CMSIS-DSP)
    include_directories(${CMAKE_CURRENT_SOURCE_DIR}/platform/arm)
    include_directories(${CMAKE_CURRENT_SOURCE_DIR}/CMSIS-DSP/Include)

    # Define ARM_MATH_CM4 for Cortex-M4
    target_compile_definitions(CMSISDSP PUBLIC ARM_MATH_CM4)

    set(LIBS CMSISDSP)
else()
    # Build for POSIX (macOS arm64)
    include_directories(${CMAKE_CURRENT_SOURCE_DIR}/platform/posix)
    include_directories(${CMAKE_CURRENT_SOURCE_DIR}/CMSIS-DSP/Include)

    # Compile CMSIS-DSP source files
    file(GLOB_RECURSE CMSISDSP_SOURCES
            ${CMAKE_CURRENT_SOURCE_DIR}/CMSIS-DSP/Source/*.c
    )

    # Exclude ARM-specific source files if necessary
    # list(FILTER CMSISDSP_SOURCES EXCLUDE REGEX ".*arm_*_q15.c$")
    # Add more exclusions as needed

    add_library(CMSISDSP STATIC ${CMSISDSP_SOURCES})

    # Define ARM_MATH_AUTOVECTORIZE for generic C implementations
    target_compile_definitions(CMSISDSP PUBLIC ARM_MATH_AUTOVECTORIZE)

    set(LIBS CMSISDSP)
endif()

# Add common source files
set(SOURCES
        src/fast_iir.cpp
        # Add other source files as needed
)

add_executable(fast-iir src/fast_iir.cpp ${SOURCES})

# Link CMSIS-DSP library
target_link_libraries(fast-iir ${LIBS})