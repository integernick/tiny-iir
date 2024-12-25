cmake_minimum_required(VERSION 3.21)
project(tiny-iir)

set(CMAKE_CXX_STANDARD 20)

find_package(GTest REQUIRED)

# Option to select platform
option(BUILD_FOR_ARM "Build for ARM platform using CMSIS-DSP" ON)

# Include common headers
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)

if(BUILD_FOR_ARM)
    # Build for ARM Cortex-M
    add_subdirectory(CMSIS-DSP)
    include_directories(${CMAKE_CURRENT_SOURCE_DIR}/platform/arm)
    include_directories(${CMAKE_CURRENT_SOURCE_DIR}/CMSIS-DSP/Include)

    # Define ARM_MATH_CM4 for Cortex-M4
    #target_compile_definitions(CMSISDSP PUBLIC ARM_MATH_CM4)

    set(LIBS CMSISDSP)
else()
    # Build for POSIX (macOS arm64)
    include_directories(${CMAKE_CURRENT_SOURCE_DIR}/platform/posix)
    include_directories(${CMAKE_CURRENT_SOURCE_DIR}/CMSIS-DSP/Include)

    # Compile CMSIS-DSP source files
    file(GLOB_RECURSE CMSISDSP_SOURCES
            ${CMAKE_CURRENT_SOURCE_DIR}/CMSIS-DSP/Source/*.c
    )

    add_library(CMSISDSP STATIC ${CMSISDSP_SOURCES})

    # Define ARM_MATH_AUTOVECTORIZE for generic C implementations
    target_compile_definitions(CMSISDSP PUBLIC ARM_MATH_AUTOVECTORIZE)

    set(LIBS CMSISDSP)
endif()

# Add common source files
set(SOURCES
        src/tiny_iir.cpp
)

add_executable(tiny-iir-tests
        tests/test_filters.cpp
)

target_link_libraries(tiny-iir-tests
        PRIVATE
        GTest::gtest
        GTest::gtest_main
        ${LIBS}
)

include(GoogleTest)
gtest_discover_tests(tiny-iir-tests)