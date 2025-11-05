# Collect sources
file(GLOB COMMON_TEST_SOURCES
        ${CMAKE_CURRENT_SOURCE_DIR}/src/common/test_*.cpp
)
set(TEST_SOURCES ${COMMON_TEST_SOURCES})

if (TINY_IIR_BUILD_WITH_CMSIS)
    file(GLOB FIXED_POINT_TEST_SOURCES
            ${CMAKE_CURRENT_SOURCE_DIR}/src/fixed-point/test_*.cpp
    )
    list(APPEND TEST_SOURCES ${FIXED_POINT_TEST_SOURCES})
endif()

add_executable(tiny_iir_tests ${TEST_SOURCES})

# Test includes (tests' own helpers)
target_include_directories(tiny_iir_tests PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/include/common
)

# Platform-specific includes for tests
if (TINY_IIR_BUILD_WITH_CMSIS)
    target_include_directories(tiny_iir_tests PRIVATE
            ${CMAKE_CURRENT_SOURCE_DIR}/../platform/arm
            ${CMAKE_CURRENT_SOURCE_DIR}/include/fixed-point
    )
else()
    target_include_directories(tiny_iir_tests PRIVATE
            ${CMAKE_CURRENT_SOURCE_DIR}/../platform/generic
    )
endif()

target_compile_definitions(tiny_iir_tests PRIVATE TINY_IIR_CASCADE_FILTER_DEBUG=1)
target_compile_features(tiny_iir_tests PRIVATE cxx_std_20)

# ---- GTest: find system copy, fall back to FetchContent if missing
include(FetchContent)
find_package(GTest QUIET)

if (NOT GTest_FOUND)
    FetchContent_Declare(
            googletest
            URL https://github.com/google/googletest/archive/refs/tags/v1.14.0.zip
    )
    set(BUILD_GMOCK OFF CACHE BOOL "" FORCE)
    set(INSTALL_GTEST OFF CACHE BOOL "" FORCE)
    FetchContent_MakeAvailable(googletest)
endif()

# Imported target names differ depending on how it’s provided
if (TARGET GTest::gtest_main)
    set(GTEST_MAIN_TGT GTest::gtest_main)
    set(GTEST_TGT      GTest::gtest)
elseif (TARGET GTest::Main)
    set(GTEST_MAIN_TGT GTest::Main)
    set(GTEST_TGT      GTest::GTest)
else()
    message(FATAL_ERROR "No GTest targets available (system & FetchContent both missing)")
endif()

target_link_libraries(tiny_iir_tests PRIVATE tiny_iir_core ${GTEST_TGT} ${GTEST_MAIN_TGT})

# Test discovery
include(GoogleTest)
gtest_discover_tests(tiny_iir_tests)