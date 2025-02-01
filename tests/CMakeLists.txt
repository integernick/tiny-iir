file(GLOB COMMON_TEST_SOURCES
        ${CMAKE_CURRENT_SOURCE_DIR}/src/common/test_*.cpp
)

set(TEST_SOURCES ${COMMON_TEST_SOURCES})

if(BUILD_FOR_ARM)
    file(GLOB FIXED_POINT_TEST_SOURCES
            ${CMAKE_CURRENT_SOURCE_DIR}/src/fixed-point/test_*.cpp
    )
    list(APPEND TEST_SOURCES ${FIXED_POINT_TEST_SOURCES})
endif()

add_executable(tiny_iir_tests ${TEST_SOURCES})

target_include_directories(tiny_iir_tests PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/include/common
)

if(BUILD_FOR_ARM)
    target_include_directories(tiny_iir_tests PRIVATE
            ${CMAKE_CURRENT_SOURCE_DIR}/../platform/arm
            ${CMAKE_CURRENT_SOURCE_DIR}/include/fixed-point
    )
else()
    target_include_directories(tiny_iir_tests PRIVATE
            ${CMAKE_CURRENT_SOURCE_DIR}/../platform/generic
    )
endif()

target_link_libraries(tiny_iir_tests
        PRIVATE
        tiny_iir_core
        GTest::gtest
        GTest::gtest_main
)

include(GoogleTest)
gtest_discover_tests(tiny_iir_tests)