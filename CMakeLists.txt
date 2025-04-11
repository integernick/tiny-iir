cmake_minimum_required(VERSION 3.21)
project(tiny-iir LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 20)

option(BUILD_FOR_ARM "Build for ARM platform using CMSIS-DSP" ON)
option(TINY_IIR_BUILD_TESTS "Build tests" ON)

if(BUILD_FOR_ARM)
    add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/CMSIS-DSP CMSIS-DSP-build)
endif()

add_library(tiny_iir_core INTERFACE)

target_include_directories(tiny_iir_core
        INTERFACE
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:include>
)

if(BUILD_FOR_ARM)
    message(STATUS "Building for ARM platform")
    target_include_directories(tiny_iir_core
            INTERFACE
            ${CMAKE_CURRENT_SOURCE_DIR}/platform/arm
            ${CMAKE_CURRENT_SOURCE_DIR}/CMSIS-DSP/Include
    )
    target_link_libraries(tiny_iir_core INTERFACE CMSISDSP)
    target_compile_definitions(tiny_iir_core INTERFACE USE_CMSIS_DSP)
else()
    message(STATUS "Building for generic platform")
    target_include_directories(tiny_iir_core
            INTERFACE
            ${CMAKE_CURRENT_SOURCE_DIR}/../platform/generic
    )
endif()

if(TINY_IIR_BUILD_TESTS)
    find_package(GTest REQUIRED)
    add_subdirectory(tests)
endif()