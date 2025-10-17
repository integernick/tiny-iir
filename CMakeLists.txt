cmake_minimum_required(VERSION 3.21)
project(tiny-iir LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 20)

option(TINY_IIR_BUILD_WITH_CMSIS "Build tiny-iir using CMSIS-DSP" ON)
option(TINY_IIR_BUILD_TESTS "Build tiny-iir tests" ON)
option(TINY_IIR_BUILD_TOOLS "Build tiny-iir CLI utility" ON)
set(TINY_IIR_CHUNK_SIZE 32 CACHE STRING "Chunk size for processing")

if(TINY_IIR_BUILD_WITH_CMSIS)
    add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/CMSIS-DSP CMSIS-DSP-build)
endif()

add_library(tiny_iir_core INTERFACE)

target_include_directories(tiny_iir_core
        INTERFACE
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:include>
)

set(TINY_IIR_CMSIS_CORE_DIR "" CACHE PATH
        "Directory that contains cmsis_compiler.h (CMSIS Core headers)")

if(TINY_IIR_BUILD_WITH_CMSIS)
    message(STATUS "Building tiny-iir with CMSIS-DSP")
    target_include_directories(tiny_iir_core
            INTERFACE
            ${CMAKE_CURRENT_SOURCE_DIR}/platform/arm
            ${CMAKE_CURRENT_SOURCE_DIR}/CMSIS-DSP/Include
    )
    target_link_libraries(tiny_iir_core INTERFACE CMSISDSP)
    target_include_directories(CMSISDSP PRIVATE ${TINY_IIR_CMSIS_CORE_DIR})
else()
    message(STATUS "Building tiny-iir for generic platform")
    target_include_directories(tiny_iir_core
            INTERFACE
            ${CMAKE_CURRENT_SOURCE_DIR}/platform/generic
    )
endif()

target_compile_definitions(tiny_iir_core INTERFACE TINY_IIR_CHUNK_SIZE=${TINY_IIR_CHUNK_SIZE})

if(TINY_IIR_BUILD_TESTS)
    find_package(GTest REQUIRED)
    add_subdirectory(tests)
endif()

if(TINY_IIR_BUILD_TOOLS)
    add_subdirectory(tools/designer-cli)
endif()
