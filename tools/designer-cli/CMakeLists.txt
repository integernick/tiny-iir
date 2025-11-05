find_package(nlohmann_json QUIET)

if (NOT nlohmann_json_FOUND)
    include(FetchContent)
    FetchContent_Declare(
            nlohmann_json
            URL https://github.com/nlohmann/json/releases/download/v3.11.3/json.tar.xz
            DOWNLOAD_EXTRACT_TIMESTAMP TRUE
    )
    FetchContent_MakeAvailable(nlohmann_json)
endif ()

find_package(cxxopts CONFIG REQUIRED)

add_executable(tiny-iir-designer-cli designer_cli.cpp)

target_link_libraries(tiny-iir-designer-cli
        PRIVATE
        tiny_iir_core
        nlohmann_json::nlohmann_json
        cxxopts::cxxopts
)

install(TARGETS tiny-iir-designer-cli RUNTIME DESTINATION bin)
