find_package(nlohmann_json CONFIG REQUIRED)
find_package(cxxopts CONFIG REQUIRED)

add_executable(tiny-iir-designer-cli designer_cli.cpp)

target_link_libraries(tiny-iir-designer-cli
        PRIVATE
        tiny_iir_core
        nlohmann_json::nlohmann_json
        cxxopts::cxxopts
)

install(TARGETS tiny-iir-designer-cli RUNTIME DESTINATION bin)
