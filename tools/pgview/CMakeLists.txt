cmake_minimum_required(VERSION 3.14)
project(pgview CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

include(FetchContent)

# FTXUI — C++ TUI framework
FetchContent_Declare(ftxui
    GIT_REPOSITORY https://github.com/ArthurSonzogni/FTXUI.git
    GIT_TAG        v6.1.9
)

# nlohmann/json — JSON parsing
FetchContent_Declare(nlohmann_json
    GIT_REPOSITORY https://github.com/nlohmann/json.git
    GIT_TAG        v3.11.3
)

FetchContent_MakeAvailable(ftxui nlohmann_json)

add_executable(pgview main.cpp)
target_link_libraries(pgview PRIVATE
    ftxui::screen
    ftxui::dom
    ftxui::component
    nlohmann_json::nlohmann_json
)

install(TARGETS pgview RUNTIME DESTINATION bin)
