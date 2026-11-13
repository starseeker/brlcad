cmake_minimum_required(VERSION 3.22)
project(libbobol_installed_consumer LANGUAGES CXX)

find_package(BRLCAD CONFIG REQUIRED COMPONENTS BObol)

add_executable(libbobol_installed_consumer main.cpp)
target_link_libraries(libbobol_installed_consumer PRIVATE BRLCAD::BObol)
