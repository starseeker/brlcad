cmake_minimum_required(VERSION 3.20)
project(fltkcad_package_consumer LANGUAGES CXX)

find_package(BRLCAD CONFIG REQUIRED COMPONENTS fltkcad)

add_executable(package_consumer main.cpp)
target_link_libraries(package_consumer PRIVATE BRLCAD::fltkcad)
