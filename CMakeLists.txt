cmake_minimum_required(VERSION 3.25)

project(NuclearSim VERSION 1.0)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_FLAGS "-O3 -lncurses")

file(GLOB_RECURSE SOURCES src/*.cpp)

add_executable(NuclearSim ${SOURCES})
target_link_libraries(NuclearSim PUBLIC stdc++ m)

