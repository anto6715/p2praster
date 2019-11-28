cmake_minimum_required(VERSION 3.15)
project(p2praster)

set(CMAKE_CXX_STANDARD 14)

add_executable(p2praster main.cpp raster.cpp raster.h error.cpp error.h statistic.cpp statistic.h)
target_link_libraries(p2praster igraph)
FIND_PACKAGE(Boost COMPONENTS program_options REQUIRED)
INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIRS})
TARGET_LINK_LIBRARIES(p2praster ${Boost_LIBRARIES})