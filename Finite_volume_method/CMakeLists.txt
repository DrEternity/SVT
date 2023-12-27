cmake_minimum_required(VERSION 3.2)
project(volumes)

find_package(inmost REQUIRED)

link_directories(${INMOST_LIBRARY_DIRS})
include_directories(${INMOST_INCLUDE_DIRS})
add_definitions(${INMOST_DEFINITIONS})

add_executable(main main.cpp)
target_link_libraries(main ${INMOST_LIBRARIES})
