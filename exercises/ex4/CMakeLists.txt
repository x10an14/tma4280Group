project(exercise4)
cmake_minimum_required(VERSION 2.6)

add_subdirectory(../../examples/common common)

add_definitions(${COMMON_DEFINITIONS})
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${COMMON_C_FLAGS} -std=c99 -Wall")
include_directories(${COMMON_INCLUDE_PATHS})

add_executable(ex4 ex4.c)
target_link_libraries(ex4 ${COMMON_LIBRARIES})
