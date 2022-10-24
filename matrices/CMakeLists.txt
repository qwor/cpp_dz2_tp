cmake_minimum_required(VERSION 3.20)

project("MatrixLib")

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} "-Wall -Wextra -Wpedantic -Werror")

file(GLOB_RECURSE HEADERS matrices/include/*.hpp matrices/include/*.h)

message("HEADERS = ${HEADERS}")

set(LIB_NAME matrices)

add_library(${LIB_NAME} ${HEADERS} INTERFACE)
target_include_directories(${LIB_NAME} INTERFACE include)

set(MatrixLib_LIB ${LIB_NAME})
set(MatrixLib_LIB ${MatrixLib_LIB} PARENT_SCOPE)
set(MatrixLib_INCLUDE_DIRS ${LIB_NAME} include)
set(MatrixLib_INCLUDE_DIRS ${MatrixLib_INCLUDE_DIRS} PARENT_SCOPE)

message("MatrixLib_LIB = ${MatrixLib_LIB}")
message("MatrixLib_INCLUDE_DIRS = ${MatrixLib_INCLUDE_DIRS}")

enable_testing()
add_subdirectory(tests)
