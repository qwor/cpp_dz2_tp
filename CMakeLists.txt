cmake_minimum_required(VERSION 3.20)

set(PROJECT_NAME "MatrixSample")

project(${PROJECT_NAME})

set(CMAKE_CXX_STANDARD 20)

enable_testing()
add_subdirectory(matrices)

message("MatrixLib_LIB = ${MatrixLib_LIB}")
message("MatrixLib_INCLUDE_DIRS = ${MatrixLib_INCLUDE_DIRS}")

add_executable(${PROJECT_NAME} src/main.cpp)
target_link_libraries(${PROJECT_NAME} PRIVATE ${MatrixLib_LIB})
target_include_directories(${PROJECT_NAME} PUBLIC ${MatrixLib_INCLUDE_DIRS})