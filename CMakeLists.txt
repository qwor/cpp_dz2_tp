cmake_minimum_required(VERSION 3.18)

set(PROJECT_NAME "MatrixSample")

project(${PROJECT_NAME})

set(CMAKE_CXX_STANDARD 14)

enable_testing()
add_subdirectory(matrices)

add_executable(${PROJECT_NAME} src/main.cpp)
target_link_libraries(${PROJECT_NAME} INTERFACE ${MatricesLib_LIB})
target_include_directories(${PROJECT_NAME} PUBLIC ${MatricesLib_INCLUDE_DIRS})
