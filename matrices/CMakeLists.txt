cmake_minimum_required(VERSION 3.18)
project(MatricesLib)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} "-Wall -Wextra -Wpedantic -Werror")

file(GLOB_RECURSE HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/include/*.hpp ${CMAKE_CURRENT_SOURCE_DIR}/include/*.h)
message("HEADERS = ${HEADERS}")

set(LIB_NAME matrices)

# Collects the names of all the source files in the specified directory
add_library(${LIB_NAME} INTERFACE)
target_sources(${LIB_NAME} INTERFACE ${HEADERS})

target_include_directories(${LIB_NAME} INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/include>
        )

set(MatricesLib_LIB ${LIB_NAME})
set(MatricesLib_LIB ${MatricesLib_LIB} PARENT_SCOPE)
set(MatricesLib_INCLUDE_DIRS ${LIB_NAME} ${CMAKE_CURRENT_SOURCE_DIR}/include)
set(MatricesLib_INCLUDE_DIRS ${MatricesLib_INCLUDE_DIRS} PARENT_SCOPE)

enable_testing()
add_subdirectory(tests)
