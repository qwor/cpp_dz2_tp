cmake_minimum_required(VERSION 3.20)
project(matrices_test)

# GoogleTest requires at least C++14
set(CMAKE_CXX_STANDARD 20)

include(FetchContent)
FetchContent_Declare(
        googletest
        GIT_REPOSITORY https://github.com/google/googletest.git
        GIT_TAG release-1.12.1
)
# For Windows: Prevent overriding the parent matrices's compiler/linker settings
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)

include(GoogleTest)

# Collects the names of all the source files in the specified directory
aux_source_directory(${CMAKE_CURRENT_SOURCE_DIR} SRC)
# file(GLOB TEST_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp)

add_executable(${PROJECT_NAME} ${SRC})

target_link_libraries(${PROJECT_NAME} PRIVATE ${MatricesLib_lib} GTest::gtest_main)
target_include_directories(${PROJECT_NAME} PUBLIC ${MatricesLib_INCLUDE_DIRS})

gtest_discover_tests(${PROJECT_NAME})