cmake_minimum_required(VERSION 3.0.2)

option(test "Build all tests." OFF)
if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -stdlib=libc++ -O0 -g3 -Wall -pthread")
else()
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14 -O0 -g3 -Wall -pthread")
endif()

set(PROJECT_NAME Mathx)
project (${PROJECT_NAME})

# The version number.
set (${PROJECT_NAME}_VERSION_MAJOR 0)
set (${PROJECT_NAME}_VERSION_MINOR 1)

include_directories(./headers ./sources ./tests ./gtest-1.8.0)

set(SOURCES ./sources)
set(HEADER_FILES ./headers)
set(GTEST_DIR /home/dustin/Documents/googletest/googletest)

###########
# Testing #
###########
if(test)
  include_directories(${gtest_SOURCE_DIR}/include ${gtest_SOURCE_DIR} ./tests)
  set(TESTS ./tests)
  add_subdirectory(./gtest-1.8.0)

  enable_testing()

  ########
  # Unit #
  ########
  add_executable(${PROJECT_NAME}Tests ${TESTS}/tests.cpp)
  target_link_libraries(${PROJECT_NAME}Tests gtest gtest_main)
  add_test(NAME ${PROJECT_NAME}-tests COMMAND ${PROJECT_NAME}Tests)
endif()

# add the executable
add_executable(${PROJECT_NAME} main.cpp)
