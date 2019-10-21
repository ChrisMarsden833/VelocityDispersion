cmake_minimum_required(VERSION 3.13)

include_directories(${VelocityDispersion_SOURCE_DIR}/src)
link_directories(${VelocityDispersion_BINARY_DIR}/src)

FIND_PACKAGE(Boost COMPONENTS system filesystem unit_test_framework REQUIRED )

INCLUDE_DIRECTORIES( ${Boost_INCLUDE_DIR} )
set(CMAKE_CXX_STANDARD 14)

set (lib_SOURCES
        ../src/utillity.h
        ../src/utillity.cpp
        ../src/Integration.h
        ../src/Integration.cpp
        ../src/desmond.h
        ../src/desmond.cpp
        external.cpp)

add_library(foo MODULE ${lib_SOURCES})