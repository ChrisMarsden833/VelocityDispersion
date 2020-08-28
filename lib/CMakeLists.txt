cmake_minimum_required(VERSION 3.9)

include_directories(${VelocityDispersion_SOURCE_DIR}/src)
link_directories(${VelocityDispersion_BINARY_DIR}/src)

FIND_PACKAGE(Boost COMPONENTS system filesystem unit_test_framework REQUIRED )
INCLUDE_DIRECTORIES( ${Boost_INCLUDE_DIR} )

find_package(OpenMP REQUIRED)


set(CMAKE_CXX_STANDARD 14)

set (lib_SOURCES
        ../src/utillity.h
        ../src/utillity.cpp
        ../src/integration.h
        ../src/integration.cpp
		../src/dark_matter.h
        ../src/dark_matter.cpp
	../src/galaxy.h
	../src/galaxy.cpp
        external.cpp)

add_library(sigma MODULE ${lib_SOURCES})

target_link_libraries(sigma PRIVATE  OpenMP::OpenMP_CXX)

