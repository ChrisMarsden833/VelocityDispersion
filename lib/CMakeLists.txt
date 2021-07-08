
cmake_minimum_required(VERSION 3.9)

include_directories(${VelocityDispersion_SOURCE_DIR}/src)
link_directories(${VelocityDispersion_BINARY_DIR}/src)

FIND_PACKAGE(Boost COMPONENTS system filesystem unit_test_framework REQUIRED )
INCLUDE_DIRECTORIES( ${Boost_INCLUDE_DIR} )

find_package(OpenMP REQUIRED)
find_package(GSL REQUIRED)

set(CMAKE_CXX_STANDARD 14)

set (lib_SOURCES
        ../src/utillity.h
        ../src/utillity.cpp
        ../src/integration.h
        ../src/integration.cpp
        ../src/sigmalib.h
        ../src/sigmalib.cpp
        ../src/bulge.h
        ../src/bulge.cpp
        ../src/disk.h
        ../src/disk.cpp
        ../src/halo.h
        ../src/halo.cpp
        ../src/blackhole.h
        ../src/blackhole.cpp
        external.cpp)

add_library(sigma MODULE ${lib_SOURCES})

target_link_libraries(sigma GSL::gsl GSL::gslcblas OpenMP::OpenMP_CXX)

