cmake_minimum_required(VERSION 3.9)

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
        ../src/dark_matter.h
        ../src/dark_matter.cpp
        external.cpp)

add_library(foo MODULE ${lib_SOURCES})

OPTION (USE_OpenMP "Use OpenMP to enamble <omp.h>" ON)

if(APPLE AND USE_OpenMP)
    if(CMAKE_C_COMPILER_ID MATCHES "Clang")
        set(OpenMP_C "${CMAKE_C_COMPILER}")
        set(OpenMP_C_FLAGS "-fopenmp=libomp -Wno-unused-command-line-argument")
        set(OpenMP_C_LIB_NAMES "libomp" "libgomp" "libiomp5")
        set(OpenMP_libomp_LIBRARY ${OpenMP_C_LIB_NAMES})
        set(OpenMP_libgomp_LIBRARY ${OpenMP_C_LIB_NAMES})
        set(OpenMP_libiomp5_LIBRARY ${OpenMP_C_LIB_NAMES})
    endif()
    if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
        set(OpenMP_CXX "${CMAKE_CXX_COMPILER}")
        set(OpenMP_CXX_FLAGS "-fopenmp=libomp -Wno-unused-command-line-argument")
        set(OpenMP_CXX_LIB_NAMES "libomp" "libgomp" "libiomp5")
        set(OpenMP_libomp_LIBRARY ${OpenMP_CXX_LIB_NAMES})
        set(OpenMP_libgomp_LIBRARY ${OpenMP_CXX_LIB_NAMES})
        set(OpenMP_libiomp5_LIBRARY ${OpenMP_CXX_LIB_NAMES})
    endif()
endif()

find_package(OpenMP REQUIRED)

target_link_libraries(foo PRIVATE OpenMP::OpenMP_CXX)


