cmake_minimum_required(VERSION 3.10)

SET(CMAKE_BUILD_TYPE Debug)

include_directories(${VelocityDispersion_SOURCE_DIR}/src)
link_directories(${VelocityDispersion_BINARY_DIR}/src)

FIND_PACKAGE(Boost COMPONENTS system filesystem unit_test_framework REQUIRED )

INCLUDE_DIRECTORIES( ${Boost_INCLUDE_DIR} )
set(CMAKE_CXX_STANDARD 14)

set ( VelocityDispersion_SOURCES
		integration.cpp
        integration.h
        main.cpp
        main.h
        utillity.cpp
        utillity.h
        dark_matter.cpp
        dark_matter.h
	galaxy.h
	galaxy.cpp)

set ( Executable_Name DirectExecution )

add_executable( ${Executable_Name} ${VelocityDispersion_SOURCES})

TARGET_LINK_LIBRARIES(${Executable_Name} LINK_PUBLIC
        ${Boost_LIBRARIES}
        ${Boost_FILESYSTEM_LIBRARY}
        ${Boost_SYSTEM_LIBRARY}
        )

