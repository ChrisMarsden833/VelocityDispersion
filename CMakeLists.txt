cmake_minimum_required(VERSION 3.13)
project(VelocityDispersion)

enable_testing()

set(CMAKE_CXX_STANDARD 14)

FIND_PACKAGE(Boost COMPONENTS system filesystem unit_test_framework REQUIRED )
INCLUDE_DIRECTORIES( ${Boost_INCLUDE_DIR} )
LINK_DIRECTORIES ( ${Boost_LIBRARY_DIRS} )

add_executable(VelocityDispersion
        desmond.cpp
        desmond.h
        integration.cpp
        integration.h
        main.cpp
        main.h
        utillity.cpp
        utillity.h
        Testing/Testing.cpp)

TARGET_LINK_LIBRARIES(VelocityDispersion LINK_PUBLIC
        ${Boost_LIBRARIES}
        ${Boost_FILESYSTEM_LIBRARY}
        ${Boost_SYSTEM_LIBRARY}
        ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
        #boost_unit_test_framework )





