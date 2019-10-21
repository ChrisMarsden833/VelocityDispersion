cmake_minimum_required(VERSION 3.13)

include_directories(${VelocityDispersion_SOURCE_DIR}/Testing)
link_directories(${VelocityDispersion_BINARY_DIR}/Testing)

FIND_PACKAGE(Boost COMPONENTS system filesystem unit_test_framework REQUIRED )

INCLUDE_DIRECTORIES( ${Boost_INCLUDE_DIR} )
set(CMAKE_CXX_STANDARD 14)


set ( Testing_SOURCES
        ../src/utillity.h
        ../src/utillity.cpp
        ../src/Integration.h
        ../src/Integration.cpp
        testmain.cpp
        testutil.cpp
        testutil.h
        testinteg.cpp
        testinteg.h)

set ( Executable_Name TestModule )

add_executable( ${Executable_Name} ${Testing_SOURCES})

TARGET_LINK_LIBRARIES(${Executable_Name} LINK_PUBLIC
        ${Boost_LIBRARIES}
        ${Boost_FILESYSTEM_LIBRARY}
        ${Boost_SYSTEM_LIBRARY}
        ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY}
        )

#target_link_libraries(${Executable_Name} ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})


enable_testing()

add_test(NAME VelocityDispersionTests COMMAND VelocityDispersionTests)
