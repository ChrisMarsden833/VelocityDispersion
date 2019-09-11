cmake_minimum_required(VERSION 3.13)
project(VelocityDispersion)

FIND_PACKAGE(Boost COMPONENTS system filesystem unit_test_framework REQUIRED )
INCLUDE_DIRECTORIES( ${Boost_INCLUDE_DIR} )

set(CMAKE_CXX_STANDARD 14)

#LINK_DIRECTORIES ( ${Boost_LIBRARY_DIRS} )

include_directories(src)

add_executable(VelocityDispersion
        src/desmond.cpp
        src/desmond.h
        src/integration.cpp
        src/integration.h
        src/main.cpp
        src/main.h
        src/utillity.cpp
        src/utillity.h)

add_executable(VelocityDispersionTests
        Testing/testmain.cpp
        src/utillity.h
        src/utillity.cpp
        src/Integration.h
        src/Integration.cpp
        Testing/testutil.cpp
        Testing/testutil.h
        Testing/testinteg.cpp
        Testing/testinteg.h)

TARGET_LINK_LIBRARIES(VelocityDispersion LINK_PUBLIC
        ${Boost_LIBRARIES}
        ${Boost_FILESYSTEM_LIBRARY}
        ${Boost_SYSTEM_LIBRARY}
       )

target_link_libraries(VelocityDispersionTests ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
enable_testing()
add_test(NAME VelocityDispersionTests COMMAND VelocityDispersionTests)






