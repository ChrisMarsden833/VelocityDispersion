# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.14

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:


# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list


# Suppress display of executed commands.
$(VERBOSE).SILENT:


# A target that is always out of date.
cmake_force:

.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /Applications/CLion.app/Contents/bin/cmake/mac/bin/cmake

# The command to remove a file.
RM = /Applications/CLion.app/Contents/bin/cmake/mac/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /Users/christophermarsden/Documents/Astronomy/VelocityDispersion

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /Users/christophermarsden/Documents/Astronomy/VelocityDispersion/cmake-build-debug

# Include any dependencies generated for this target.
include CMakeFiles/VelocityDispersion.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/VelocityDispersion.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/VelocityDispersion.dir/flags.make

CMakeFiles/VelocityDispersion.dir/src/desmond.cpp.o: CMakeFiles/VelocityDispersion.dir/flags.make
CMakeFiles/VelocityDispersion.dir/src/desmond.cpp.o: ../src/desmond.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/christophermarsden/Documents/Astronomy/VelocityDispersion/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/VelocityDispersion.dir/src/desmond.cpp.o"
	/Library/Developer/CommandLineTools/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/VelocityDispersion.dir/src/desmond.cpp.o -c /Users/christophermarsden/Documents/Astronomy/VelocityDispersion/src/desmond.cpp

CMakeFiles/VelocityDispersion.dir/src/desmond.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/VelocityDispersion.dir/src/desmond.cpp.i"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/christophermarsden/Documents/Astronomy/VelocityDispersion/src/desmond.cpp > CMakeFiles/VelocityDispersion.dir/src/desmond.cpp.i

CMakeFiles/VelocityDispersion.dir/src/desmond.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/VelocityDispersion.dir/src/desmond.cpp.s"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/christophermarsden/Documents/Astronomy/VelocityDispersion/src/desmond.cpp -o CMakeFiles/VelocityDispersion.dir/src/desmond.cpp.s

CMakeFiles/VelocityDispersion.dir/src/integration.cpp.o: CMakeFiles/VelocityDispersion.dir/flags.make
CMakeFiles/VelocityDispersion.dir/src/integration.cpp.o: ../src/integration.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/christophermarsden/Documents/Astronomy/VelocityDispersion/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object CMakeFiles/VelocityDispersion.dir/src/integration.cpp.o"
	/Library/Developer/CommandLineTools/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/VelocityDispersion.dir/src/integration.cpp.o -c /Users/christophermarsden/Documents/Astronomy/VelocityDispersion/src/integration.cpp

CMakeFiles/VelocityDispersion.dir/src/integration.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/VelocityDispersion.dir/src/integration.cpp.i"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/christophermarsden/Documents/Astronomy/VelocityDispersion/src/integration.cpp > CMakeFiles/VelocityDispersion.dir/src/integration.cpp.i

CMakeFiles/VelocityDispersion.dir/src/integration.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/VelocityDispersion.dir/src/integration.cpp.s"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/christophermarsden/Documents/Astronomy/VelocityDispersion/src/integration.cpp -o CMakeFiles/VelocityDispersion.dir/src/integration.cpp.s

CMakeFiles/VelocityDispersion.dir/src/main.cpp.o: CMakeFiles/VelocityDispersion.dir/flags.make
CMakeFiles/VelocityDispersion.dir/src/main.cpp.o: ../src/main.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/christophermarsden/Documents/Astronomy/VelocityDispersion/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building CXX object CMakeFiles/VelocityDispersion.dir/src/main.cpp.o"
	/Library/Developer/CommandLineTools/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/VelocityDispersion.dir/src/main.cpp.o -c /Users/christophermarsden/Documents/Astronomy/VelocityDispersion/src/main.cpp

CMakeFiles/VelocityDispersion.dir/src/main.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/VelocityDispersion.dir/src/main.cpp.i"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/christophermarsden/Documents/Astronomy/VelocityDispersion/src/main.cpp > CMakeFiles/VelocityDispersion.dir/src/main.cpp.i

CMakeFiles/VelocityDispersion.dir/src/main.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/VelocityDispersion.dir/src/main.cpp.s"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/christophermarsden/Documents/Astronomy/VelocityDispersion/src/main.cpp -o CMakeFiles/VelocityDispersion.dir/src/main.cpp.s

CMakeFiles/VelocityDispersion.dir/src/utillity.cpp.o: CMakeFiles/VelocityDispersion.dir/flags.make
CMakeFiles/VelocityDispersion.dir/src/utillity.cpp.o: ../src/utillity.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/christophermarsden/Documents/Astronomy/VelocityDispersion/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Building CXX object CMakeFiles/VelocityDispersion.dir/src/utillity.cpp.o"
	/Library/Developer/CommandLineTools/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/VelocityDispersion.dir/src/utillity.cpp.o -c /Users/christophermarsden/Documents/Astronomy/VelocityDispersion/src/utillity.cpp

CMakeFiles/VelocityDispersion.dir/src/utillity.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/VelocityDispersion.dir/src/utillity.cpp.i"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/christophermarsden/Documents/Astronomy/VelocityDispersion/src/utillity.cpp > CMakeFiles/VelocityDispersion.dir/src/utillity.cpp.i

CMakeFiles/VelocityDispersion.dir/src/utillity.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/VelocityDispersion.dir/src/utillity.cpp.s"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/christophermarsden/Documents/Astronomy/VelocityDispersion/src/utillity.cpp -o CMakeFiles/VelocityDispersion.dir/src/utillity.cpp.s

# Object files for target VelocityDispersion
VelocityDispersion_OBJECTS = \
"CMakeFiles/VelocityDispersion.dir/src/desmond.cpp.o" \
"CMakeFiles/VelocityDispersion.dir/src/integration.cpp.o" \
"CMakeFiles/VelocityDispersion.dir/src/main.cpp.o" \
"CMakeFiles/VelocityDispersion.dir/src/utillity.cpp.o"

# External object files for target VelocityDispersion
VelocityDispersion_EXTERNAL_OBJECTS =

VelocityDispersion: CMakeFiles/VelocityDispersion.dir/src/desmond.cpp.o
VelocityDispersion: CMakeFiles/VelocityDispersion.dir/src/integration.cpp.o
VelocityDispersion: CMakeFiles/VelocityDispersion.dir/src/main.cpp.o
VelocityDispersion: CMakeFiles/VelocityDispersion.dir/src/utillity.cpp.o
VelocityDispersion: CMakeFiles/VelocityDispersion.dir/build.make
VelocityDispersion: CMakeFiles/VelocityDispersion.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/Users/christophermarsden/Documents/Astronomy/VelocityDispersion/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_5) "Linking CXX executable VelocityDispersion"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/VelocityDispersion.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/VelocityDispersion.dir/build: VelocityDispersion

.PHONY : CMakeFiles/VelocityDispersion.dir/build

CMakeFiles/VelocityDispersion.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/VelocityDispersion.dir/cmake_clean.cmake
.PHONY : CMakeFiles/VelocityDispersion.dir/clean

CMakeFiles/VelocityDispersion.dir/depend:
	cd /Users/christophermarsden/Documents/Astronomy/VelocityDispersion/cmake-build-debug && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/christophermarsden/Documents/Astronomy/VelocityDispersion /Users/christophermarsden/Documents/Astronomy/VelocityDispersion /Users/christophermarsden/Documents/Astronomy/VelocityDispersion/cmake-build-debug /Users/christophermarsden/Documents/Astronomy/VelocityDispersion/cmake-build-debug /Users/christophermarsden/Documents/Astronomy/VelocityDispersion/cmake-build-debug/CMakeFiles/VelocityDispersion.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/VelocityDispersion.dir/depend

