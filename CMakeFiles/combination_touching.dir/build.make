# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.9

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
CMAKE_COMMAND = /home/rasa/Documents/clion-2017.3.1/bin/cmake/bin/cmake

# The command to remove a file.
RM = /home/rasa/Documents/clion-2017.3.1/bin/cmake/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/rasa/NC-codes

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/rasa/NC-codes

# Include any dependencies generated for this target.
include CMakeFiles/combination_touching.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/combination_touching.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/combination_touching.dir/flags.make

CMakeFiles/combination_touching.dir/touching.cpp.o: CMakeFiles/combination_touching.dir/flags.make
CMakeFiles/combination_touching.dir/touching.cpp.o: touching.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/rasa/NC-codes/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/combination_touching.dir/touching.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/combination_touching.dir/touching.cpp.o -c /home/rasa/NC-codes/touching.cpp

CMakeFiles/combination_touching.dir/touching.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/combination_touching.dir/touching.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/rasa/NC-codes/touching.cpp > CMakeFiles/combination_touching.dir/touching.cpp.i

CMakeFiles/combination_touching.dir/touching.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/combination_touching.dir/touching.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/rasa/NC-codes/touching.cpp -o CMakeFiles/combination_touching.dir/touching.cpp.s

CMakeFiles/combination_touching.dir/touching.cpp.o.requires:

.PHONY : CMakeFiles/combination_touching.dir/touching.cpp.o.requires

CMakeFiles/combination_touching.dir/touching.cpp.o.provides: CMakeFiles/combination_touching.dir/touching.cpp.o.requires
	$(MAKE) -f CMakeFiles/combination_touching.dir/build.make CMakeFiles/combination_touching.dir/touching.cpp.o.provides.build
.PHONY : CMakeFiles/combination_touching.dir/touching.cpp.o.provides

CMakeFiles/combination_touching.dir/touching.cpp.o.provides.build: CMakeFiles/combination_touching.dir/touching.cpp.o


# Object files for target combination_touching
combination_touching_OBJECTS = \
"CMakeFiles/combination_touching.dir/touching.cpp.o"

# External object files for target combination_touching
combination_touching_EXTERNAL_OBJECTS =

combination_touching: CMakeFiles/combination_touching.dir/touching.cpp.o
combination_touching: CMakeFiles/combination_touching.dir/build.make
combination_touching: CMakeFiles/combination_touching.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/rasa/NC-codes/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable combination_touching"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/combination_touching.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/combination_touching.dir/build: combination_touching

.PHONY : CMakeFiles/combination_touching.dir/build

CMakeFiles/combination_touching.dir/requires: CMakeFiles/combination_touching.dir/touching.cpp.o.requires

.PHONY : CMakeFiles/combination_touching.dir/requires

CMakeFiles/combination_touching.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/combination_touching.dir/cmake_clean.cmake
.PHONY : CMakeFiles/combination_touching.dir/clean

CMakeFiles/combination_touching.dir/depend:
	cd /home/rasa/NC-codes && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/rasa/NC-codes /home/rasa/NC-codes /home/rasa/NC-codes /home/rasa/NC-codes /home/rasa/NC-codes/CMakeFiles/combination_touching.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/combination_touching.dir/depend
