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
include CMakeFiles/check.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/check.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/check.dir/flags.make

CMakeFiles/check.dir/open_space_improved_version.cpp.o: CMakeFiles/check.dir/flags.make
CMakeFiles/check.dir/open_space_improved_version.cpp.o: open_space_improved_version.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/rasa/NC-codes/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/check.dir/open_space_improved_version.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/check.dir/open_space_improved_version.cpp.o -c /home/rasa/NC-codes/open_space_improved_version.cpp

CMakeFiles/check.dir/open_space_improved_version.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/check.dir/open_space_improved_version.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/rasa/NC-codes/open_space_improved_version.cpp > CMakeFiles/check.dir/open_space_improved_version.cpp.i

CMakeFiles/check.dir/open_space_improved_version.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/check.dir/open_space_improved_version.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/rasa/NC-codes/open_space_improved_version.cpp -o CMakeFiles/check.dir/open_space_improved_version.cpp.s

CMakeFiles/check.dir/open_space_improved_version.cpp.o.requires:

.PHONY : CMakeFiles/check.dir/open_space_improved_version.cpp.o.requires

CMakeFiles/check.dir/open_space_improved_version.cpp.o.provides: CMakeFiles/check.dir/open_space_improved_version.cpp.o.requires
	$(MAKE) -f CMakeFiles/check.dir/build.make CMakeFiles/check.dir/open_space_improved_version.cpp.o.provides.build
.PHONY : CMakeFiles/check.dir/open_space_improved_version.cpp.o.provides

CMakeFiles/check.dir/open_space_improved_version.cpp.o.provides.build: CMakeFiles/check.dir/open_space_improved_version.cpp.o


# Object files for target check
check_OBJECTS = \
"CMakeFiles/check.dir/open_space_improved_version.cpp.o"

# External object files for target check
check_EXTERNAL_OBJECTS =

check: CMakeFiles/check.dir/open_space_improved_version.cpp.o
check: CMakeFiles/check.dir/build.make
check: /usr/lib/x86_64-linux-gnu/libboost_python.so
check: /usr/lib/libvtkGenericFiltering.so.5.10.1
check: /usr/lib/libvtkGeovis.so.5.10.1
check: /usr/lib/libvtkCharts.so.5.10.1
check: /usr/lib/libvtkViews.so.5.10.1
check: /usr/lib/libvtkInfovis.so.5.10.1
check: /usr/lib/libvtkWidgets.so.5.10.1
check: /usr/lib/libvtkVolumeRendering.so.5.10.1
check: /usr/lib/libvtkHybrid.so.5.10.1
check: /usr/lib/libvtkParallel.so.5.10.1
check: /usr/lib/libvtkRendering.so.5.10.1
check: /usr/lib/libvtkImaging.so.5.10.1
check: /usr/lib/libvtkGraphics.so.5.10.1
check: /usr/lib/libvtkIO.so.5.10.1
check: /usr/lib/libvtkFiltering.so.5.10.1
check: /usr/lib/libvtkCommon.so.5.10.1
check: /usr/lib/libvtksys.so.5.10.1
check: CMakeFiles/check.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/rasa/NC-codes/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable check"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/check.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/check.dir/build: check

.PHONY : CMakeFiles/check.dir/build

CMakeFiles/check.dir/requires: CMakeFiles/check.dir/open_space_improved_version.cpp.o.requires

.PHONY : CMakeFiles/check.dir/requires

CMakeFiles/check.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/check.dir/cmake_clean.cmake
.PHONY : CMakeFiles/check.dir/clean

CMakeFiles/check.dir/depend:
	cd /home/rasa/NC-codes && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/rasa/NC-codes /home/rasa/NC-codes /home/rasa/NC-codes /home/rasa/NC-codes /home/rasa/NC-codes/CMakeFiles/check.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/check.dir/depend
