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
include CMakeFiles/final_only_track.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/final_only_track.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/final_only_track.dir/flags.make

CMakeFiles/final_only_track.dir/final_only_track.cpp.o: CMakeFiles/final_only_track.dir/flags.make
CMakeFiles/final_only_track.dir/final_only_track.cpp.o: final_only_track.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/rasa/NC-codes/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/final_only_track.dir/final_only_track.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/final_only_track.dir/final_only_track.cpp.o -c /home/rasa/NC-codes/final_only_track.cpp

CMakeFiles/final_only_track.dir/final_only_track.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/final_only_track.dir/final_only_track.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/rasa/NC-codes/final_only_track.cpp > CMakeFiles/final_only_track.dir/final_only_track.cpp.i

CMakeFiles/final_only_track.dir/final_only_track.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/final_only_track.dir/final_only_track.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/rasa/NC-codes/final_only_track.cpp -o CMakeFiles/final_only_track.dir/final_only_track.cpp.s

CMakeFiles/final_only_track.dir/final_only_track.cpp.o.requires:

.PHONY : CMakeFiles/final_only_track.dir/final_only_track.cpp.o.requires

CMakeFiles/final_only_track.dir/final_only_track.cpp.o.provides: CMakeFiles/final_only_track.dir/final_only_track.cpp.o.requires
	$(MAKE) -f CMakeFiles/final_only_track.dir/build.make CMakeFiles/final_only_track.dir/final_only_track.cpp.o.provides.build
.PHONY : CMakeFiles/final_only_track.dir/final_only_track.cpp.o.provides

CMakeFiles/final_only_track.dir/final_only_track.cpp.o.provides.build: CMakeFiles/final_only_track.dir/final_only_track.cpp.o


# Object files for target final_only_track
final_only_track_OBJECTS = \
"CMakeFiles/final_only_track.dir/final_only_track.cpp.o"

# External object files for target final_only_track
final_only_track_EXTERNAL_OBJECTS =

final_only_track: CMakeFiles/final_only_track.dir/final_only_track.cpp.o
final_only_track: CMakeFiles/final_only_track.dir/build.make
final_only_track: /usr/lib/x86_64-linux-gnu/libboost_python.so
final_only_track: /usr/lib/libvtkGenericFiltering.so.5.10.1
final_only_track: /usr/lib/libvtkGeovis.so.5.10.1
final_only_track: /usr/lib/libvtkCharts.so.5.10.1
final_only_track: /usr/lib/libvtkViews.so.5.10.1
final_only_track: /usr/lib/libvtkInfovis.so.5.10.1
final_only_track: /usr/lib/libvtkWidgets.so.5.10.1
final_only_track: /usr/lib/libvtkVolumeRendering.so.5.10.1
final_only_track: /usr/lib/libvtkHybrid.so.5.10.1
final_only_track: /usr/lib/libvtkParallel.so.5.10.1
final_only_track: /usr/lib/libvtkRendering.so.5.10.1
final_only_track: /usr/lib/libvtkImaging.so.5.10.1
final_only_track: /usr/lib/libvtkGraphics.so.5.10.1
final_only_track: /usr/lib/libvtkIO.so.5.10.1
final_only_track: /usr/lib/libvtkFiltering.so.5.10.1
final_only_track: /usr/lib/libvtkCommon.so.5.10.1
final_only_track: /usr/lib/libvtksys.so.5.10.1
final_only_track: CMakeFiles/final_only_track.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/rasa/NC-codes/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable final_only_track"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/final_only_track.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/final_only_track.dir/build: final_only_track

.PHONY : CMakeFiles/final_only_track.dir/build

CMakeFiles/final_only_track.dir/requires: CMakeFiles/final_only_track.dir/final_only_track.cpp.o.requires

.PHONY : CMakeFiles/final_only_track.dir/requires

CMakeFiles/final_only_track.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/final_only_track.dir/cmake_clean.cmake
.PHONY : CMakeFiles/final_only_track.dir/clean

CMakeFiles/final_only_track.dir/depend:
	cd /home/rasa/NC-codes && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/rasa/NC-codes /home/rasa/NC-codes /home/rasa/NC-codes /home/rasa/NC-codes /home/rasa/NC-codes/CMakeFiles/final_only_track.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/final_only_track.dir/depend
