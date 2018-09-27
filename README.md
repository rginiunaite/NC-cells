# NC-cells

# About

This is the supplementary code repository for the paper:

Rebecca McLennan, Mary C. McKinney, Jessica M. Teddy, Jason A. Morrison, Jennifer C. Kasemeier-Kulesa, Dennis A. Ridenour, Craig Manthe, Rasa Giniunaite, Martin Robinson Ruth E. Baker, Philip K. Maini, Paul M. Kulesa, "Aquaporin-1 promotes cranial neural crest migration by stabilizing filopodia and influencing extracellular matrix degradation"

# Pre-requisites
Requires a C++14 compiler, Boost v1.65, Eigen v3, and CMake v2.8.

For example, these can be installed on Ubuntu 18.04 using apt:

$ sudo apt install build-essential libboost-all-dev cmake libeigen3-dev

# Installation
First clone this repository:

$ git clone --recurse-submodules https://github.com/rginiunaite/NC-cells.git
Then create a build directory build under the main source directory:

$ cd NC-cells
$ mkdir build
$ cd build
Then configure and compile the C++ module

$ cmake ..
$ make

