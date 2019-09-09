# CMakeLists for the e2a analysis repository

# Set up
cmake_minimum_required (VERSION 2.6)
set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake_modules")
find_package(ROOT REQUIRED)

# Project name
project (e2a)

# Create an install directory
if (DEFINED ENV{E2A_INSTALL})
   message(STATUS "Found E2A_INSTALL set to $ENV{E2A_INSTALL}")
else()
   message(FATAL_ERROR "E2A_INSTALL variable was not defined. It needs to be set. You could, for example, add to your .bashrc script:\n\
   	    export E2A_INSTALL=\${HOME}/.e2a")
endif()

# ADD SUBDIRECTORIES
add_subdirectory(cut_library)
add_subdirectory(calibration_data)
add_subdirectory(maps)
add_subdirectory(cut_production)
add_subdirectory(andrew_analysis)
add_subdirectory(Adin_analysis)
add_subdirectory(Axel_analysis)
add_subdirectory(Rey_neutrons)

# If on the farm, compile programs that require CLASTOOL
option(Build_Farm "Builds programs that require ifarm libraries" OFF )

if(Build_Farm)
	add_subdirectory(write_tree)
endif()

include_directories(${CMAKE_SOURCE_DIR}/include ${ROOT_INCLUDE_DIR})

# Add executable for making skim code
add_executable(skim_tree skim_tree.cpp) 
target_link_libraries(skim_tree ${ROOT_LIBRARIES} e2aCuts)
