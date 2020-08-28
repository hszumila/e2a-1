# CMakeLists for various codes for Erin's analysis

include_directories(${CMAKE_SOURCE_DIR}/include ${ROOT_INCLUDE_DIR})

# Add executable for event counting
add_executable(helloworld helloworld.cpp) 
target_link_libraries(helloworld ${ROOT_LIBRARIES})


