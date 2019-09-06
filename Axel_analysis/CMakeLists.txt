# CMakeLists for various codes for Axel's analysis

include_directories(${CMAKE_SOURCE_DIR}/include ${ROOT_INCLUDE_DIR})

# Add executable for event counting
add_executable(example_analysis example_analysis.cpp) 
target_link_libraries(example_analysis ${ROOT_LIBRARIES})


# Add executable for program to look at empty target vertex positions
add_executable(study_empty study_empty.cpp)
target_link_libraries(study_empty ${ROOT_LIBRARIES})
