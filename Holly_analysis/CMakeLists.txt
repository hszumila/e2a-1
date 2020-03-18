# CMakeLists for Holly's neutron analysis

include_directories(${CMAKE_SOURCE_DIR}/include ${ROOT_INCLUDE_DIR})

# Add executable for acceptance corrections
add_executable(plotter plotter.cpp) 
target_link_libraries(plotter ${ROOT_LIBRARIES} e2aCuts)
