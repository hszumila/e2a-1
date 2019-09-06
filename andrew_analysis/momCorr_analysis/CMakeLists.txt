# CMakeLists for andrew's  analysis
include_directories(${CMAKE_SOURCE_DIR}/include ${ROOT_INCLUDE_DIR})

# Add executable for makeing cuts on the lead
add_executable(nMass nMass.cpp)
target_link_libraries(nMass ${ROOT_LIBRARIES})
