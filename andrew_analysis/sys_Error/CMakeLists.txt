# CMakeLists for andrew's  analysis

include_directories(${CMAKE_SOURCE_DIR}/include ${ROOT_INCLUDE_DIR} ${andrew_include})

# Add executable for calculating histograms
add_executable(get_xB_Dis get_xB_Dis.cpp)
target_link_libraries(get_xB_Dis ${ROOT_LIBRARIES} target_Info event_Info)
