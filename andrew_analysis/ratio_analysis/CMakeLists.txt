# CMakeLists for andrew's  analysis

include_directories(${CMAKE_SOURCE_DIR}/include ${ROOT_INCLUDE_DIR} ${andrew_include})

# Add executable for makeing cuts on the lead
add_executable(leadCut leadCut.cpp)
target_link_libraries(leadCut ${ROOT_LIBRARIES} event_Info srcCut_Info)

# Add executable for calculating histograms
add_executable(get_hist get_hist.cpp)
target_link_libraries(get_hist ${ROOT_LIBRARIES} target_Info event_Info)

# Add executable for dividing histograms
add_executable(make_ratio make_ratio.cpp)
target_link_libraries(make_ratio ${ROOT_LIBRARIES})
