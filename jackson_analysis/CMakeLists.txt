# CMakeLists for jackson's  analysis

include_directories(${CMAKE_SOURCE_DIR}/include ${ROOT_INCLUDE_DIR} ${andrew_include})

# Add executable for makeing cuts on the lead
add_executable(apply_SRC_Cut apply_SRC_Cut.cpp)
target_link_libraries(apply_SRC_Cut ${ROOT_LIBRARIES} event_Info srcCut_Info)

# Add executable for calculating histograms
add_executable(get_mMiss_Hist get_mMiss_Hist.cpp)
target_link_libraries(get_mMiss_Hist ${ROOT_LIBRARIES} target_Info event_Info)
