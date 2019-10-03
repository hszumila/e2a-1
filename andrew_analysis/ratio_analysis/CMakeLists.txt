# CMakeLists for andrew's  analysis

include_directories(${CMAKE_SOURCE_DIR}/include ${ROOT_INCLUDE_DIR} ${andrew_include})

# Add library for target info
add_library(target_Info target_Info.cpp)
target_link_libraries(target_Info Acceptance ${ROOT_LIBRARIES})

# Add library for event info
add_library(event_Info event_Info.cpp)
target_link_libraries(event_Info ${ROOT_LIBRARIES})

# Add executable for makeing cuts on the lead
add_executable(leadCut leadCut.cpp)
target_link_libraries(leadCut ${ROOT_LIBRARIES} event_Info)

# Add executable for calculating histograms
add_executable(get_hist get_hist.cpp)
target_link_libraries(get_hist ${ROOT_LIBRARIES} target_Info event_Info)

# Add executable for calculating histograms for acceptance
add_executable(checkAcc checkAcc.cpp)
target_link_libraries(checkAcc ${ROOT_LIBRARIES} target_Info)

# Add executable for calculating values for acceptance
add_executable(queryLowAcc queryLowAcc.cpp)
target_link_libraries(queryLowAcc ${ROOT_LIBRARIES} target_Info)

# Add executable for dividing histograms
add_executable(make_ratio make_ratio.cpp)
target_link_libraries(make_ratio ${ROOT_LIBRARIES})
