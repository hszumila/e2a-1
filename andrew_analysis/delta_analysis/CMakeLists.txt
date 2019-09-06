# CMakeLists for andrew's  analysis
include_directories(${CMAKE_SOURCE_DIR}/include ${ROOT_INCLUDE_DIR} ${andrew_include})

# Add library for event info
add_library(event_Info event_Info.cpp)
target_link_libraries(event_Info ${ROOT_LIBRARIES})

# Add executable for makeing cuts on the lead
#add_executable(findDeltas findDeltas.cpp event_Info.cpp)
#target_link_libraries(findDeltas ${ROOT_LIBRARIES})
