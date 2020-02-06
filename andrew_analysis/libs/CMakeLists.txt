# CMakeLists for andrew's  analysis

include_directories(${CMAKE_SOURCE_DIR}/include ${ROOT_INCLUDE_DIR} ${andrew_include})

# Add library for target info
add_library(target_Info target_Info.cpp)
target_link_libraries(target_Info Acceptance e2aCuts ${ROOT_LIBRARIES})

# Add library for particle info
add_library(part_Info part_Info.cpp)
target_link_libraries(part_Info ${ROOT_LIBRARIES})

# Add library for event info
add_library(event_Info event_Info.cpp)
target_link_libraries(event_Info ${ROOT_LIBRARIES} part_Info)

# Add library for SRC Cut info
add_library(srcCut_Info srcCut_Info.cpp)
target_link_libraries(srcCut_Info ${ROOT_LIBRARIES})

