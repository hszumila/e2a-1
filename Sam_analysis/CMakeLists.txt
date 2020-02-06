

set(INCLUDES ${ROOT_INCLUDE_DIR} ${CMAKE_SOURCE_DIR}/include)
set(andrew_include ${CMAKE_SOURCE_DIR}/andrew_analysis/include ${CMAKE_SOURCE_DIR}/andrew_analysis/ratio_analysis)
include_directories(${INCLUDES} ${andrew_include})

# Add executable for event counting
add_executable(find_mass_combined find_mass_combined.cpp) 
target_link_libraries(find_mass_combined ${ROOT_LIBRARIES} target_Info)

# Add executable for CLAS simulator
add_executable(simulator simulator.cpp )
target_link_libraries(simulator ${ROOT_LIBRARIES} target_Info Acceptance e2aCuts)

# Normalize SRC and Mean Field DATA
add_executable(scale_weight scale_weight.cpp)
target_link_libraries(scale_weight ${ROOT_LIBRARIES})
