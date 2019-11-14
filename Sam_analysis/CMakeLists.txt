

set(INCLUDES ${ROOT_INCLUDE_DIR} ${CMAKE_SOURCE_DIR}/include)
set(andrew_include ${CMAKE_SOURCE_DIR}/andrew_analysis/include ${CMAKE__SOURCE_DIR}/andrew_analysis/ratio_analysis)
include_directories(${INCLUDES} ${andrew_include})

# Add executable for event counting
add_executable(find_mass_combined find_mass_combined.cpp) 
target_link_libraries(find_mass_combined ${ROOT_LIBRARIES} target_Info)
