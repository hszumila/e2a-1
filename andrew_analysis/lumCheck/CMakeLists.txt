# CMakeLists for andrew's  analysis
include_directories(${CMAKE_SOURCE_DIR}/include ${ROOT_INCLUDE_DIR} ${andrew_include})

# Add executable for makeing cuts on the lead
add_executable(lumCheck lumCheck.cpp )
target_link_libraries(lumCheck ${ROOT_LIBRARIES} target_Info)

# Add executable for vertex
add_executable(vtxCheck vtxCheck.cpp)
target_link_libraries(vtxCheck ${ROOT_LIBRARIES} event_Info target_Info)
