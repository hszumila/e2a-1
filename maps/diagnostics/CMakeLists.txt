# CMakeLists for diagnostic tools for testing maps

include_directories(${CMAKE_SOURCE_DIR}/include ${ROOT_INCLUDE_DIR})

add_executable(basic_test basic_test.cpp)
target_link_libraries(basic_test ${ROOT_LIBRARIES})

add_executable(tof_test tof_test.cpp)
target_link_libraries(tof_test ${ROOT_LIBRARIES} e2aCuts)

add_executable(map_query map_query.cpp)
target_link_libraries(map_query ${ROOT_LIBRARIES} Acceptance)
