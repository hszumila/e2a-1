# CMakeLists for Adin's fiducial cuts

include_directories(${CMAKE_SOURCE_DIR}/include ${ROOT_INCLUDE_DIR})

# Add executable for acceptance corrections
add_executable(count count.cpp ) 
target_link_libraries(count ${ROOT_LIBRARIES} e2aCuts)

# Add executable for acceptance corrections
add_executable(fiducial fiducial.cpp)
target_link_libraries(fiducial ${ROOT_LIBRARIES})
