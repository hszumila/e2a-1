# CMakeLists for andrew's  analysis

include_directories(${CMAKE_SOURCE_DIR}/include ${ROOT_INCLUDE_DIR} ${andrew_include})

# Add executable for calculating histograms for acceptance
add_executable(checkAcc checkAcc.cpp)
target_link_libraries(checkAcc ${ROOT_LIBRARIES} target_Info)

# Add executable for calculating values for acceptance
add_executable(queryLowAcc queryLowAcc.cpp)
target_link_libraries(queryLowAcc ${ROOT_LIBRARIES} target_Info)
