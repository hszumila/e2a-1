# CMakeLists for auxiliary scripts

include_directories(${CMAKE_SOURCE_DIR}/include ${ROOT_INCLUDE_DIR})

# Add executable for acceptance corrections
add_executable(sim_pid_electron sim_pid_electron.cpp )
target_link_libraries(sim_pid_electron ${ROOT_LIBRARIES} e2aCuts)
