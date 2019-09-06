# CMakeLists for Rey's neutron analysis

include_directories(${CMAKE_SOURCE_DIR}/include ${ROOT_INCLUDE_DIR})

# Add executable for acceptance corrections
add_executable(n_det_eff n_det_eff.cxx) 
target_link_libraries(n_det_eff ${ROOT_LIBRARIES})
