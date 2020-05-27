# CMakeLists for andrew's  analysis

set(andrew_include ${CMAKE_SOURCE_DIR}/andrew_analysis/include )

add_subdirectory(libs)
add_subdirectory(ratio_analysis)
add_subdirectory(momCorr_analysis)
add_subdirectory(lumCheck)
add_subdirectory(delta_analysis)
add_subdirectory(sys_Error)
add_subdirectory(acc_Check)
