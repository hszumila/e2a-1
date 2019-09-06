# CMakeLists for creating maps and other programs
#   that require CLASTOOL


set (CLASTOOL ${CMAKE_SOURCE_DIR}/include/ClasTool)
set (ANALYZER ${CMAKE_SOURCE_DIR}/include/analysis_lib )
set (OSNAME Linux64RHEL7 )

include_directories(${CMAKE_SOURCE_DIR}/include ${CLASTOOL}/include ${ANALYZER}/include ${ROOT_INCLUDE_DIR} )
link_directories( ${ANALYZER}/slib ${CLASTOOL}/slib/${OSNAME} )

add_executable(convert_uniform convert_uniform.cxx)
target_link_libraries(convert_uniform ${ROOT_LIBRARIES})

add_executable(mapmaker mapmaker.cpp)
target_link_libraries(mapmaker ${ROOT_LIBRARIES} e2aCuts)

add_executable(mapmaker_eg2 mapmaker_eg2.cpp)
target_link_libraries(mapmaker_eg2 ${ROOT_LIBRARIES} e2aCuts)

add_executable(mapmaker_eg2_fine_bin mapmaker_eg2_fine_bin.cpp)
target_link_libraries(mapmaker_eg2_fine_bin ${ROOT_LIBRARIES} e2aCuts)

add_executable(generator_non_electrons generator_non_electrons.cxx)
target_link_libraries(generator_non_electrons ${ROOT_LIBRARIES})

add_executable(generator_electrons generator_electrons.cxx)
target_link_libraries(generator_electrons ${ROOT_LIBRARIES})
