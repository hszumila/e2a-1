# CMakeLists for write_tree_e2a and other programs
#   that require CLASTOOL
include("${CMAKE_SOURCE_DIR}/cmake_modules/FindROOT.cmake")

include_directories(${ROOT_INCLUDE_DIRS})

set (CLASTOOL ${CMAKE_SOURCE_DIR}/include/ClasTool)
set (ANALYZER ${CMAKE_SOURCE_DIR}/include/analysis_lib )
set (OSNAME Linux64RHEL7 )

include_directories( ${CLASTOOL}/include ${ANALYZER}/include ${ROOT_INCLUDE_DIR} )
link_directories( ${ANALYZER}/slib ${CLASTOOL}/slib/${OSNAME} )

add_executable(write_tree_e2a write_tree_e2a.cxx)
target_link_libraries(write_tree_e2a ${ROOT_LIBRARIES} ClasTool ClasBanks TIdentificator)
