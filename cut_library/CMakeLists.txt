# CMakeLists for the e2a analysis repository

SET(LIBNAME "e2aCuts")

include_directories( ${CMAKE_SOURCE_DIR}/include ${ROOT_INCLUDE_DIR} )

# Library for all of the cuts
add_library(${LIBNAME} Fiducial.cpp Run_dependent.cpp)
target_link_libraries(${LIBNAME} ${ROOT_LIBRARIES} )

add_executable(fiducial_test fiducial_test.cpp)
target_link_libraries(fiducial_test ${LIBNAME})

add_executable(fiducial_print fiducial_print.cpp)
target_link_libraries(fiducial_print ${LIBNAME})
