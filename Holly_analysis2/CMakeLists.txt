# CMakeLists for Holly's neutron analysis

set(holly_include ${CMAKE_SOURCE_DIR}/Holly_analysis2/include)

include_directories(${CMAKE_SOURCE_DIR}/include ${ROOT_INCLUDE_DIR} ${holly_include})

add_executable(calcNCorr calcNCorr.cpp)
target_link_libraries(calcNCorr ${ROOT_LIBRARIES} e2aCuts)

add_executable(plotNeff plotNeff.cpp)
target_link_libraries(plotNeff ${ROOT_LIBRARIES} e2aCuts)

add_executable(findSRCcuts findSRCcuts.cpp)
target_link_libraries(findSRCcuts ${ROOT_LIBRARIES} e2aCuts)

add_executable(findMFcuts findMFcuts.cpp)
target_link_libraries(findMFcuts ${ROOT_LIBRARIES} e2aCuts)

add_executable(ratioNP_SRC ratioNP_SRC.cpp)
target_link_libraries(ratioNP_SRC ${ROOT_LIBRARIES} e2aCuts Acceptance)

add_executable(ratioNP_MF ratioNP_MF.cpp)
target_link_libraries(ratioNP_MF ${ROOT_LIBRARIES} e2aCuts Acceptance)

