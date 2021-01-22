# CMakeLists for Holly's neutron analysis

include_directories(${CMAKE_SOURCE_DIR}/include ${ROOT_INCLUDE_DIR})

# Add executable for acceptance corrections
add_executable(plotMF plotMF.cpp)
target_link_libraries(plotMF ${ROOT_LIBRARIES} e2aCuts)

add_executable(plotSRC plotSRC.cpp)
target_link_libraries(plotSRC ${ROOT_LIBRARIES} e2aCuts)

add_executable(compareSRC compareSRC.cpp)
target_link_libraries(compareSRC ${ROOT_LIBRARIES} e2aCuts)

add_executable(compareMF compareMF.cpp)
target_link_libraries(compareMF ${ROOT_LIBRARIES} e2aCuts)

add_executable(calcMomRes calcMomRes.cpp)
target_link_libraries(calcMomRes ${ROOT_LIBRARIES} e2aCuts)

add_executable(calcMomRes2 calcMomRes2.cpp)
target_link_libraries(calcMomRes2 ${ROOT_LIBRARIES} e2aCuts)

add_executable(plot_neff plot_neff.cpp)
target_link_libraries(plot_neff ${ROOT_LIBRARIES} e2aCuts)

add_executable(calcNP_SRC calcNP_SRC.cpp)
target_link_libraries(calcNP_SRC ${ROOT_LIBRARIES} e2aCuts Acceptance)

add_executable(calcNP_MF calcNP_MF.cpp)
target_link_libraries(calcNP_MF ${ROOT_LIBRARIES} e2aCuts Acceptance)

add_executable(plot_neff_2261 plot_neff_2261.cpp)
target_link_libraries(plot_neff_2261 ${ROOT_LIBRARIES} e2aCuts)

