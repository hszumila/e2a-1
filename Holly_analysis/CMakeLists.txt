# CMakeLists for Holly's neutron analysis

include_directories(${CMAKE_SOURCE_DIR}/include ${ROOT_INCLUDE_DIR})

# Add executable for acceptance corrections
add_executable(plotter plotter.cpp) 
target_link_libraries(plotter ${ROOT_LIBRARIES} e2aCuts)

add_executable(plotMF plotMF.cpp)
target_link_libraries(plotMF ${ROOT_LIBRARIES} e2aCuts)

add_executable(plotter_3He plotter_3He.cpp)
target_link_libraries(plotter_3He ${ROOT_LIBRARIES} e2aCuts)

add_executable(plotSRC plotSRC.cpp)
target_link_libraries(plotSRC ${ROOT_LIBRARIES} e2aCuts)
