# CMakeLists for e2a maps

include(downloader)

# Set project name
PROJECT(Acceptance)

add_subdirectory(creation)
add_subdirectory(diagnostics)

include_directories( ${CMAKE_SOURCE_DIR}/include ${ROOT_INCLUDE_DIR} )

add_library(Acceptance Acceptance.cpp )
target_link_libraries(Acceptance ${ROOT_LIBRARIES} )

#################################################
#
# Download and install the maps
#
#################################################

download (e2a_He3_1161_750_e
http://mit.edu/~schmidta/www/e2a_maps/e2a_He3_1161_750_e.root maps/e2a_He3_1161_750_e.root)
download (e2a_He3_1161_750_p
http://mit.edu/~schmidta/www/e2a_maps/e2a_He3_1161_750_p.root maps/e2a_He3_1161_750_p.root)
download (e2a_He3_1161_750_pip
http://mit.edu/~schmidta/www/e2a_maps/e2a_He3_1161_750_pip.root maps/e2a_He3_1161_750_pip.root)
download (e2a_He3_1161_1500_e
http://mit.edu/~schmidta/www/e2a_maps/e2a_He3_1161_1500_e.root maps/e2a_He3_1161_1500_e.root)
download (e2a_He3_1161_1500_p
http://mit.edu/~schmidta/www/e2a_maps/e2a_He3_1161_1500_p.root maps/e2a_He3_1161_1500_p.root)
download (e2a_He3_1161_1500_pip
http://mit.edu/~schmidta/www/e2a_maps/e2a_He3_1161_1500_pip.root maps/e2a_He3_1161_1500_pip.root)
download (e2a_He3_2261_2250_e
http://mit.edu/~schmidta/www/e2a_maps/e2a_He3_2261_2250_e.root maps/e2a_He3_2261_2250_e.root)
download (e2a_He3_2261_2250_p
http://mit.edu/~schmidta/www/e2a_maps/e2a_He3_2261_2250_p.root maps/e2a_He3_2261_2250_p.root)
download (e2a_He3_2261_2250_pip
http://mit.edu/~schmidta/www/e2a_maps/e2a_He3_2261_2250_pip.root maps/e2a_He3_2261_2250_pip.root)
download (e2a_He3_4461_2250_e
http://mit.edu/~schmidta/www/e2a_maps/e2a_He3_4461_2250_e.root maps/e2a_He3_4461_2250_e.root)
download (e2a_He3_4461_2250_p
http://mit.edu/~schmidta/www/e2a_maps/e2a_He3_4461_2250_p.root maps/e2a_He3_4461_2250_p.root)
download (e2a_He3_4461_2250_pip
http://mit.edu/~schmidta/www/e2a_maps/e2a_He3_4461_2250_pip.root maps/e2a_He3_4461_2250_pip.root)
download (e2a_He4_2261_2250_e
http://mit.edu/~schmidta/www/e2a_maps/e2a_He4_2261_2250_e.root maps/e2a_He4_2261_2250_e.root)
download (e2a_He4_2261_2250_p
http://mit.edu/~schmidta/www/e2a_maps/e2a_He4_2261_2250_p.root maps/e2a_He4_2261_2250_p.root)
download (e2a_He4_2261_2250_pip
http://mit.edu/~schmidta/www/e2a_maps/e2a_He4_2261_2250_pip.root maps/e2a_He4_2261_2250_pip.root)
download (e2a_He4_4461_2250_e
http://mit.edu/~schmidta/www/e2a_maps/e2a_He4_4461_2250_e.root maps/e2a_He4_4461_2250_e.root)
download (e2a_He4_4461_2250_p
http://mit.edu/~schmidta/www/e2a_maps/e2a_He4_4461_2250_p.root maps/e2a_He4_4461_2250_p.root)
download (e2a_He4_4461_2250_pip
http://mit.edu/~schmidta/www/e2a_maps/e2a_He4_4461_2250_pip.root maps/e2a_He4_4461_2250_pip.root)
download (e2a_solid_1161_750_e
http://mit.edu/~schmidta/www/e2a_maps/e2a_solid_1161_750_e.root maps/e2a_solid_1161_750_e.root)
download (e2a_solid_1161_750_p
http://mit.edu/~schmidta/www/e2a_maps/e2a_solid_1161_750_p.root maps/e2a_solid_1161_750_p.root)
download (e2a_solid_1161_750_pip
http://mit.edu/~schmidta/www/e2a_maps/e2a_solid_1161_750_pip.root maps/e2a_solid_1161_750_pip.root)
download (e2a_solid_1161_1500_e
http://mit.edu/~schmidta/www/e2a_maps/e2a_solid_1161_1500_e.root maps/e2a_solid_1161_1500_e.root)
download (e2a_solid_1161_1500_p
http://mit.edu/~schmidta/www/e2a_maps/e2a_solid_1161_1500_p.root maps/e2a_solid_1161_1500_p.root)
download (e2a_solid_1161_1500_pip
http://mit.edu/~schmidta/www/e2a_maps/e2a_solid_1161_1500_pip.root maps/e2a_solid_1161_1500_pip.root)
download (e2a_solid_2261_2250_e
http://mit.edu/~schmidta/www/e2a_maps/e2a_solid_2261_2250_e.root maps/e2a_solid_2261_2250_e.root)
download (e2a_solid_2261_2250_p
http://mit.edu/~schmidta/www/e2a_maps/e2a_solid_2261_2250_p.root maps/e2a_solid_2261_2250_p.root)
download (e2a_solid_2261_2250_pip
http://mit.edu/~schmidta/www/e2a_maps/e2a_solid_2261_2250_pip.root maps/e2a_solid_2261_2250_pip.root)
download (e2a_solid_4461_2250_e
http://mit.edu/~schmidta/www/e2a_maps/e2a_solid_4461_2250_e.root maps/e2a_solid_4461_2250_e.root)
download (e2a_solid_4461_2250_p
http://mit.edu/~schmidta/www/e2a_maps/e2a_solid_4461_2250_p.root maps/e2a_solid_4461_2250_p.root)
download (e2a_solid_4461_2250_pip
http://mit.edu/~schmidta/www/e2a_maps/e2a_solid_4461_2250_pip.root maps/e2a_solid_4461_2250_pip.root)

add_dependencies(Acceptance e2a_He3_1161_750_e e2a_He3_1161_750_p e2a_He3_1161_750_pip 
			    e2a_He3_1161_1500_e e2a_He3_1161_1500_p e2a_He3_1161_1500_pip
			    e2a_He3_2261_2250_e e2a_He3_2261_2250_p e2a_He3_2261_2250_pip
			    e2a_He3_4461_2250_e e2a_He3_4461_2250_p e2a_He3_4461_2250_pip
			    e2a_He4_2261_2250_e e2a_He4_2261_2250_p e2a_He4_2261_2250_pip
			    e2a_He4_4461_2250_e e2a_He4_4461_2250_p e2a_He4_4461_2250_pip
			    e2a_solid_1161_750_e e2a_solid_1161_750_p e2a_solid_1161_750_pip 
			    e2a_solid_1161_1500_e e2a_solid_1161_1500_p e2a_solid_1161_1500_pip
			    e2a_solid_2261_2250_e e2a_solid_2261_2250_p e2a_solid_2261_2250_pip
			    e2a_solid_4461_2250_e e2a_solid_4461_2250_p e2a_solid_4461_2250_pip )

# Install the completed maps so other programs can find them
INSTALL ( DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/maps/" 
	  DESTINATION $ENV{HOME}/.e2a/maps
	  FILES_MATCHING PATTERN "e2a_*.root"
	  PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ GROUP_EXECUTE GROUP_READ )

