################################################################################
include(DownloadProject)

# With CMake 3.8 and above, we can hide warnings about git being in a
# detached head by passing an extra GIT_CONFIG option
if(NOT (${CMAKE_VERSION} VERSION_LESS "3.8.0"))
	set(PE_EXTRA_OPTIONS "GIT_CONFIG advice.detachedHead=false")
else()
	set(PE_EXTRA_OPTIONS "")
endif()

# Shortcut function
function(pe_download_project name)
	download_project(
		PROJ         ${name}
		SOURCE_DIR   ${THIRD_PARTY_DIR}/${name}
		DOWNLOAD_DIR ${THIRD_PARTY_DIR}/.cache/${name}
		QUIET
		${PE_EXTRA_OPTIONS}
		${ARGN}
	)
endfunction()

################################################################################

## libigl
function(pe_download_libigl)
	pe_download_project(libigl
		GIT_REPOSITORY https://github.com/libigl/libigl.git
		GIT_TAG        682e4b9685d2737215f6629ecafcb318d714d556
	)
endfunction()

## tbb
function(pe_download_tbb)
    pe_download_project(tbb
    GIT_REPOSITORY https://github.com/wjakob/tbb.git
    GIT_TAG        20357d83871e4cb93b2c724fe0c337cd999fd14f
  )
endfunction()
