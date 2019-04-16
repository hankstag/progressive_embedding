# Umfpack lib usually requires linking to a blas library.
# It is up to the user of this module to find a BLAS and link to it.

if (UMFPACK_INCLUDES AND UMFPACK_LIBRARIES)
  set(UMFPACK_FIND_QUIETLY TRUE)
endif (UMFPACK_INCLUDES AND UMFPACK_LIBRARIES)

find_path(UMFPACK_INCLUDES
  NAMES
  umfpack.h
  PATHS
  $ENV{UMFPACKDIR}
  ${INCLUDE_INSTALL_DIR}
  PATH_SUFFIXES
  Include
)
message(STATUS "ihave it aalllllll ${UMFPACK_INCLUDES}")

find_library(UMFPACK_LIBRARIES umfpack
	HINTS
	${UMFPACK_INCLUDES}/../Lib
	PATHS 
	$ENV{UMFPACKDIR} 
	${LIB_INSTALL_DIR})
message(STATUS "ifdsfadsfasd ${UMFPACK_LIBRARIES}")

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(UMFPACK DEFAULT_MSG
                                  UMFPACK_INCLUDES UMFPACK_LIBRARIES)

mark_as_advanced(UMFPACK_INCLUDES UMFPACK_LIBRARIES AMD_LIBRARY COLAMD_LIBRARY)
add_library(umfpack INTERFACE)
target_link_libraries(umfpack INTERFACE ${UMFPACK_LIBRARIES})