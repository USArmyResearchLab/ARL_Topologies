# Created by the script cgal_create_cmake_script_noqt
# This is the CMake script for compiling a CGAL application.


project( femofv )

cmake_minimum_required(VERSION 2.6.2)
if("${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}" VERSION_GREATER 2.6)
  if("${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}.${CMAKE_PATCH_VERSION}" VERSION_GREATER 2.8.3)
    cmake_policy(VERSION 2.8.4)
  else()
    cmake_policy(VERSION 2.6)
  endif()
endif()

add_definitions(-DCGAL_MESH_3_NO_DEPRECATED_SURFACE_INDEX
                -DCGAL_MESH_3_NO_DEPRECATED_C3T3_ITERATORS)
find_package(CGAL QUIET COMPONENTS Core )

if ( CGAL_FOUND )
# Definitions
if(DEFINE_USE_MPI)
	add_definitions(-DUSE_MPI)
endif(DEFINE_USE_MPI)
# Library locations
	set(TOPOLOGIES_SRC ${TOPOLOGIES_ROOT}/src)
	set(TOPOLOGIESCORE_LIB ${TOPOLOGIES_ROOT}/lib/libtopologiescore.so)
	set(NLOPT_LIB_DIR ${NLOPT_ROOT}/lib )
	set(NLOPT_LIB ${NLOPT_LIB_DIR}/libnlopt.so)
	set(CATCH_INC_DIR ${CATCH_ROOT}/single_include)
	set (PUGIXML_INC_DIR ${PUGIXML_ROOT}/include)
	set (PUGIXML_LIB_DIR ${PUGIXML_ROOT}/lib64)
# Set include directories
	include( ${CGAL_USE_FILE} )
	include_directories(BEFORE ${CMAKE_CURRENT_BINARY_DIR})
	include_directories (BEFORE ${TOPOLOGIES_SRC})
	include_directories (BEFORE ${TOPOLOGIES_SRC}/IO)
	include_directories (BEFORE ${TOPOLOGIES_SRC}/OBJFUN)
	include_directories (BEFORE ${TOPOLOGIES_SRC}/OPT)
	include_directories (BEFORE ${TOPOLOGIES_SRC}/REP)
	include_directories (BEFORE ${TOPOLOGIES_SRC}/UTIL)
	include_directories (BEFORE "UTIL/")
	include_directories (BEFORE "FEM/")
	include_directories (BEFORE ${EIGEN_ROOT})
	include_directories (BEFORE ${PUGIXML_INC_DIR})
# Set source files
	file (GLOB FEMSOURCE "FEM/*.cpp")
	file (GLOB UTILSOURCE "UTIL/*.cpp")
	set (FULL_SOURCE ${FEMSOURCE} ${UTILSOURCE})
# Libraries
	set (PUGIXML_LIB ${PUGIXML_LIB_DIR}/libpugixml.so)
	link_directories(${PUGIXML_LIB_DIR})
# Set executable
	add_library (femofv SHARED ${FULL_SOURCE})
	set_property(TARGET femofv PROPERTY INSTALL_RPATH_USE_LINK_PATH TRUE)
	set(FEM_LIB femofv)
	add_to_cached_list( CGAL_EXECUTABLE_TARGETS femofv )
# Link the executable to CGAL and third-party libraries
	target_link_libraries(femofv ${PUGIXML_LIB} ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES} )
# Set subdirectories, contain unit tests
	link_directories(${NLOPT_LIB_DIR})
	include_directories (BEFORE ${CATCH_INC_DIR})
if(COMPILE_TESTS)
	enable_testing()
	add_subdirectory (FEM)
	add_subdirectory (UTIL)
endif(COMPILE_TESTS)
# Install
	install (TARGETS femofv LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX}/lib)
else()
  
    message(STATUS "NOTICE: This demo requires the CGAL library, and will not be compiled.")
  
endif()

