# Created by the script cgal_create_cmake_script_noqt
# This is the CMake script for compiling a CGAL application.


project( topologies )

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
# Compiler flags
if(DEFINE_USE_MPI)
	add_definitions(-DUSE_MPI)
endif(DEFINE_USE_MPI)
if(NOT DEFINED BOOST_ROOT)
	set(BOOST_ROOT $ENV{BOOST_ROOT})
endif(NOT DEFINED BOOST_ROOT)
# Directory definitions for libraries and headers
	set (NLOPT_INC_DIR ${NLOPT_ROOT}/include/)
	set (NLOPT_LIB_DIR ${NLOPT_ROOT}/lib/ )
	set (BOOST_LIB_DIR ${BOOST_ROOT}/lib)
	set (CATCH_INC_DIR ${CATCH_ROOT}/single_include)
	set (PUGIXML_INC_DIR ${PUGIXML_ROOT}/include)
	set (PUGIXML_LIB_DIR ${PUGIXML_ROOT}/lib64)
# Set include directories
	include( ${CGAL_USE_FILE} )
	include_directories (BEFORE ${CMAKE_CURRENT_BINARY_DIR})
	include_directories (BEFORE "IO/")
	include_directories (BEFORE "OBJFUN/")
	include_directories (BEFORE "OPT/")
	include_directories (BEFORE "REP/")
	include_directories (BEFORE "UTIL/")
	include_directories (BEFORE ${NLOPT_INC_DIR})
	include_directories (BEFORE ${PUGIXML_INC_DIR})
# Set source files
	add_subdirectory (IO)
	add_subdirectory (OBJFUN)
	add_subdirectory (OPT)
	add_subdirectory (REP)
	add_subdirectory (UTIL)
# Required libraries
	set (BOOST_REQ_LIBS ${BOOST_LIB_DIR}/libboost_thread.so ${BOOST_LIB_DIR}/libboost_chrono.so ${BOOST_LIB_DIR}/libboost_system.so)
	set (NLOPT_LIB ${NLOPT_LIB_DIR}/libnlopt.so)
	set (PUGIXML_LIB ${PUGIXML_LIB_DIR}/libpugixml.so)
	set (DL_LIB libdl.so)
# RPATH for install
	set (CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
# Set library
	link_directories(${NLOPT_LIB_DIR})
	link_directories(${BOOST_LIB_DIR})
	link_directories(${PUGIXML_LIB_DIR})
	add_library(topologiescore SHARED ${TOLIB_SOURCE})
	set_property(TARGET topologiescore PROPERTY CXX_STANDARD 11)
	set_property(TARGET topologiescore PROPERTY CXX_STANDARD_REQUIRED ON)
	set_property(TARGET topologiescore PROPERTY INSTALL_RPATH_USE_LINK_PATH TRUE)
	set(topologies_LIB topologiescore)
# Set executable
	add_subdirectory (MAIN)
	add_executable(topologies ${MAIN_SOURCE})
	set_property(TARGET topologies PROPERTY CXX_STANDARD 11)
	set_property(TARGET topologies PROPERTY CXX_STANDARD_REQUIRED ON)
	set_property(TARGET topologies PROPERTY INSTALL_RPATH_USE_LINK_PATH TRUE)
	add_to_cached_list( CGAL_EXECUTABLE_TARGETS topologies )
# Link the executable to CGAL and third-party libraries
	set (REQ_LIBS ${topologies_LIB} ${NLOPT_LIB} ${PUGIXML_LIB} ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES} ${BOOST_REQ_LIBS} ${DL_LIB})
	target_link_libraries(topologies ${REQ_LIBS})
# Set up unit and regression testing
if(COMPILE_TESTS)
	enable_testing()
	include_directories (BEFORE ${CATCH_INC_DIR})
	add_subdirectory (IO/unit_test)
	add_subdirectory (OBJFUN/unit_test)
	add_subdirectory (OPT/unit_test)
	add_subdirectory (REP/unit_test)
	add_subdirectory (UTIL/unit_test)
	# Note that the following tests may take a long time to complete:
	add_subdirectory (TEST)
endif(COMPILE_TESTS)
# Install
	install (TARGETS topologiescore LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX}/lib)
	install (TARGETS topologies RUNTIME DESTINATION ${CMAKE_INSTALL_PREFIX}/bin)
# Install example files
	add_subdirectory (EX_INPUT)
else()
  
    message(STATUS "NOTICE: The CGAL library is required and will not be compiled.")
  
endif()

