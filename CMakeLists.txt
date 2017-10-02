project( arl_topologies )

cmake_minimum_required(VERSION 2.6.2)
if("${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}" VERSION_GREATER 2.6)
  if("${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}.${CMAKE_PATCH_VERSION}" VERSION_GREATER 2.8.3)
    cmake_policy(VERSION 2.8.4)
  else()
    cmake_policy(VERSION 2.6)
  endif()
endif()

# Enable tests
if(COMPILE_TESTS)
  enable_testing()
endif(COMPILE_TESTS)

# Add topologies project and fem project
add_subdirectory(topologies/src)
add_subdirectory(fem_ofv/src)

