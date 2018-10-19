set (OPTLIB_SOURCE ${OPTLIB_SOURCE}
	${CMAKE_CURRENT_SOURCE_DIR}/topopt.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/geneticoperators.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/topoptfactory.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/topoptgd.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/topoptuniverse.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/topoptchain.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/topoptga.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/topoptnlopt.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/topoptoc.cpp
	PARENT_SCOPE )
#add_subdirectory(unit_test)
