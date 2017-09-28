set (TOLIB_SOURCE ${TOLIB_SOURCE} 
	${CMAKE_CURRENT_SOURCE_DIR}/exotxtmeshloader.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/gmshtxtloader.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/stltxtloader.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/inputloaderopt.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/inputloadertopopt.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/outputhandler.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/postprocess.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/inputloader.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/inputloaderrep.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/mpihandler.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/outputwriter.cpp
	PARENT_SCOPE )
#add_subdirectory(unit_test)
