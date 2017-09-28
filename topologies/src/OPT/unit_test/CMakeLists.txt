# TopOpt
add_executable(TopOptTest test_topopt.cpp)
target_link_libraries(TopOptTest ${REQ_LIBS})
add_test(TopOptTest TopOptTest)
# Copy files
configure_file(./inputfiles/gdinput.xml ${CMAKE_CURRENT_BINARY_DIR}/gdinput.xml COPYONLY)
configure_file(./inputfiles/ocinput.xml ${CMAKE_CURRENT_BINARY_DIR}/ocinput.xml COPYONLY)
configure_file(./inputfiles/bfgsinput.xml ${CMAKE_CURRENT_BINARY_DIR}/bfgsinput.xml COPYONLY)
configure_file(./inputfiles/mmainput.xml ${CMAKE_CURRENT_BINARY_DIR}/mmainput.xml COPYONLY)
configure_file(./inputfiles/gainput.xml ${CMAKE_CURRENT_BINARY_DIR}/gainput.xml COPYONLY)
configure_file(./inputfiles/chaininput.xml ${CMAKE_CURRENT_BINARY_DIR}/chaininput.xml COPYONLY)
