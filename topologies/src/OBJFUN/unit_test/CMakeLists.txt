# Unit tests
include_directories (BEFORE ${CMAKE_CURRENT_BINARY_DIR})
# TopOptObjFun
add_executable(TopOptObjFunTest test_topoptobjfun.cpp)
target_link_libraries(TopOptObjFunTest ${REQ_LIBS})
add_test(TopOptObjFunTest TopOptObjFunTest)
