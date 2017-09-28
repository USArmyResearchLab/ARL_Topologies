# Unit tests
# TOMeshProcessing
add_executable(TOMeshProcessingTest test_tomeshprocessing.cpp)
target_link_libraries(TOMeshProcessingTest ${REQ_LIBS})
add_test(TOMeshProcessingTest TOMeshProcessingTest)
# Filter2D
add_executable(Filter2DTest test_filter2d.cpp)
target_link_libraries(Filter2DTest ${REQ_LIBS})
add_test(Filter2DTest Filter2DTest)
# Filter3D
add_executable(Filter3DTest test_filter3d.cpp)
target_link_libraries(Filter3DTest ${REQ_LIBS})
add_test(Filter3DTest Filter3DTest)
# Helper
add_executable(HelperTest test_helper.cpp)
target_link_libraries(HelperTest ${REQ_LIBS})
add_test(HelperTest HelperTest)
# GenericMaterial
add_executable(GenericMaterialTest test_genericmaterial.cpp)
target_link_libraries(GenericMaterialTest ${REQ_LIBS})
add_test(GenericMaterialTest GenericMaterialTest)
# copy input files
configure_file(./inputfiles/filter2dmesh.txt ${CMAKE_CURRENT_BINARY_DIR}/filter2dmesh.txt COPYONLY)
configure_file(./inputfiles/filter3dmesh.txt ${CMAKE_CURRENT_BINARY_DIR}/filter3dmesh.txt COPYONLY)
