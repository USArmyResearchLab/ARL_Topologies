# ExoTxtMeshLoader
add_executable(ExoTxtMeshLoaderTest test_exotxtmeshloader.cpp)
target_link_libraries(ExoTxtMeshLoaderTest ${REQ_LIBS})
add_test(ExoTxtMeshLoaderTest ExoTxtMeshLoaderTest)
# GMSHTxtMeshLoader
add_executable(GMSHTxtMeshLoaderTest test_gmshtxtloader.cpp)
target_link_libraries(GMSHTxtMeshLoaderTest ${REQ_LIBS})
add_test(GMSHTxtMeshLoaderTest GMSHTxtMeshLoaderTest)
# STLTxtLoader
add_executable(STLTxtLoaderTest test_stltxtloader.cpp)
target_link_libraries(STLTxtLoaderTest ${REQ_LIBS})
add_test(STLTxtLoaderTest STLTxtLoaderTest)
# PostProcess
add_executable(PostProcessTest test_postprocess.cpp)
target_link_libraries(PostProcessTest ${REQ_LIBS})
add_test(PostProcessTest PostProcessTest)
# TORGenericVolume
add_executable(TORGenericVolumeTest test_torgenericvolume.cpp)
target_link_libraries(TORGenericVolumeTest ${REQ_LIBS})
add_test(TORGenericVolumeTest TORGenericVolumeTest)
# TORGenericMesh
add_executable(TORGenericMeshTest test_torgenericmesh.cpp)
target_link_libraries(TORGenericMeshTest ${REQ_LIBS})
add_test(TORGenericMeshTest TORGenericMeshTest)
# TORCSG
add_executable(TORCSGTest test_torcsg.cpp)
target_link_libraries(TORCSGTest ${REQ_LIBS})
add_test(TORCSGTest TORCSGTest)
# TOOGeneric
add_executable(TOOGenericTest test_toogeneric.cpp)
target_link_libraries(TOOGenericTest ${REQ_LIBS})
add_test(TOOGenericTest TOOGenericTest)
# TOOGA
add_executable(TOOGATest test_tooga.cpp)
target_link_libraries(TOOGATest ${REQ_LIBS})
add_test(TOOGATest TOOGATest)
# TOOChain
add_executable(TOOChainTest test_toochain.cpp)
target_link_libraries(TOOChainTest ${REQ_LIBS})
add_test(TOOChainTest TOOChainTest)
# NodeInfo
add_executable(NodeInfoTest test_nodeinfo.cpp)
target_link_libraries(NodeInfoTest ${REQ_LIBS})
add_test(NodeInfoTest NodeInfoTest)
# PugiXMLWrapper
add_executable(PugiXMLWrapperTest test_pugixmlwrapper.cpp)
target_link_libraries(PugiXMLWrapperTest ${REQ_LIBS})
add_test(PugiXMLWrapperTest PugiXMLWrapperTest)
# MainInput
add_executable(MainInputTest test_maininput.cpp)
target_link_libraries(MainInputTest ${REQ_LIBS})
add_test(MainInputTest MainInputTest)
# copy input files
configure_file(./inputfiles/testtri.txt ${CMAKE_CURRENT_BINARY_DIR}/testtri.txt COPYONLY)
configure_file(./inputfiles/testtri.msh ${CMAKE_CURRENT_BINARY_DIR}/testtri.msh COPYONLY)
configure_file(./inputfiles/testquad.txt ${CMAKE_CURRENT_BINARY_DIR}/testquad.txt COPYONLY)
configure_file(./inputfiles/testtet.txt ${CMAKE_CURRENT_BINARY_DIR}/testtet.txt COPYONLY)
configure_file(./inputfiles/testtet.msh ${CMAKE_CURRENT_BINARY_DIR}/testtet.msh COPYONLY)
configure_file(./inputfiles/testhex.txt ${CMAKE_CURRENT_BINARY_DIR}/testhex.txt COPYONLY)
configure_file(./inputfiles/test2blocks.txt ${CMAKE_CURRENT_BINARY_DIR}/test2blocks.txt COPYONLY)
configure_file(./inputfiles/test2blocks.msh ${CMAKE_CURRENT_BINARY_DIR}/test2blocks.msh COPYONLY)
configure_file(./inputfiles/testpix.xml ${CMAKE_CURRENT_BINARY_DIR}/testpix.xml COPYONLY)
configure_file(./inputfiles/testmesh2d.xml ${CMAKE_CURRENT_BINARY_DIR}/testmesh2d.xml COPYONLY)
configure_file(./inputfiles/testpugixml.xml ${CMAKE_CURRENT_BINARY_DIR}/testpugixml.xml COPYONLY)
configure_file(./inputfiles/nodeinfo.xml ${CMAKE_CURRENT_BINARY_DIR}/nodeinfo.xml COPYONLY)
configure_file(./inputfiles/testcsg.xml ${CMAKE_CURRENT_BINARY_DIR}/testcsg.xml COPYONLY)
configure_file(./inputfiles/testtoo.xml ${CMAKE_CURRENT_BINARY_DIR}/testtoo.xml COPYONLY)
configure_file(./inputfiles/testga.xml ${CMAKE_CURRENT_BINARY_DIR}/testga.xml COPYONLY)
configure_file(./inputfiles/testchain.xml ${CMAKE_CURRENT_BINARY_DIR}/testchain.xml COPYONLY)
configure_file(./inputfiles/testtopopt.xml ${CMAKE_CURRENT_BINARY_DIR}/testtopopt.xml COPYONLY)
configure_file(./inputfiles/testoutput.xml ${CMAKE_CURRENT_BINARY_DIR}/testoutput.xml COPYONLY)
configure_file(./inputfiles/testinitialguess.xml ${CMAKE_CURRENT_BINARY_DIR}/testinitialguess.xml COPYONLY)
configure_file(./inputfiles/cube.stl ${CMAKE_CURRENT_BINARY_DIR}/cube.stl COPYONLY)
