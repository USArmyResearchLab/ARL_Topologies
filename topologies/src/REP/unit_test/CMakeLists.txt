# Unit tests
# Needed for meshtestns.h
include_directories (BEFORE "REP/unit_test")
# CartesianMesherTest
add_executable(CartesianMesherTest test_cartesianmesher.cpp)
target_link_libraries(CartesianMesherTest ${REQ_LIBS})
add_test(CartesianMesherTest CartesianMesherTest)
# GeometryTranslationTest
add_executable(GeometryTranslationTest test_geometrytranslation.cpp)
target_link_libraries(GeometryTranslationTest ${REQ_LIBS})
add_test(GeometryTranslationTest GeometryTranslationTest)
# TOMeshTest
add_executable(TOMeshTest test_tomesh.cpp)
target_link_libraries(TOMeshTest ${REQ_LIBS})
add_test(TOMeshTest TOMeshTest)
# PixelRepTest
add_executable(PixelRepTest test_pixel.cpp)
target_link_libraries(PixelRepTest ${REQ_LIBS})
add_test(PixelRepTest PixelRepTest)
# VolMesh2DTest
add_executable(VolMesh2DTest test_volmesh2d.cpp)
target_link_libraries(VolMesh2DTest ${REQ_LIBS})
add_test(VolMesh2DTest VolMesh2DTest)
# VoxelRepTest
add_executable(VoxelRepTest test_voxel.cpp)
target_link_libraries(VoxelRepTest ${REQ_LIBS})
add_test(VoxelRepTest VoxelRepTest)
# VolMesh3DTest
add_executable(VolMesh3DTest test_volmesh3d.cpp)
target_link_libraries(VolMesh3DTest ${REQ_LIBS})
add_test(VolMesh3DTest VolMesh3DTest)
# CSGTest
add_executable(CSGTest test_csg.cpp)
target_link_libraries(CSGTest ${REQ_LIBS})
add_test(CSGTest CSGTest)
# Copy files
configure_file(./inputfiles/testpix.xml ${CMAKE_CURRENT_BINARY_DIR}/testpix.xml COPYONLY)
configure_file(./inputfiles/testheavi2.xml ${CMAKE_CURRENT_BINARY_DIR}/testheavi2.xml COPYONLY)
configure_file(./inputfiles/testmesh2d.xml ${CMAKE_CURRENT_BINARY_DIR}/testmesh2d.xml COPYONLY)
configure_file(./inputfiles/testmesh2d_fixed.xml ${CMAKE_CURRENT_BINARY_DIR}/testmesh2d_fixed.xml COPYONLY)
configure_file(./inputfiles/testmesh3d_fixed.xml ${CMAKE_CURRENT_BINARY_DIR}/testmesh3d_fixed.xml COPYONLY)
configure_file(./inputfiles/testquad.txt ${CMAKE_CURRENT_BINARY_DIR}/testquad.txt COPYONLY)
configure_file(./inputfiles/testheavimesh2d.xml ${CMAKE_CURRENT_BINARY_DIR}/testheavimesh2d.xml COPYONLY)
configure_file(./inputfiles/testvox.xml ${CMAKE_CURRENT_BINARY_DIR}/testvox.xml COPYONLY)
configure_file(./inputfiles/testheavi3.xml ${CMAKE_CURRENT_BINARY_DIR}/testheavi3.xml COPYONLY)
configure_file(./inputfiles/testmesh3d.xml ${CMAKE_CURRENT_BINARY_DIR}/testmesh3d.xml COPYONLY)
configure_file(./inputfiles/rectsolid.stl ${CMAKE_CURRENT_BINARY_DIR}/rectsolid.stl COPYONLY)
configure_file(./inputfiles/cylinder.txt ${CMAKE_CURRENT_BINARY_DIR}/cylinder.txt COPYONLY)
configure_file(./inputfiles/testheavimesh3d.xml ${CMAKE_CURRENT_BINARY_DIR}/testheavimesh3d.xml COPYONLY)
configure_file(./inputfiles/testcsg.xml ${CMAKE_CURRENT_BINARY_DIR}/testcsg.xml COPYONLY)
configure_file(./inputfiles/testcsg_point.xml ${CMAKE_CURRENT_BINARY_DIR}/testcsg_point.xml COPYONLY)
