# Regression tests
# Test2D
add_executable(Test2D test2d.cpp)
target_link_libraries(Test2D ${REQ_LIBS})
add_test(Test2D Test2D)
# Test3D
add_executable(Test3D test3d.cpp)
target_link_libraries(Test3D ${REQ_LIBS})
add_test(Test3D Test3D)
# Copy and configure files
# Pixel files
configure_file(./inputfiles/test2d/testpix.xml ${CMAKE_CURRENT_BINARY_DIR}/testpix.xml)
configure_file(./inputfiles/test2d/tof_geo_2loads.xml ${CMAKE_CURRENT_BINARY_DIR}/tof_geo_2loads.xml COPYONLY)
configure_file(./inputfiles/test2d/pixres.txt ${CMAKE_CURRENT_BINARY_DIR}/pixres.txt COPYONLY)
# Pixel with tri mesh
configure_file(./inputfiles/test2d/testpixtri.xml ${CMAKE_CURRENT_BINARY_DIR}/testpixtri.xml)
configure_file(./inputfiles/test2d/tof_geo_1load.xml ${CMAKE_CURRENT_BINARY_DIR}/tof_geo_1load.xml COPYONLY)
configure_file(./inputfiles/test2d/pixrestri.txt ${CMAKE_CURRENT_BINARY_DIR}/pixrestri.txt COPYONLY)
# Heaviside pixel
configure_file(./inputfiles/test2d/testheavi2.xml ${CMAKE_CURRENT_BINARY_DIR}/testheavi2.xml)
configure_file(./inputfiles/test2d/heavi2res.txt ${CMAKE_CURRENT_BINARY_DIR}/heavi2res.txt COPYONLY)
# Mesh, 2 blocks
configure_file(./inputfiles/test2d/testmesh_2blocks.xml ${CMAKE_CURRENT_BINARY_DIR}/testmesh_2blocks.xml)
configure_file(./inputfiles/test2d/tof_exo.xml ${CMAKE_CURRENT_BINARY_DIR}/tof_exo.xml COPYONLY)
configure_file(./inputfiles/test2d/test2blocks_coarse.txt ${CMAKE_CURRENT_BINARY_DIR}/test2blocks_coarse.txt COPYONLY)
configure_file(./inputfiles/test2d/mesh2blockres.txt ${CMAKE_CURRENT_BINARY_DIR}/mesh2blockres.txt COPYONLY)
# Heaviside mesh, 2 blocks
configure_file(./inputfiles/test2d/testheavimesh2_2blocks.xml ${CMAKE_CURRENT_BINARY_DIR}/testheavimesh2_2blocks.xml)
configure_file(./inputfiles/test2d/heavimesh2blockres.txt ${CMAKE_CURRENT_BINARY_DIR}/heavimesh2blockres.txt COPYONLY)
# GMSH input, 2 blocks
configure_file(./inputfiles/test2d/testgmsh_2blocks.xml ${CMAKE_CURRENT_BINARY_DIR}/testgmsh_2blocks.xml)
configure_file(./inputfiles/test2d/tof_gmsh.xml ${CMAKE_CURRENT_BINARY_DIR}/tof_gmsh.xml COPYONLY)
configure_file(./inputfiles/test2d/test2d.msh ${CMAKE_CURRENT_BINARY_DIR}/test2d.msh COPYONLY)
configure_file(./inputfiles/test2d/gmsh2blockres.txt ${CMAKE_CURRENT_BINARY_DIR}/gmsh2blockres.txt COPYONLY)
# Element-density auto mesh, tris
configure_file(./inputfiles/test2d/testmesh2auto.xml ${CMAKE_CURRENT_BINARY_DIR}/testmesh2auto.xml)
configure_file(./inputfiles/test2d/mesh2autores.txt ${CMAKE_CURRENT_BINARY_DIR}/mesh2autores.txt COPYONLY)
# Voxel files
configure_file(./inputfiles/test3d/testvox.xml ${CMAKE_CURRENT_BINARY_DIR}/testvox.xml)
configure_file(./inputfiles/test3d/tof_geo_3d.xml ${CMAKE_CURRENT_BINARY_DIR}/tof_geo_3d.xml COPYONLY)
configure_file(./inputfiles/test3d/voxres.txt ${CMAKE_CURRENT_BINARY_DIR}/voxres.txt COPYONLY)
# Voxel with tetrahedral elemeents
configure_file(./inputfiles/test3d/testvoxtet.xml ${CMAKE_CURRENT_BINARY_DIR}/testvoxtet.xml)
configure_file(./inputfiles/test3d/voxrestet.txt ${CMAKE_CURRENT_BINARY_DIR}/voxrestet.txt COPYONLY)
# Heaviside voxel
configure_file(./inputfiles/test3d/testheavi3.xml ${CMAKE_CURRENT_BINARY_DIR}/testheavi3.xml)
configure_file(./inputfiles/test3d/heavi3res.txt ${CMAKE_CURRENT_BINARY_DIR}/heavi3res.txt COPYONLY)
# 3D tet mesh, 2 blocks
configure_file(./inputfiles/test3d/testmesh3_2blocks.xml ${CMAKE_CURRENT_BINARY_DIR}/testmesh3_2blocks.xml)
configure_file(./inputfiles/test3d/tof_exo_3d.xml ${CMAKE_CURRENT_BINARY_DIR}/tof_exo_3d.xml COPYONLY)
configure_file(./inputfiles/test3d/test2blocks3d.txt ${CMAKE_CURRENT_BINARY_DIR}/test2blocks3d.txt COPYONLY)
configure_file(./inputfiles/test3d/mesh3_2blocksres.txt ${CMAKE_CURRENT_BINARY_DIR}/mesh3_2blocksres.txt COPYONLY)
# 3d tet mesh, heaviside, 2 blocks
configure_file(./inputfiles/test3d/testheavimesh3_2blocks.xml ${CMAKE_CURRENT_BINARY_DIR}/testheavimesh3_2blocks.xml)
configure_file(./inputfiles/test3d/heavimesh3_2blocksres.txt ${CMAKE_CURRENT_BINARY_DIR}/heavimesh3_2blocksres.txt COPYONLY)
