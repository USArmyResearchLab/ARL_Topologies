enable_testing()
# Unit tests
# LinearTriangle
add_executable(LinTriTest test_lintri.cpp)
target_link_libraries(LinTriTest ${FEM_LIB} ${NLOPT_LIB} ${TOPOLOGIESCORE_LIB} )
add_test(LinTriTest LinTriTest)
# LinearQuadrilateral
add_executable(LinQuadTest test_linquad.cpp)
target_link_libraries(LinQuadTest ${FEM_LIB} ${NLOPT_LIB} ${TOPOLOGIESCORE_LIB} )
add_test(LinQuadTest LinQuadTest)
# LinTetra
add_executable(LinTetraTest test_lintetra.cpp)
target_link_libraries(LinTetraTest ${FEM_LIB} ${NLOPT_LIB} ${TOPOLOGIESCORE_LIB} )
add_test(LinTetraTest LinTetraTest)
# LinHex
add_executable(LinHexTest test_linhex.cpp)
target_link_libraries(LinHexTest ${FEM_LIB} ${NLOPT_LIB} ${TOPOLOGIESCORE_LIB} )
add_test(LinHexTest LinHexTest)
# Element, tests connectivity functions
add_executable(ElementTest test_element.cpp)
target_link_libraries(ElementTest ${FEM_LIB} ${NLOPT_LIB} ${TOPOLOGIESCORE_LIB} )
add_test(ElementTest ElementTest)
# ElemEdge, tests connectivity functions
add_executable(ElemEdgeTest test_elemedge.cpp)
target_link_libraries(ElemEdgeTest ${FEM_LIB} ${NLOPT_LIB} ${TOPOLOGIESCORE_LIB} )
add_test(ElemEdgeTest ElemEdgeTest)
# Mesh2D
add_executable(Mesh2DTest test_mesh2d.cpp)
target_link_libraries(Mesh2DTest ${FEM_LIB} ${NLOPT_LIB} ${TOPOLOGIESCORE_LIB} )
add_test(Mesh2DTest Mesh2DTest)
# Mesh3D
add_executable(Mesh3DTest test_mesh3d.cpp)
target_link_libraries(Mesh3DTest ${FEM_LIB} ${NLOPT_LIB} ${TOPOLOGIESCORE_LIB} )
add_test(Mesh3DTest Mesh3DTest)
# BoundaryCondition
add_executable(BoundaryConditionTest test_boundarycondition.cpp)
target_link_libraries(BoundaryConditionTest ${FEM_LIB} ${NLOPT_LIB} ${TOPOLOGIESCORE_LIB} )
add_test(BoundaryConditionTest BoundaryConditionTest)
# LoadCondition
add_executable(LoadConditionTest test_loadcondition.cpp)
target_link_libraries(LoadConditionTest ${FEM_LIB} ${NLOPT_LIB} ${TOPOLOGIESCORE_LIB} )
add_test(LoadConditionTest LoadConditionTest)
# FEMProblem
add_executable(FEMProblemTest test_femproblem.cpp)
target_link_libraries(FEMProblemTest ${FEM_LIB} ${NLOPT_LIB} ${TOPOLOGIESCORE_LIB} )
add_test(FEMProblemTest FEMProblemTest)
# TOFEMObjFun
add_executable(TOFEMObjFunTest test_tofemobjfun.cpp)
target_link_libraries(TOFEMObjFunTest ${FEM_LIB} ${NLOPT_LIB} ${TOPOLOGIESCORE_LIB} )
add_test(TOFEMObjFunTest TOFEMObjFunTest)
# copy input files
configure_file(./inputfiles/bctesttrimesh.txt ${CMAKE_CURRENT_BINARY_DIR}/bctesttrimesh.txt COPYONLY)
configure_file(./inputfiles/testtri.msh ${CMAKE_CURRENT_BINARY_DIR}/testtri.msh COPYONLY)
configure_file(./inputfiles/hexelemmat.txt ${CMAKE_CURRENT_BINARY_DIR}/hexelemmat.txt COPYONLY)
configure_file(./inputfiles/cylinder.txt ${CMAKE_CURRENT_BINARY_DIR}/cylinder.txt COPYONLY)
configure_file(./inputfiles/rectsolid.stl ${CMAKE_CURRENT_BINARY_DIR}/rectsolid.stl COPYONLY)
configure_file(./inputfiles/testheavimesh2d.xml ${CMAKE_CURRENT_BINARY_DIR}/testheavimesh2d.xml COPYONLY)
configure_file(./inputfiles/testheavimesh3d.xml ${CMAKE_CURRENT_BINARY_DIR}/testheavimesh3d.xml COPYONLY)
configure_file(./inputfiles/testmesh2d.xml ${CMAKE_CURRENT_BINARY_DIR}/testmesh2d.xml COPYONLY)
configure_file(./inputfiles/testmesh2dGMSH.xml ${CMAKE_CURRENT_BINARY_DIR}/testmesh2dGMSH.xml COPYONLY)
configure_file(./inputfiles/testmesh3d.xml ${CMAKE_CURRENT_BINARY_DIR}/testmesh3d.xml COPYONLY)
configure_file(./inputfiles/testquad.txt ${CMAKE_CURRENT_BINARY_DIR}/testquad.txt COPYONLY)
configure_file(./inputfiles/tof_exo_2.xml ${CMAKE_CURRENT_BINARY_DIR}/tof_exo_2.xml COPYONLY)
configure_file(./inputfiles/tof_gmsh_2.xml ${CMAKE_CURRENT_BINARY_DIR}/tof_gmsh_2.xml COPYONLY)
configure_file(./inputfiles/tof_exo_heavi2.xml ${CMAKE_CURRENT_BINARY_DIR}/tof_exo_heavi2.xml COPYONLY)
configure_file(./inputfiles/tof_exo_3.xml ${CMAKE_CURRENT_BINARY_DIR}/tof_exo_3.xml COPYONLY)
configure_file(./inputfiles/tof_exo_heavi3.xml ${CMAKE_CURRENT_BINARY_DIR}/tof_exo_heavi3.xml COPYONLY)
