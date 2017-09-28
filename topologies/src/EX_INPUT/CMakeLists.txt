# Example input files
# Copy and configure files (sets absolute paths for libraries)
# 2D Mesh, separate input files for rep and opt
configure_file(./inpfiles/quad_mesh_hole.txt ${CMAKE_CURRENT_BINARY_DIR}/inpfiles/quad_mesh_hole.txt COPYONLY)
configure_file(./inpfiles/tof_geo.xml ${CMAKE_CURRENT_BINARY_DIR}/inpfiles/tof_geo.xml)
configure_file(./inpfiles/tooinp.xml ${CMAKE_CURRENT_BINARY_DIR}/inpfiles/tooinp.xml)
configure_file(./inpfiles/topopt.xml ${CMAKE_CURRENT_BINARY_DIR}/inpfiles/topopt.xml)
configure_file(./inpfiles/torinp.xml ${CMAKE_CURRENT_BINARY_DIR}/inpfiles/torinp.xml)
# 2D Pixel rep
configure_file(./inpfiles_pix/topopt.xml ${CMAKE_CURRENT_BINARY_DIR}/inpfiles_pix/topopt.xml)
configure_file(./inpfiles_pix/tof_geo.xml ${CMAKE_CURRENT_BINARY_DIR}/inpfiles_pix/tof_geo.xml)
# 3D Voxel rep
configure_file(./inpfiles_vox/topopt.xml ${CMAKE_CURRENT_BINARY_DIR}/inpfiles_vox/topopt.xml)
configure_file(./inpfiles_vox/tof_geo.xml ${CMAKE_CURRENT_BINARY_DIR}/inpfiles_vox/tof_geo.xml)
configure_file(./inpfiles_vox/topopt_pp.xml ${CMAKE_CURRENT_BINARY_DIR}/inpfiles_vox/topopt_pp.xml)
# Fixed block examples
configure_file(./inpfiles_2blocks/test2blocks.txt ${CMAKE_CURRENT_BINARY_DIR}/inpfiles_2blocks/test2blocks.txt COPYONLY)
configure_file(./inpfiles_2blocks/test2blocks3d.txt ${CMAKE_CURRENT_BINARY_DIR}/inpfiles_2blocks/test2blocks3d.txt COPYONLY)
configure_file(./inpfiles_2blocks/tof_exo.xml ${CMAKE_CURRENT_BINARY_DIR}/inpfiles_2blocks/tof_exo.xml)
configure_file(./inpfiles_2blocks/tof_exo_3d.xml ${CMAKE_CURRENT_BINARY_DIR}/inpfiles_2blocks/tof_exo_3d.xml)
configure_file(./inpfiles_2blocks/tof_geo_3d.xml ${CMAKE_CURRENT_BINARY_DIR}/inpfiles_2blocks/tof_geo_3d.xml)
configure_file(./inpfiles_2blocks/topopt.xml ${CMAKE_CURRENT_BINARY_DIR}/inpfiles_2blocks/topopt.xml)
configure_file(./inpfiles_2blocks/topopt_3d.xml ${CMAKE_CURRENT_BINARY_DIR}/inpfiles_2blocks/topopt_3d.xml)
# GMSH example
configure_file(./inpfiles_gmsh/trimesh.msh ${CMAKE_CURRENT_BINARY_DIR}/inpfiles_gmsh/trimesh.msh COPYONLY)
configure_file(./inpfiles_gmsh/tof_exo.xml ${CMAKE_CURRENT_BINARY_DIR}/inpfiles_gmsh/tof_exo.xml)
configure_file(./inpfiles_gmsh/topopt.xml ${CMAKE_CURRENT_BINARY_DIR}/inpfiles_gmsh/topopt.xml)
# Heaviside method
configure_file(./inpfiles_heavi2/tof_geo.xml ${CMAKE_CURRENT_BINARY_DIR}/inpfiles_heavi2/tof_geo.xml)
configure_file(./inpfiles_heavi2/topopt.xml ${CMAKE_CURRENT_BINARY_DIR}/inpfiles_heavi2/topopt.xml)
# Install to examples directory
# 2D Mesh
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/inpfiles/quad_mesh_hole.txt 
	${CMAKE_CURRENT_BINARY_DIR}/inpfiles/tof_geo.xml
	${CMAKE_CURRENT_BINARY_DIR}/inpfiles/tooinp.xml
	${CMAKE_CURRENT_BINARY_DIR}/inpfiles/topopt.xml
	${CMAKE_CURRENT_BINARY_DIR}/inpfiles/torinp.xml
	DESTINATION ${CMAKE_INSTALL_PREFIX}/examples/inpfiles)
# 2D Pixel
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/inpfiles_pix/topopt.xml  
	${CMAKE_CURRENT_BINARY_DIR}/inpfiles_pix/tof_geo.xml
	DESTINATION ${CMAKE_INSTALL_PREFIX}/examples/inpfiles_pix)
# 3D Voxel
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/inpfiles_vox/topopt.xml  
	${CMAKE_CURRENT_BINARY_DIR}/inpfiles_vox/tof_geo.xml
	${CMAKE_CURRENT_BINARY_DIR}/inpfiles_vox/topopt_pp.xml
	DESTINATION ${CMAKE_INSTALL_PREFIX}/examples/inpfiles_vox)
# Fixed block
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/inpfiles_2blocks/test2blocks.txt
	${CMAKE_CURRENT_BINARY_DIR}/inpfiles_2blocks/test2blocks3d.txt
	${CMAKE_CURRENT_BINARY_DIR}/inpfiles_2blocks/tof_exo.xml
	${CMAKE_CURRENT_BINARY_DIR}/inpfiles_2blocks/tof_exo_3d.xml
	${CMAKE_CURRENT_BINARY_DIR}/inpfiles_2blocks/tof_geo_3d.xml
	${CMAKE_CURRENT_BINARY_DIR}/inpfiles_2blocks/topopt.xml
	${CMAKE_CURRENT_BINARY_DIR}/inpfiles_2blocks/topopt_3d.xml
	DESTINATION ${CMAKE_INSTALL_PREFIX}/examples/inpfiles_2blocks)
# GMSH
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/inpfiles_gmsh/trimesh.msh
	${CMAKE_CURRENT_BINARY_DIR}/inpfiles_gmsh/tof_exo.xml
	${CMAKE_CURRENT_BINARY_DIR}/inpfiles_gmsh/topopt.xml
	DESTINATION ${CMAKE_INSTALL_PREFIX}/examples/inpfiles_gmsh)
# Heaviside
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/inpfiles_heavi2/tof_geo.xml
	${CMAKE_CURRENT_BINARY_DIR}/inpfiles_heavi2/topopt.xml
	DESTINATION ${CMAKE_INSTALL_PREFIX}/examples/inpfiles_heavi2)

