# ARL_Topologies

ARL_Topologies is an extensible topology optimization program meant to be both a research platform and usable in general.  To achieve this, it separates concepts from topology optimization into three parts: Representation, optimizer, and objective function.  

Currently, the following is implemented:

* **Representation**
	* 2D and 3D element-based density discretization with SIMP interpolation
	* 2D and 3D nodal-based density discretization using Heaviside projection and filtering
	* General mesh input for both of the above (ExodusII or GMSH, both in ASCII) or pixel/voxel rectangular regions
* **Optimizer**
	* Optimality criterion (with or without filtered sensitivities)
	* Method of moving asymptotes (via NLOpt, with or without filtered sensitivities)
	* BFGS (via NLOpt, with or without filtered sensitivities)
	* Genetic algorithm
* **Objective function**
	* Static, structural finite element

In addition, post-processing to STL is available for density-based representations in both 2D and 3D.  Output file formats include VTK, GMSH, and Matlab. 

## Prerequisites

The following libraries are required to compile ARL_Topologies:

* A C++11 compatible compiler
* [CGAL, version 4.6 or later](http://cgal.org)
  * [Boost, version 1.59 or later](http://www.boost.org/)
  * [GMP, version 6 or later](https://gmplib.org/)
  * [MPFR, version 3.1.3 or later](http://www.mpfr.org/)
* [PugiXML, version 1.7 or later](https://github.com/zeux/pugixml)
* [NLOpt, version 2.4.2 or later](http://ab-initio.mit.edu/wiki/index.php/NLopt)
* [Eigen, version 3.2.10 or later](http://eigen.tuxfamily.org) (required for the finite element solver)
* [Catch, version 1.5.8 or later](https://github.com/philsquared/Catch) (required for testing)

## Installation

1. Install all libraries listed above:
	1. CGAL requires Boost, GMP, and MPFR and recent versions are header-only
	2. PugiXML must be compiled as a shared library (set CMake option `BUILD_SHARED_LIBS` to `ON`)
	3. NLOpt must be compiled as a shared library
	4. Eigen is header-only and does not require compilation
	5. Catch is header-only
2. Build ARL_Topologies (ARL_Topologies uses the CMake build system and an out-of-source build is recommended)
	1. Create a build directory
	2. Edit the file `cmake_script_config`
		1. Change the paths to point toward the required libraries
		2. Change the last argument to point towards the `src` directory of ARL_Topologies
		3. Edit the compilers if necessary
		4. Edit the installation path
		5. By default, ARL_Topologies compiles against MPI, this can be removed by deleting `-DUSE_MPI` from the configuration command
	3. Run `cmake_script_config` in the build directory
	4. Run `make`
	5. Run `make install`
3. Build the FEM objective function shared library
	1. Follow the same steps as above, but in the objective function directory
	2. It is recommended that the same installation path is used for the FEM library and ARL_Topologies
	3. Run `make test` to execute the unit tests
4. Test the installation
	1. After installing the FEM library, if desired, run `make test` in the ARL_Topologies build directory
	* Note that all tests may take up to 5 minutes to complete

## Usage

ARL_Topologies takes only 1 or 2 inputs:
```
./topologies [-mpi] input.xml
```
where `-mpi` is to be used when MPI is enabled and useful (note that the included FEM library does not use MPI) and `input.xml` is the input file to be executed.  

### Input file format

ARL_Topologies (and the FEM library) uses an XML file format, a basic example of which is:
```
<?xml version="1.0" encoding="utf-8"?>
<topologies>
  <representation type="pixel" tag="rep1"/>
  <optimizer type="oc" tag="opt1"/>
  <objective_function shared_library="topologies/lib/libfemofv.so" input_file="tof_geo.xml"/>
  <initial_guess type="constant">
    <constant_val>0.5</constant_val>
  </initial_guess>
  <output type="volume">
    <file_name>mmvol</file_name>
    <file_format>vtk</file_format>
    <output_period>5</output_period>
    <overwrite>false</overwrite>
    <write_periodic_results>true</write_periodic_results>
  </output>
</topologies>
<!-- Representation definition, can also be defined in separate file -->
<pixel tag="rep1">
  <!-- Lists need an xml node as a delimiter, used <d/> below but could be anything -->
  <dimensions>1.<d/>1.</dimensions>
  <discretization_size>30<d/>30</discretization_size>
  <mesh_element_type>quad</mesh_element_type>
  <minimum_density>1.e-3</minimum_density>
  <penalty_power>3.</penalty_power>
  <threshold>0.5</threshold>
</pixel>
<!-- Optimizer definition, can also be defined in separate file -->
<oc tag="opt1">
  <filter_size>0.101</filter_size>
  <stop_tol>0.01</stop_tol>
  <max_iterations>30</max_iterations>
</oc>
```

with a basic FEM input given by

```
<?xml version="1.0" encoding="utf-8"?>
<fem>
  <dimension>2</dimension>
  <volume_fraction>0.5</volume_fraction>
  <youngs_modulus>1.</youngs_modulus>
  <poissons_ratio>0.25</poissons_ratio>
  <geo_bc>
    <x_support>true</x_support>
    <y_support>true</y_support>
    <geometry type="v_line">
      <intercept>0.</intercept>
    </geometry>
  </geo_bc>
  <load_case>
    <geo_lc>
      <load_vector>0.<d/>-1.</load_vector>
      <geometry type="point_2">
        <point_2 x="1." y="0."/>
      </geometry>
    </geo_lc>
  </load_case>
  <load_case>
    <geo_lc>
      <load_vector>0.<d/>1.</load_vector>
      <geometry type="point_2">
        <point_2 x="1." y="1."/>
      </geometry>
    </geo_lc>
  </load_case>
</fem>
```

Running `make install` will install several example input files in an `examples` directory. Please see the files contained in that directory for more information on the input file format

### Objective function shared library

Note that ARL_Topologies loads the objective function from a shared library, so that it is easy to extend to different problems.  The path to the shared library to load is specified in the `shared_library` attribute of the `objective_function` XML node.

## Contributing

This project is a work of the United States government and is not subject to domestic copyright protection under 17 USC § 105. Any outside contributions must be accompanied by a [Contributor License Agreement (CLA)](https://github.com/USArmyResearchLab/ARL-Open-Source-Guidance-and-Instructions#D3DC705AC3C411E6BBB4003EE1B763F8), which can be found [here](https://github.com/USArmyResearchLab/ARL-Open-Source-Guidance-and-Instructions/blob/master/ARL%20Form%20-%20266.pdf).

Please see [ARL's open source guidance](https://github.com/USArmyResearchLab/ARL-Open-Source-Guidance-and-Instructions) for more information.  

