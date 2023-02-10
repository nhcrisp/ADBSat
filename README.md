# ADBSat (Aerodynamic Database for Satellites)

ADBSat is a MATLAB toolkit for computing aerodynamic coefficient databases for satellite geometries in free-molecular flow (FMF) conditions.

The main features of ADBSat are:
- Import of body geometry from .stl or .obj mesh files
- Aerodynamic force and moment coefficient calculation for different gas-surface interaction models
- Solar radiation pressure force and moment coefficient calculation
- Consideration of multiple surface/material characteristics

ADBSat is a panel-method tool that is able to calculate aerodynamic or solar force and moment coefficient sets for satellite geometries by applying analytical (closed-form) expressions for the interactions to discrete flat-plate mesh elements. The panel method of ADBSat assumes FMF conditions.

ADBSat features a basic shadowing analysis to identify panels that are shielded from the flow by other parts of the body and will therefore not experience any surface interactions. However, this method is dependent on the refinement of the input mesh and can be sensitive to the orientation and arrangement of the mesh elements with respect to the oncoming flow direction.

ADBSat may therefore lose accuracy when:
- Parts of the geometry shadow others from the oncoming flow. 
- Secondary and multiple reflections between surface elements are expected to occur (thus it may not be suitable for the analysis of concave bodies).
- When the flow has a very low molecular speed ratio (hypothermal flow conditions) such that the oncoming flow may experience significant interactions with backwards facing surfaces.

## Installation
- Clone the git repository to the desired location (usually MATLAB folder in user area)
- In MATLAB navigate to the install folder
- Run 'ADBSat_install.m' from MATLAB command line to add the toolkit folders to MATLAB path

## Dependencies
- MATLAB 
- Aerospace Toolbox (for computation of environmental parameters if not supplied by the user)
- ~~meshlabserver (if automated conversion from .stl to .obj files is required)~~ **FUNCTIONALITY DEPRECIATED IN LATEST MESHLAB**

## Usage
- Requires Wavefront .obj mesh files as input geometry (generated using Blender or from .stl using meshlab)
- ~~Ensure that meshlabserver is available on the system PATH (if automated conversion from .stl to .obj files is required)~~
- ADBSatImport generates the geometric database file (output in models directory)
- ADBSatFcn generates the output force and moment coefficient database (output in results directory)

## Getting Started
The examples provided demonstrate the basic capabilities of ADBSat.
- example_import.m will 
- example_cube.m
- example_plate.m
- example_GSI.m will plot the variation in CD and CL with angle of incidence for two different GSI models (e.g. Sentman and CLL)

For most applications, the user will subsequently use the ADBSatImport.m and ADBSatFcn.m functions within their own scripts.

## Directory Structure
```bash
. 
├── examples                                    #Directory containing the example cube case.
|     ├── example_ADBSatFcn.m                   #Example run producing an aerodynamic database for the cube geometry for multiple AoA and AoS.
|     ├── example_ADBSatImport.m                #Import of the example cube object from a .obj to a .mat file.
|     ├── example_cube.m                        #Example import and aerodynamic analysis of the cube geometry at a single AoA and AoS.
|     ├── example_GSI.m                         #Plot example GSI models for surface inclination at varying accommodation coefficient, AoA/AoS.
|     └── example_plate.m                       #Example import and aerodynamic analysis of the plate geometry at a single AoA and AoS.
├── inou                                        #Directory containing the model inputs from the program.
|     ├── obj_files                             #Directory containing example .obj models (CAD inputs to ADBSat).
|     |     ├── cube.obj                        #CAD model of a cube used in the example cases.
|     |     └── plate.obj                       #CAD model of a plate used in the example cases.
|     └── stl_files                             #Directory containing example .stl models.
|           ├── cube.STL                        #CAD model of a cube used in the example cases (alternative).
├── install                                     #Directory containing the installation files
|     ├── ADBSat_dynpath.m                      #Creates the path to the ADBSat base folder
|     └── ADBSat_install.m                      #Installs ADBSat and adds necessary folders and paths to the functions.
├── toolbox                                     #Directory containing the program functions.
|     ├── calc                                  #Directory containing general calculation functions.
|     |     ├── ADBSatConstants.m               #Creates useful constants.
|     |     ├── calc_coeff.m                    #Calculates local and global aerodynamic and solar coefficients.
|     |     ├── environment.m                   #Interfaces with the MATLAB version of the NRLMSISE-00 atmospheric model.
|     |     ├── insidetri.m                     #Checks if a point is inside a triangle.
|     |     └── shadowAnaly.m                   #Checks if some mesh panels are shadowed by others.
|     ├── fmf_eq                                #Directory containing GSI model equations to calculate the drag coefficient.
|     |     ├── coeff_CLL.m                     #Cercignani-Lampis-Lord model.
|     |     ├── coeff_cook.m                    #Cook model.
|     |     ├── coeff_DRIA.m                    #Diffuse Reflection with Incomplete Accommodation model.
|     |     ├── coeff_maxwell.m                 #Maxwell model.
|     |     ├── coeff_newton.m                  #Newton hard-sphere model.
|     |     ├── coeff_schaaf.m                  #Schaaf and Chambre model.
|     |     ├── coeff_sentman.m                 #Sentman model.
|     |     ├── coeff_storchHyp.m               #Storch hyperthermal model.
|     |     └── mainCoeff.m                     #Passes the input parameters to the selected GSI model.
|     ├── import                                #Directory containing functions used in the model import.
|     |     ├── importobjtri.m                  #Imports a triangular mesh from a .obj file.
|     |     ├── meshlab_reset_origin.mlx        #Resets bounding box when converting CAD file from .stl to .obj.
|     |     ├── obj_fileTri2patch.m             #Reads a .obj file and outputs lists of vertices, faces, coordinates, and material identifiers. 
|     |     ├── stl2obj.m                       #Converts a .stl file to a .uobj file using meshlabserver (meshlab CLI).
|     |     └── surfaceNormals.m                #Calculates surface normals, areas, and barycentres of all elements of a triangular mesh.
|     ├── postpro                               #Directory containing post-processing functions.
|     |     ├── mergeAEDB.m                     #Merges coefficients calculated for different AOA/AOS into a single structure.
|     |     ├── plotNormals.m                   #Plots the result of the .obj import
|     |     └── plot_surfq.m                    #Plots the result of the ADBSat run, colour-coded to a chosen parameter
|     └── srp_eq                                #Directory containing solar radiation pressure model equations to calculate solar coefficients.
|           └── coeff_solar.m                   #Luthcke SRP model.
├── ADBSatFcn.m                                 #Function to run ADBSat.
├── ADBSatImport.m                              #Function to import a .obj CAD file and produce an output .mat file.
├── CHANGELOG.txt                               #Log of changes to the program.
├── LICENSE.txt                                 #License file.
└── README.md                                   #README file.
```

## References
Sinpetru, L. A., Crisp, N. H., Mostaza-Prieto, D., Livadiotti, S., and Roberts, P. C. E. “ADBSat: Methodology of a Novel Panel Method Tool for Aerodynamic Analysis of Satellites.” Computer Physics Communications, Vol. 275, 2022, p. 108326. <https://doi.org/10.1016/j.cpc.2022.108326>.

Sinpetru, L. A., Crisp, N. H., Roberts, P. C. E., Sulliotti-Linner, V., Hanessian, V., Herdrich, G. H., Romano, F., Garcia-Almiñana, D., Rodríguez-Donaire, S., and Seminari, S. “ADBSat: Verification and Validation of a Novel Panel Method for Quick Aerodynamic Analysis of Satellites.” Computer Physics Communications, Vol. 275, 2022, p. 108327. <https://doi.org/10.1016/j.cpc.2022.108327>.

Mostaza-Prieto, D., Characterisation and Applications of Aerodynamic Torques on Satellites (PhD Thesis), The University of Manchester, 2017.

## Licence
This project is licensed under the GPL-3.0 License - see the LICENSE file for details.