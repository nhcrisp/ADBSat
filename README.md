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
- Secondary and multiple reflections between surface elements are expected to occur (this it may not be suitable for the analysis of concave bodies).
- When the flow has a very low molecular speed ratio (hypothermal flow conditions) such that the oncoming flow may experience significant interactions with backwards facing surfaces.

## Installation
- Clone the git repository to the desired location (usually MATLAB folder in user area)
- In MATLAB navigate to the install folder
- Run 'ADBSat_install.m' from MATLAB command line to add the toolkit folders to MATLAB path

## Dependencies
- MATLAB 
- Aerospace Toolbox
- meshlabserver (if automated conversion from .stl to .obj files is required)

## Usage
- Requires Wavefront .obj mesh files as input geometry (generated using Blender or from .stl using meshlab)
- Ensure that meshlabserver is available on the system PATH
- ADBSatImport generates the geometric database file (output in models directory)
- ADBSatFcn generates the output force and moment coefficient database (output in results directory)

## References

D. Mostaza-Prieto, Characterisation and Applications of Aerodynamic Torques on Satellites (PhD Thesis), The University of Manchester, 2017.