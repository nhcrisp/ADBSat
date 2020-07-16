# ADBSat (Aerodynamic Database for Satellites)

ADBSat is a MATLAB toolkit for computing the aerodynamic database sets for different satellite geometries.

The main features of ADBSat are:
- Aerodynamic force and torque database generation from mesh geometry
- Solar radiation pressure force and torque database generation from mesh geometry
- Consideration of multiple surface/material characteristics on different faces 

## Installation
- Place unzipped ADBSat directory in desired location (usually MATLAB folder in user area)
- In MATLAB navigate to the install folder
- Run 'ADBSat_install.m' from MATLAB command line to add folders to MATLAB path

## Dependencies
- MATLAB
- Aerospace Toolbox
- meshlabserver (unless manually generated .obj files are available)

## Usage
- Requires Wavefront .obj mesh files as input geometry (generated using Blender or from .stl using meshlab)
- Ensure that meshlabserver is available on the system PATH
- MainImportFcn generates the geometry database (output in models directory)
- ADBSatMain generates the aerodynamic database (output in results directory)