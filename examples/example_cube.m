% Cube Example
%
%--- Copyright notice ---%
% Copyright (C) 2021 The University of Manchester
% Written by David Mostaza Prieto,  Nicholas H. Crisp, Luciana Sinpetru and Sabrina Livadiotti
%
% This file is part of the ADBSat toolkit.
%
% This program is free software: you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or (at
% your option) any later version.
%
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
% Public License for more details.
%
% You should have received a copy of the GNU General Public License along
% with this program. If not, see <http://www.gnu.org/licenses/>.
%% ========================================================================
% CHANGELOG
% List of functions that have been updated:
% example_plate.m - now allows for easy selection of input file type and increased mesh processing control
% ADBSatImport.m - is able to handle all file types seamlessly
% importobjtri.m - validates, subdivides and "grades" mesh
% obj_fileTri2Patch.m - streamlined .obj file scanning to speed up process
%
% calc_coeff.m - now can account for hyper and hypo thermal flow conditions
% 
% List of functions that have been introduced:
% plotMeshQuality.m - shows .obj mesh quality plots: a complete heatmap and with problematic elements isolated
% stpToMesh.m - requires PDE Toolbox - converts .stp AP242 files into a mesh compatible with ADBSat
% plotMesh.m - shows the original and subdivided mesh for .obj files
% validateMesh.m - ensures that the obj mesh is watertight, orientable etc
%
% shadowAnalyHypo.m - performs path tracing to determine if panels are shadowed in hypothermal conditions
% 
% Other changes:
% Merged all input file locations into one folder
%% ========================================================================
%------------ BEGIN CODE ----------%
clear

modName = 'cube';
fileType = 'obj';
% Path to model file
ADBSat_path = ADBSat_dynpath;
% Change this line to suit the file location and type according to your
% needs. .stl does not contain material or appearance data, .obj contains 
% material data and .stp contains appearance data
modIn = fullfile(ADBSat_path,'inou','input_files',[modName,'.',fileType]);
modOut = fullfile(ADBSat_path,'inou','models');
resOut = fullfile(ADBSat_path,'inou','results',modName);

% Input conditions
alt = 200; %km
inc = 51.6; %deg
env = [alt*1e3, inc/2, 0, 106, 0, 65, 65, ones(1,7)*3, 0]; % Environment variables

aoa_deg = 0; % Angle of attack
aos_deg = 0; % Angle of sideslip

% Model parameters
shadow = 1;
inparam.gsi_model = 'sentman';

% Turning this off requires either the Sentman or Schaaf and Chambre GSIM.
% It will also DRASTICALLY increase computation time. You have been warned
inparam.hyperthermal = 0;
inparam.rays = 40;

% inparam.alpha represents the accomodation coefficients for each material.
% Ensure that it is ordered correctly using the debug plots provided
inparam.alpha = [1 1 1]; % Accommodation (altitude dependent)
inparam.Tw = 300; % Wall Temperature [K]

solar = 1;
inparam.sol_cR = 0.15; % Specular Reflectivity
inparam.sol_cD = 0.25; % Diffuse Reflectivity

% Mesh Parameters
meshParam.subdivisions = 0; % for .obj and .stl only
meshParam.elementSize = 0.002; % for .stp only
meshParam.translation = [0 0 0]; % X, Y, Z
meshParam.rotationAngles = [90 0 0]; % Pitch, Roll, Yaw
meshParam.centering = 1; % automatically centres the model around (0,0,0)

verb = 1;
del = 0;

% Import model
[modOut] = ADBSatImport(modIn, modOut, verb, meshParam);

% Environment Calculations
inparam = environment(inparam, env(1),env(2),env(3),env(4),env(5),env(6),env(7),env(8:14),env(15));

% Coefficient Calculation
fileOut = calc_coeff(modOut, resOut, deg2rad(aoa_deg), deg2rad(aos_deg), inparam, shadow, solar, 1, 0); 

% Plot surface distribution
if verb && ~del
    plot_surfq(fileOut, modOut, aoa_deg(1), aos_deg(1), 'cd');

    % Remove edges on plotted patches
    ax = gca;
    patches = findobj(ax, 'Type', 'Patch');
    for k = 1:numel(patches)
        set(patches(k), 'EdgeColor', 'none');
    end  
end

%------------ END CODE -----------%
