% Sphere Example
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
%------------ BEGIN CODE ----------%
%% Set up 
clear
clc 

% Path to model file and to result file
modName = 'sphere'; % model name 
ADBSat_path = ADBSat_dynpath;
modIn = fullfile(ADBSat_path,'inou','obj_files',[modName,'.obj']);
modOut = fullfile(ADBSat_path,'inou','models');
resOut = fullfile(ADBSat_path,'inou','results',modName);

% Input orbital conditions
alt = 200; % orbital altitude [km]
inc = 51.6; % orbital inclination [deg]
% Environment conditions 
env = [alt*1e3, inc/2, 0, 106, 0, 65, 65, ones(1,7)*3, 0]; % Environment variables
% Attitude conditions 
aoa_deg = 10; % Angle of attack [deg]
aos_deg = 10; % Angle of sideslip [deg]

% Flag for shading algorithm
shadow = 1;

% Gas surface interaction model parameters
inparam.gsi_model = 'sentman';
inparam.alpha = 1; % Accommodation Coefficient [-]
inparam.Tw = 300; % Wall Temperature [K]

% Flag for solar radiation model
solar = 1;
% Solar Radiation Model parameters
inparam.sol_cR = 0.15; % Specular Reflectivity
inparam.sol_cD = 0.25; % Diffuse Reflectivity

% Flags for verbose mode and delete temp files
verb = 1;
del = 0;

%% Main Calculations 
% Import model
[modOut] = ADBSatImport(modIn, modOut, verb);

% Environment Calculations
inparam = environment(inparam, env(1),env(2),env(3),env(4),env(5),env(6),env(7),env(8:14),env(15));

% Coefficient Calculation
fileOut = calc_coeff(modOut, resOut, deg2rad(aoa_deg), deg2rad(aos_deg), inparam, shadow, solar, 1, 0); 

%% Post processing 
% Plot surface distribution
if verb && ~del
    plot_surfq(fileOut, modOut, aoa_deg(1), aos_deg(1), 'cp');
end
%------------ END CODE -----------%
