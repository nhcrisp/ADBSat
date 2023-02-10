%EXAMPLE_ADB
% Creates an example aerodynamic database (in AoA/AoS) using ADBSat
%
% Author: Nicholas Crisp
% The University of Manchester
% November 2022
%
%--- Copyright notice ---%
% Copyright (C) 2022 The University of Manchester
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
%------------- BEGIN CODE --------------

clear

modName = 'cube';
% Path to model file
ADBSat_path = ADBSat_dynpath;
modIn = fullfile(ADBSat_path,'inou','obj_files',[modName,'.obj']);
modOut = fullfile(ADBSat_path,'inou','models');
resOut = fullfile(ADBSat_path,'inou','results');

%Input conditions
alt = 200; %km
inc = 51.6; %deg
env = [alt*1e3, inc/2, 0, 106, 0, 65, 65, ones(1,7)*3, 0]; % Environment variables

aoa = -45:1:45; % Angle of attack
aos = -45:1:45; % Angle of sideslip

% Model parameters
shadow = 1;
inparam.gsi_model = 'cook';
inparam.alpha = 1; % Accommodation (altitude dependent)
inparam.Tw = 300; % Wall Temperature [K]

solar = 1;
inparam.sol_cR = 0.15; % Specular Reflectivity
inparam.sol_cD = 0.25; % Diffuse Reflectivity

verb = 0; % Verbose
del = 1; % Delete temp

% Import model
modOut = ADBSatImport(modIn, modOut, verb);

% Coefficient Calculation
pathOut = ADBSatFcn( modOut, resOut, inparam, aoa, aos, shadow, solar, env, del, verb );

% Plot surface distribution
figure
load(pathOut,'aedb')
contourf((aedb.aos).*(180/pi), aedb.aoa.*(180/pi), aedb.aero.Cf_wX)
colorbar
xlabel('Angle of Sideslip')
ylabel('Angle of Attack')
title('Drag Force Coefficient')

%------------ END CODE -----------%

