% Comparison of GSI Example
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
clc 
clear 

% Load the two gas-surface interaction model functions
GSI_model1 = @coeff_sentman;
GSI_model2 = @coeff_schaaf;

% Setting up the solver parameters manually
param_eq = struct;

% Set up the orbital and environmental variables
h = 300e3; % Altitude [m]
lat = 0;   % Geodetic Latitude [deg]
lon = 0;   % Longitude [deg]
dayOfYear = 1; % Day of the Year [1-365] 
UTseconds = 0; % Seconds of the day [s]
f107Average = 140; % Solar Activtiy 
f107Daily = 140;   
magneticIndex = ones([1,7])*15; % Geomagnetic Activity 
AnO = 1; % Flag for anomalous oxygen

% Incidence Angle to the flow (for a flat plate) [rad]
delta = deg2rad(0:1:90); 

%% Main Calculations
% Calculate the environmental variables
param_eq = environment(param_eq, h, lat, lon, dayOfYear, UTseconds, f107Average, f107Daily, magneticIndex, AnO);

% Set up of the Sentman GSI
param_eq.alpha = 1; % Accommodation Coefficient [-] 
param_eq.Tw = 300; % Wall temperature [K]
% Calculate the coefficients using Sentman GSI
[cp1, ctau1, cd1, cl1] = GSI_model1(param_eq, delta);

% Set up Schaaf and Chambre model
param_eq.sigmaN = 1; % Normal Accom Coefficient [-] 
param_eq.sigmaT = 1; % Tangential Accom Coefficient [-] 
% Calculate the coefficients using Schaaf and Chambre GSI
[cp2, ctau2, cd2, cl2] = GSI_model2(param_eq, delta);

%% Post Processing 
% Plotting the Cd and Cl distribution
hold on
p_1 = plot(delta*(180/pi),cd1);
set(p_1, 'Color', 'k', 'LineStyle','--', 'LineWidth',1.25)
p_2 = plot(delta*(180/pi),cd2);
set(p_2, 'Color', 'k', 'LineStyle','-','LineWidth',1.25)
grid on
ylabel('Drag Coefficient, C_D')

yyaxis right
set(gca, 'YColor',[0.5 0.5 0.5])
hold on
p_3 = plot(delta*(180/pi),cl1);
set(p_3, 'Color', [0.5 0.5 0.5], 'LineStyle','--','LineWidth',1.25)
p_4 = plot(delta*(180/pi),cl2);
set(p_4, 'Color', [0.5 0.5 0.5], 'LineStyle','-','LineWidth',1.25)
ylabel('Lift Coefficient, C_L')

xlabel('Incidence Angle [deg]')
legend('C_D Sentman (\alpha = 1)', 'C_D Schaaf and Chambre (\sigma_N = \sigma_T = 1)', 'C_L Sentman (\alpha = 1)', 'C_L Schaaf and Chambre (\sigma_N = \sigma_T = 1)')

set(gcf,'color','w');
