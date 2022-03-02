% Creates a structure of useful constants.
%
% Author: Nicholas H. Crisp
% The University of Manchester
% October 2019
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
%------------- BEGIN CODE -------------

function [data] = ADBSatConstants(data)

data.constants.mu_E = 3.986004418e14;  % GM Earth [m^3/s^2] (source: SMAD 3rd edition)
data.constants.R_E = 6.37813649e6; % Equatorial Earth radius [m] (source: SMAD 3rd edition)

data.constants.R = 8.31446261815324; % Molar (Universal/Ideal) Gas constant [J K^-1 mol^-1]
data.constants.kb = 1.3806503e-23; % Boltzmann constant [m^2 kg^-2 K^-1]
data.constants.NA = 6.02214076e23; % Avogadro constant [n mol^-1]

data.constants.mHe = 2; % Molecular mass of He [g mol^-1]
data.constants.mO = 15.9994; % Molecular mass of O [g mol^-1]
data.constants.mN2 = 28.0134; % Molecular mass of N2 [g mol^-1]
data.constants.mO2 = 31.9988; % Molecular mass of O2 [g mol^-1]
data.constants.mAr = 39.948; % Molecular mass of Ar [g mol^-1]
data.constants.mH = 1.0079; % Molecular mass of H [g mol^-1]
data.constants.mN = 14.0067; % Molecular mass of N [g mol^-1]
data.constants.mAnO = 15.99; % Molecular mass of anomalous oxygen [g mol^-1]

%------------- END OF CODE --------------
