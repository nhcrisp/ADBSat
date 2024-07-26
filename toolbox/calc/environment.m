% Computes a number of atmospheric parameters using the MATLAB atmosnrlmsise00 function.
% Assumes a circular orbit without atmospheric co-rotation or atmospheric winds.
%
% Inputs:
%   h               : Altitude [m]
%   lat             : Geodetic latitude [deg]
%   lon             : Longitude [deg]
%   dayOfYear       : Day of the year 1-365
%   UTseconds       : Seconds of the day
%   f107Average     : A 81 day average of F10.7 flux (centred on dayOfYear).
%   f107Daily       : Daily F10.7 flux for previous day.
%   magneticIndex   : A 1x7 array of magnetic index (AP) information
%   AnO             : Flag for anomalous oxygen
%
% Outputs:
%   vinf            : Free flow bulk velocity [m s^-1]
%   rho             : Atmospheric densities
%   s               : Speed ratio [-]
%   Rmean           : Specific gas constant [J kg^-1 K^-1]
%   Talt            : Temperature at altitude [K]
%
% Author: Nicholas H. Crisp
% The University of Manchester
% October 2018
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

function [ param_eq ] = environment( param_eq, h, lat, lon, dayOfYear, UTseconds, f107Average, f107Daily, magneticIndex, AnO )

% Constants
[data] = ADBSatConstants;

if AnO
    Oflag = 'Oxygen';
else
    Oflag = 'NoOxygen';
end

% Atmospheric properties
[T, param_eq.rho] = atmosnrlmsise00(h, lat, lon, 2002, dayOfYear, UTseconds, (UTseconds/3600 + lon/15), f107Average, f107Daily, magneticIndex, Oflag);

% Format temperature data
param_eq.Texo = T(1);
param_eq.Tinf = T(2);

% Calculate mean molecular mass [g mol^-1]
if AnO
    ND_T = (param_eq.rho(1)+param_eq.rho(2)+param_eq.rho(3)+param_eq.rho(4)+param_eq.rho(5)+param_eq.rho(7)+param_eq.rho(8)+param_eq.rho(9));
    param_eq.massConc(1,8) = param_eq.rho(8)/param_eq.rho(6) * (data.constants.mAnO/data.constants.NA/1000);

    param_eq.mmean = (data.constants.mHe * param_eq.rho(1)...
        + data.constants.mO * param_eq.rho(2)...
        + data.constants.mN2 * param_eq.rho(3)...
        + data.constants.mO2 * param_eq.rho(4)...
        + data.constants.mAr * param_eq.rho(5)...
        + data.constants.mH *  param_eq.rho(7)...
        + data.constants.mN * param_eq.rho(8)...
        + data.constants.mAnO * param_eq.rho(9))...
        /ND_T;
else
    ND_T = (param_eq.rho(1)+param_eq.rho(2)+param_eq.rho(3)+param_eq.rho(4)+param_eq.rho(5)+param_eq.rho(7)+param_eq.rho(8));

    param_eq.mmean = (data.constants.mHe * param_eq.rho(1)...
        + data.constants.mO * param_eq.rho(2)...
        + data.constants.mN2 * param_eq.rho(3)...
        + data.constants.mO2 * param_eq.rho(4)...
        + data.constants.mAr * param_eq.rho(5)...
        + data.constants.mH *  param_eq.rho(7)...
        + data.constants.mN * param_eq.rho(8))...
        /ND_T;
end

param_eq.massConc(1,1) = param_eq.rho(1)/param_eq.rho(6) * (data.constants.mHe/data.constants.NA/1000);
param_eq.massConc(1,2) = param_eq.rho(2)/param_eq.rho(6) * (data.constants.mO/data.constants.NA/1000);
param_eq.massConc(1,3) = param_eq.rho(3)/param_eq.rho(6) * (data.constants.mN2/data.constants.NA/1000);
param_eq.massConc(1,4) = param_eq.rho(4)/param_eq.rho(6) * (data.constants.mO2/data.constants.NA/1000);
param_eq.massConc(1,5) = param_eq.rho(5)/param_eq.rho(6) * (data.constants.mAr/data.constants.NA/1000);
param_eq.massConc(1,6) = param_eq.rho(7)/param_eq.rho(6) * (data.constants.mH/data.constants.NA/1000);
param_eq.massConc(1,7) = param_eq.rho(8)/param_eq.rho(6) * (data.constants.mN/data.constants.NA/1000);

% Calculate specific gas constant [J kg^-1 K^-1]
param_eq.Rmean = (data.constants.R/param_eq.mmean)*1000;

% Orbital velocity [m s^-1]
param_eq.vinf = sqrt(data.constants.mu_E/(data.constants.R_E+h));

% Thermal velocity
param_eq.vth = sqrt(2*data.constants.kb*param_eq.Tinf/(param_eq.mmean/data.constants.NA/1000));

%Speed ratio [-]
param_eq.s = param_eq.vinf/param_eq.vth;

%------------- END OF CODE --------------
