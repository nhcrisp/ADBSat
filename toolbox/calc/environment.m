function [ vinf, rhoTot, s, Rmean, Talt ] = environment( h, lat, lon, dayOfYear, UTseconds, f107Average, f107Daily, magneticIndex )
%ENVIRONMENT Computes a number of atmospheric parameters using the MATLAB atmosnrlmsise00 function
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
%   magneticIndex   : An 1x7 array of magnetic index (AP) information
%
% Outputs:
%   vinf            : Free flow bulk velocity [m s^-1]
%   rhoTot          : Total atmospheric density in [kg m^-3]
%   s               : Speed ratio [-]
%   Rmean           : Specific gas constant [J kg^-1 K^-1]
%   Talt            : Temperature at altitude [K]
%
% Author: Nicholas H. Crisp
% The University of Manchester
% October 2018
%
%------------- BEGIN CODE -------------

% Constants
[data] = astrophysicalConstants;

% Atmospheric properties
[T, rho] = atmosnrlmsise00(h, lat, lon, 0, dayOfYear, UTseconds, f107Average, f107Daily, magneticIndex);

% Format density and temperature data
Talt = T(2);
rhoTot = rho(6);

% Calculate mean molecular mass [g mol^-1]
mbar = (data.constants.mHe * rho(1)...
    + data.constants.mO * rho(2)... 
    + data.constants.mN2 * rho(3)...
    + data.constants.mO2 * rho(4)...
    + data.constants.mAr * rho(5)...
    + data.constants.mH *  rho(7)...
    + data.constants.mN * rho(8))...
    /(rho(1)+rho(2)+rho(3)+rho(4)+rho(5)+rho(7)+rho(8));

% Calculate specific gas constant [J kg^-1 K^-1]
Rmean = (data.constants.R/mbar)*1000;

% Orbital velocity [m s^-1]
vinf = sqrt(data.constants.mu/(data.constants.R_E+h));

% Thermal velocity
mmean = mbar/data.constants.NA/1000; % Mean mass [kg]
vth = sqrt(2*data.constants.kb*T(2)/mmean);

%Speed ratio [-]
s = vinf/vth;

%------------- END OF CODE --------------