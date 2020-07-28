function [data] = astrophysicalConstants(data)
%ASTROPHYSICALCONSTANTS creates a structure of useful constants
%
% Author: Nicholas H. Crisp
% The University of Manchester
% October 2019
%
%------------- BEGIN CODE -------------

data.constants.G = 6.67430e-11; % Gravitational constant [m^3/kg s^2]
data.constants.g0 = 9.80665; %Earth Standard Acceleration due to Gravity [m/s^2]
data.constants.mu = 3.986004418e14;  % GM Earth [m^3/s^2] (source: SMAD 3rd edition)
data.constants.R_E = 6.37813649e6; % Equatorial Earth radius [m] (source: SMAD 3rd edition)
data.constants.m_E = 5.9737e24; % Earth Mass [kg] (source: SMAD 3rd edition)
data.constants.w_E = 7.2921158553e-5; % Earth Angular velocity in [rad/s] (source: SMAD 3rd edition)
data.constants.J2 = 0.00108263; % Earth Geopotential Second Zonal Coefficient (source: SMAD 3rd edition)

data.constants.dRAAN_SSO = 1.991063853e-7; % Nodal Precession Rate for SSO [rad/s] (source: Vallado 2013)

data.constants.Q_S = 1367; % Solar irradiance [W/m2]
data.constants.c = 3e8; % Speed of light [m/s]
data.constants.P_SRP = 4.575e-6; % Nominal solar pressure [N/m^2]
data.constants.AU = 1.49597870700e11; % Astronomical Unit [m]

data.constants.R = 8.3144621; % Molar (Universal/Ideal) Gas constant [J K^-1 mol^-1]
data.constants.kb = 1.3806503e-23; % Boltzmann constant [m^2 kg^-2 K^-1]
data.constants.NA = 6.02214076e23; % Avogadro constant [n mol^-1]

data.constants.mHe = 2; % Molecular mass of He [g mol^-1]
data.constants.mO = 15.99; % Molecular mass of O [g mol^-1]
data.constants.mN2 = 28; % Molecular mass of N2 [g mol^-1]
data.constants.mO2 = 32; % Molecular mass of O2 [g mol^-1]
data.constants.mAr = 40; % Molecular mass of Ar [g mol^-1]
data.constants.mH = 1; % Molecular mass of H [g mol^-1]
data.constants.mN = 14; % Molecular mass of N [g mol^-1]

%------------- END OF CODE --------------