% Calculates FMF aerodynamic coefficients according to Maxwell's model.
% Assumes that the fraction of molecules which undergo diffuse reflection is param_eq.alpha,
% and the fraction which undergo specular reflection is (1-param_eq.alpha).
%
% Inputs:
%       param_eq.alpha  : Fraction of molecules which undergo diffuse reflection [-]
%       param_eq.Tw     : Wall temperature [K]
%       param_eq.Tinf   : Kinetic temperature Ti = V^2/3R 
%       delta           : Angle between the flow and the surface normal(external) [rad]
%
% Outputs:
%       cp      : Pressure coefficient [-]
%       ctau    : Shear stress coefficient [-]
%       cd      : Drag coefficent [-]
%       cl      : Local Lift coefficient [-]
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
%----------- BEGIN CODE ---------

function [cp, ctau, cd, cl] = coeff_maxwell(param_eq,  delta)

alpha = param_eq.alpha;
Tw = param_eq.Tw;
Tinf = param_eq.Tinf;
s = param_eq.s;

f = 1-alpha;
theta = pi/2-delta;

cd = 2*((1-f*cos(2*theta))./(sqrt(pi)*s)).*exp(-s^2*sin(theta).^2)...
    + (sin(theta)/s^2).*(1+2*s^2+f*(1-2*s^2*cos(2*theta))).*erf(s*sin(theta))...
    + (1-f)/s * sqrt(pi)*sin(theta).^2*sqrt(Tw/Tinf);

cl = ((4*f)/(sqrt(pi)*s))*sin(theta).*cos(theta).*exp(-s^2*sin(theta).^2)...
    + (cos(theta)/s^2).*(1+f*(1+4*s^2*sin(theta).^2)).*erf(s*sin(theta))...
    + ((1-f)/s)*sqrt(pi)*sin(theta).*cos(theta)*sqrt(Tw/Tinf);

cp   = cd.*cos(delta) + cl.*sin(delta);
ctau = cd.*sin(delta) - cl.*cos(delta);
%----------- END CODE ----------
