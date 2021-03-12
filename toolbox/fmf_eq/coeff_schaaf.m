function [cp, ctau, cd, cl] = coeff_schaaf(param_eq, delta)
%COEFF_SCHAAF calculates aerodynamic coefficients for a flat plate using Schaaf and Chambre formula
%
% Inputs:
%       param_eq.sigmaN : Normal accomodation coeffcient [-]
%       param_eq.sigmaT : Tangential accomodation coeffcient [-]
%       param_eq.s      : Thermal speed ratio [-]
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
% Author: David Mostaza-Prieto
% The University of Manchester
% September 2012
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
%------------- BEGIN CODE --------------

sigmaN = param_eq.sigmaN;
sigmaT = param_eq.sigmaT;
s = param_eq.s;
Tw = param_eq.Tw;
Tinf = param_eq.Tinf;

cp   = 1/s^2.*(((2-sigmaN).*s/sqrt(pi).*cos(delta)+sigmaN/2*(Tw/Tinf)^0.5).*exp(-s^2.*(cos(delta)).^2)+...
((2-sigmaN).*(0.5+s^2.*(cos(delta)).^2)+sigmaN/2.*(Tw/Tinf)^0.5.*sqrt(pi).*s.*cos(delta)).*(1+erf(s.*cos(delta))));

ctau = sigmaT.*sin(delta)/(s*sqrt(pi)).*(exp(-s^2.*(cos(delta)).^2)+ s.*sqrt(pi).*cos(delta).*(1+erf(s.*cos(delta))));

cd   = cp.*cos(delta) + ctau.*sin(delta);
cl   = cp.*sin(delta) - ctau.*cos(delta);

%------------- END OF CODE --------------
