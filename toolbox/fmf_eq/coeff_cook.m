function [cp, ctau, cd, cl] = coeff_cook(param_eq, delta)
% Calculates "Hyperthermal" FMF coefficients for a flat plate using Cook's formula
% Koppenwallner. "Comment on special section: New perspectives on the satellite 
% drag environments of Earth, Mars, and Venus". Journal of Spacecraft and Rockets, 2008.
% - Hyperthermal
% - Diffuse
%
% Inputs:
%       param_eq.alpha  : Energy accomodation coeffcient [-]
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

alpha = param_eq.alpha;
Tw = param_eq.Tw;
Tinf = param_eq.Tinf;

cd = 2.*(1 + 2/3.*sqrt(1 + alpha .*(Tw/Tinf-1)).*cos(delta)).*cos(delta);
cl = 4/3.*sqrt(1 + alpha .*(Tw/Tinf-1)).*sin(delta).*cos(delta);

cp   = cd.*cos(delta) + cl.*sin(delta);
ctau = cd.*sin(delta) - cl.*cos(delta);

%------------- END OF CODE --------------
