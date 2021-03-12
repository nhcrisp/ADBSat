function [cp, ctau, cd, cl] = coeff_newton(param_eq, delta)
%COEFF_NEWTON calculates aerodynamic coefficients for a flat plate using Netwon's formula
%
% Inputs:
%       param_eq    : N/A (blank)
%       delta       : Angle between the flow and the surface normal(external) [rad]
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

cp   = 2.*(cos(delta).^2);
ctau = zeros(1,length(delta));

ind = find(delta>pi/2);
cp(ind) = 0;
ctau(ind) = 0;

cd = cp.*cos(delta) + ctau.*sin(delta);
cl = cp.*sin(delta) - ctau.*cos(delta);

%------------- END OF CODE --------------

% cd = 2.*cos(delta).*(sigmaT + sigmaN.*(Vw/V).*cos(delta) + (2 - sigmaN - sigmaT).*cos(delta).^2);
% cl = 2.*cos(2.*delta).*(sigmaN.*(Vw/V) + (2 - sigmaN - sigmaT).*cos(delta));
% cp   = cd.*cos(delta) - cl.*sin(delta);
% ctau = cd.*sin(delta) + cl.*cos(delta);
