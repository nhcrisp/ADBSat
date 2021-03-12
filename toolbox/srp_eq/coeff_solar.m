function [cn, cs] = coeff_solar(delta, param_eq)
% Calculates solar coefficients for a flat plate using Luthcke et al 1997 formula
%Adsorptivity (alpha) + Specular Reflectivity (rho) + Diffuse Reflectivity (delta) = 1
%Transmissivity = 0;
%
% Inputs:
%    param_eq.cR    : Specular reflectivity component
%    param_eq.cD    : Diffuse refelectivity component
%    delta          : Angle between the flow and the surface normal(external) [rad]
%
% Outputs:
%    cn  : normal coefficient [-]
%    cs  : incident coefficient [-]
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
%------------- BEGIN CODE --------------

cn = 2 * ( (param_eq.sol_cD/3)*cos(delta) + param_eq.sol_cR * cos(delta).^2 );
cs = (1 - param_eq.sol_cR)*cos(delta);
cs(cs<0) = 0;

%------------- END OF CODE --------------
