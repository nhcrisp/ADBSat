function [cp, ctau, cd, cl] = coeff_storchHyp(param_eq, delta)
% Calculates hyperthermal aerodynamic coefficients for a flat plate using Storch's formula
% J. A. Storch. "Aerodynamic disturbances on spacecraftin free-molecular flow (monograph)". Technical report, The Aerospace Corporation, 2002.
%
% Inputs:
%       param_eq.sigmaN : Normal accomodation coeffcient [-]
%       param_eq.sigmaT : Tangential accomodation coeffcient [-]
%       param_eq.Vw     : Average normal velocity of diffusely reflected molecules [m/s]
%       param_eq.V      : Incident velocity [m/s]
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
%------------- BEGIN CODE -------------

sigmaN = param_eq.sigmaN;
sigmaT = param_eq.sigmaT;
Tw = param_eq.Tw;
Rmean = param_eq.Rmean;
Vw = sqrt((pi*Rmean*Tw)/2);
V = param_eq.vinf;

cp   = 2.*cos(delta).*( sigmaN.*(Vw/V) + (2 - sigmaN).*cos(delta));
ctau = 2.*cos(delta).*sin(delta).*sigmaT;

ind = find(delta>pi/2);
cp(ind) = 0;
ctau(ind) = 0;

cd = cp.*cos(delta) + ctau.*sin(delta);
cl = cp.*sin(delta) - ctau.*cos(delta);

%------------- END OF CODE --------------

% cd = 2.*cos(delta).*(sigmaT + sigmaN.*(Vw/V).*cos(delta) + (2 - sigmaN - sigmaT).*cos(delta).^2);
% cl = 2.*cos(2.*delta).*(sigmaN.*(Vw/V) + (2 - sigmaN - sigmaT).*cos(delta));
% 
% cp   = cd.*cos(delta) - cl.*sin(delta);
% ctau = cd.*sin(delta) + cl.*cos(delta);
