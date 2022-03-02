function [cp, ctau, cd, cl] = coeff_sentman(param_eq, delta)
% Calculates aerodynamic coefficients for a flat plate using Sentman's [1]
% formula. Assumes the Schamberg's [2] modification to account for incomplete
% accommodation and the Moe et al. [3] correction factor of sqrt(2/3).
%
% [1] Sentman LH (1961) Free molecule flow theory and its application to
% the determination of aerodynamic forces. Sunnyvale, CA
%
% [2] Schamberg R (1959) A New Analytic Representation of Surface
% Interaction for Hyperthermal Free Molecule Flow with Applications to
% Neutral-particle Drag Estimates of Satellites. Rand Corporation
%
% [3] 1. Moe K, Moe MM, Rice CJ (2004) Simultaneous Analysis of
% Multi-Instrument Satellite Measurements of Atmospheric Density. J Spacecr
% Rockets 41:849–853.https://doi.org/10.2514/1.2090
%
% Inputs:
%       param_eq.alpha  : Accomodation coefficient
%       param_eq.s      : Thermal speed ratio [-]
%       param_eq.Tw     : Wall temperature [K]
%       param_eq.Tatm   : Flow temperature [K]
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
% February 2013
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

Tinf = param_eq.Tinf;
alpha = param_eq.alpha;
s = param_eq.s;
Tw = param_eq.Tw;

Ti = 2/3*s.^2*Tinf;

cp = ((cos(delta)).^2 + 1./(2*s.^2)).*(1+erf(s.*cos(delta)))+ ...
    cos(delta)./(sqrt(pi)*s).*exp(-s.^2.*(cos(delta)).^2) + ...
    0.5*sqrt(2/3*(1+alpha*(Tw./Ti-1))).*(sqrt(pi)*cos(delta).*(1+erf(s.*cos(delta)))+1./s.*exp(-s.^2.*(cos(delta)).^2));
    
ctau = sin(delta).*cos(delta).*(1+erf(s.*cos(delta))) + sin(delta)./(s*sqrt(pi)).*(exp(-s.^2.*(cos(delta)).^2));

cd   = cp.*cos(delta) + ctau.*sin(delta);
cl   = cp.*sin(delta) - ctau.*cos(delta);

%------------- END OF CODE --------------
