function [cp, ctau, cd, cl] = coeff_DRIA(param_eq, delta)
% Calculates aerodynamic coefficients for a flat plate using DRIA (Diffuse
% Reflection with Incomplete Accommodation) 
%
% Inputs:
%       param_eq.alpha  : Accomodation coefficient
%       param_eq.s      : Thermal speed ratio [-]
%       param_eq.Tw     : Wall temperature [K]
%       delta           : Angle between the flow and the surface normal(external) [rad]
%
% Outputs:
%       cp      : Pressure coefficient [-]
%       ctau    : Shear stress coefficient [-]
%       cd      : Drag coefficent [-]
%       cl      : Local Lift coefficient [-]
%
% Author: Nicholas H. Crisp
% The University of Manchester
% July 2021
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

[data] = astrophysicalConstants;

%Tinf = param_eq.Tinf;
Vinf = param_eq.V;
alpha = param_eq.alpha;
s = param_eq.s;
Tw = param_eq.Tw;
rho = param_eq.rho;

rhoHe = rho(1);
rhoO = rho(2);
rhoN2 = rho(3);
rhoO2 = rho(4);
rhoAr = rho(5);
%rhoTot = rho(6);
rhoH = rho(7);
rhoN = rho(8);
%rhoAnO = rho(9);

% Calculate mean molecular mass
Mmean = (data.constants.mHe*rhoHe + data.constants.mO*rhoO +...
    data.constants.mN2*rhoN2 + data.constants.mO2*rhoO2 +...
    data.constants.mAr*rhoAr + data.constants.mH*rhoH +...
    data.constants.mN*rhoN)./...
    (rhoHe + rhoO + rhoN2 + rhoO2 + rhoAr + rhoH + rhoN); % [g mol^-1]

%T_ki = ((Mmean*1e-3)/data.constants.NA * Vinf^2)/(3*data.constants.kb); % Incoming kinetic temperature [K]
%T_kr  = T_ki*(1-alpha) + alpha*Tw; % Reflected kinetic temperature [K]
%V_r = Vinf * sqrt(2/3) * sqrt(1 + alpha*(Tw/T_ki - 1)); % Reflected particle velocity

% [Doornbos E (2012) Thermospheric Density and Wind Determination from Satellite Dynamics]
P = exp(-cos(delta).^2 * s.^2)./s;
G = 1/(2*s.^2);
Q = 1+G;
Z = 1+erf(cos(delta).*s);
R = data.constants.R/Mmean * 1e3;
Vratio = sqrt((1/2)*(1+alpha.*((4*R*Tw)./Vinf.^2 - 1))); % [Doornbos 2012]

cd = P./sqrt(pi) + cos(delta).*Q.*Z + 0.5*cos(delta).*Vratio.*(cos(delta).*sqrt(pi).*Z+P);
cl = sin(delta).*G.*Z + 0.5*sin(delta).*Vratio.*(cos(delta).*sqrt(pi).*Z+P);

cp   = cd.*cos(delta) + cl.*sin(delta);
ctau = cd.*sin(delta) - cl.*cos(delta);

%------------- END OF CODE --------------
