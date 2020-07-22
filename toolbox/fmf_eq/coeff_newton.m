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