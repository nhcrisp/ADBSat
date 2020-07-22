function [cp, ctau, cd, cl] = coeff_storchHyp(param_eq, delta)
%COEFF_STORCHHYP calculates hyperthermal aerodynamic coefficients for a flat plate using Storch's formula
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
%------------- BEGIN CODE -------------

sigmaN = param_eq.sigmaN;
sigmaT = param_eq.sigmaT;
Vw = param_eq.Vw;
V = param_eq.V;

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