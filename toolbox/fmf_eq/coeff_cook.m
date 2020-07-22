function [cp, ctau, cd, cl] = coeff_cook(param_eq, delta)
%COEFF_COOK calculates "Hyperthermal" FMF coefficients for a flat plate using Cook's formula
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
%------------- BEGIN CODE --------------

alpha = param_eq.alpha;
Tw = param_eq.Tw;
Tinf = param_eq.Tinf;

cd = 2.*(1 + 2/3.*sqrt(1 + alpha .*(Tw/Tinf-1)).*cos(delta)).*cos(delta);
cl = 4/3.*sqrt(1 + alpha .*(Tw/Tinf-1)).*sin(delta).*cos(delta);

cp   = cd.*cos(delta) + cl.*sin(delta);
ctau = cd.*sin(delta) - cl.*cos(delta);

%------------- END OF CODE --------------