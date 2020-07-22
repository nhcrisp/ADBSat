function [cp, ctau, cd, cl] = coeff_sentman(param_eq, delta)
%COEFF_SENTMAN calculates aerodynamic coefficients for a flat plate using Sentman's formula
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
%------------- BEGIN CODE -------------

Tinf = param_eq.Tinf;
alpha = param_eq.alpha;
s = param_eq.s;
Tw = param_eq.Tw;

Ti = 2/3*s.^2*Tinf;

cp = ((cos(delta)).^2 + 1./(2*s.^2)).*(1+erf(s.*cos(delta)))+ ...
    cos(delta)./(sqrt(pi)*s).*exp(-s.^2.*(cos(delta)).^2) + ...
    0.5*sqrt(2/3*(1+alpha*(Tw./Ti-1))).*[sqrt(pi)*cos(delta).*(1+erf(s.*cos(delta)))+1./s.*exp(-s.^2.*(cos(delta)).^2)];
    
ctau = sin(delta).*cos(delta).*(1+erf(s.*cos(delta))) + sin(delta)./(s*sqrt(pi)).*(exp(-s.^2.*(cos(delta)).^2));


cd   = cp.*cos(delta) + ctau.*sin(delta);
cl   = cp.*sin(delta) - ctau.*cos(delta);

%------------- END OF CODE --------------
