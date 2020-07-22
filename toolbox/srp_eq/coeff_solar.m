function [cn, cs] = coeff_solar(delta, param_eq)
%COEFF_SOLAR Calculates solar coefficients for a flat plate using Luthcke et al 1997 formula
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
% Cotober 2018
%
%------------- BEGIN CODE --------------

cn = 2 * ( (param_eq.sol_cD/3)*cos(delta) + param_eq.sol_cR * cos(delta).^2 );
cs = (1 - param_eq.sol_cR)*cos(delta);
cs(cs<0) = 0;

%------------- END OF CODE --------------