function [alpha] = accom_SESAM(inparam)
% Calculates the accommodation coefficient using the Semi-Empirical Satellite Accommodation Model (SESAM)
% [M. D. Pilinski, B. M. Argrow, S. E. Palo, and B. R. Bowman, “Semi-empirical ...
% satellite ac-commodation model for spherical and randomly tumbling objects,”...
% Journal of Spacecraft and Rockets, vol. 50, no. 3, pp. 556–571, 2013].
% 
% Takes atmospheric parameters from 'environment' function in ADBSat as an
% input and reports accommodation coefficient using the SESAM approach
% Inputs:
%       
%       inparam.s       : Thermal speed ratio [-]
%       inparam.vinf    : Free flow bulk velocity [m s^-1]
%       inparam.m_s     : Surface molecule mass [amu]
%       inparam.K_s     : Substrate coefficient, takes a value between 2.4 and 3.6. 
%       inparam.rho     : an array of M-by-9 values of densities.  These values are
%                         HE number density in meters^-3, O number density in meters^-3,
%                         N2 number density in meters^-3, O2 number density in meters^-3,
%                         AR number density in meters^-3, total mass density in kilogram 
%                         per meters cubed, H number density in meters^-3, N number
%                         density in meters^-3, and Anomalous oxygen number density in
%                         meters^-3.
%
% Outputs:
%       param_eq.alpha  : Accommodation coefficient using SESAM
%
% Author: Mshary R. Aldakheel
% The University of Manchester
% July 2024
%
%--- Copyright notice ---%
% Copyright (C) 2021 The University of Manchester
% Written by David Mostaza Prieto,  Nicholas H. Crisp, Luciana Sinpetru,
% Sabrina Livadiotti and Mshary Aldakheel
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

% Assign value for K_s if not provided by the user
if ~isfield(inparam,'K_s')
    inparam.K_s = 2.4;
    warning('Substrate coefficient (K_s) not provided, default value of 2.4 is used')
end

% Check if K_s is within limits
if inparam.K_s<2.4 || inparam.K_s>3.6
    inparam.K_s = 2.4;
    warning('Substrate coefficient (K_s) value outside the limits [2.4,3.6]. Default value of 2.4 is used')
end

% Assign value for m_s if not provided by the user
if ~isfield(inparam,'m_s')
    inparam.m_s = 65;
    warning('Surface atomic mass (m_s) not provided, default value of 65 amu is used')
end

%-------------------    Constants   --------------
E_b = 5.7/(6.24150974e18);  % Adsorption energy (J)
T_ab = 93.31; % Transition temperature (K)  
K_Lo = 5e06; % Initial Langmuir parameter (Torr^-1)
K_Lf = 3e04; % Final Langmuir parameter (Torr^-1)
data=ADBSatConstants; % Calls ADBSat constants

%-------------------       --------------
rho_O = inparam.rho(2)*(data.constants.mO*1e-03/data.constants.NA); % Computes the density of oxygen 

P_O= 0.5*rho_O*inparam.vinf^2*(((2*inparam.s^2+1)/(sqrt(pi)*inparam.s^3))*exp(-inparam.s^2) ...
    + ((4*inparam.s^4 + 4*inparam.s^2-1) / (2*inparam.s^4))*erf(inparam.s))/(101325/760); % Partial pressure of oxygen in torr

E_r = 0.5*(data.constants.mO*1e-03/data.constants.NA)*inparam.vinf^2; % Computes the incident kinetic energy

zeta = exp(2*sqrt(E_b*E_r)/(data.constants.kb*T_ab)); % Check for
if zeta==inf
    zeta = realmax;
end

s_o = (sqrt(pi*data.constants.kb*T_ab*E_r)*(erf((sqrt(E_b)-sqrt(E_r))/sqrt(data.constants.kb*T_ab))+erf(sqrt(E_r/(data.constants.kb*T_ab))))...
         +data.constants.kb*T_ab*exp(-(E_b+E_r)/(data.constants.kb*T_ab))*(exp(E_b/(data.constants.kb*T_ab))-zeta))...
         /(sqrt(pi*data.constants.kb*T_ab.*E_r)*(erf(sqrt(E_r/(data.constants.kb*T_ab)))+1) + data.constants.kb*T_ab*exp(-E_r/(data.constants.kb*T_ab)));

K_L = s_o*K_Lo+K_Lf;


m_bar = ((inparam.rho(1)*data.constants.mHe+inparam.rho(2)*data.constants.mO+inparam.rho(3)*data.constants.mN2+inparam.rho(4)*data.constants.mO2+ ...
        inparam.rho(5)*data.constants.mAr+inparam.rho(7)*data.constants.mH+inparam.rho(8)*data.constants.mN+inparam.rho(9)*data.constants.mAnO)*1e-03/(data.constants.NA))/ ...
        (inparam.rho(1)+inparam.rho(2)+inparam.rho(3)+inparam.rho(4)+inparam.rho(5)+inparam.rho(7)+inparam.rho(8)+inparam.rho(9)); 

mu = m_bar/(inparam.m_s*1.6605e-027);
alpha_s = inparam.K_s*mu/(1+mu)^2;
theta_ = K_L*P_O/(1+K_L*P_O);

alpha = (1-theta_)*alpha_s+theta_;

end