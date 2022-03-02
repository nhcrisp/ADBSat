function [cp, ctau, cd, cl] = coeff_CLL(param_eq, delta)
%COEFF_CLL computes the aerodynamic coefficients using CLL model
%Analytical implementation provided by Walker et al. from the S&C/CLL model
%
% Inputs:
%       param_eq.alphaN : normal thermal energy accommodation coefficient;
%       param_eq.sigmaT : tangential momentum accommodation coefficient;
%       param_eq.Tw     : wall temperature [K];
%       param_eq.Tinf   : free stream temperature [K];
%       param_eq.rho    : structure containing density value [provided by atmosnrlmsise00 model for the termosphere]
%                         param_eq.rho(i) with i = 1, 2, 3, 4, 5, 7, 8, 9 in [1/m3]
%                         param_eq.rho(6) in [kg/m3]
%       delta           : angle between surface normal and the flow
%
% Outputs:
%       cp     : normal pressure coefficient;
%       ctau   : shear stress coefficient;
%       cd     : drag coefficient;
%       cl     : lift coefficient
%
% Author: Sabrina Livadiotti
% The University of Manchester
% April 2020
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

% Constants
[data] = ADBSatConstants;

M_j = [data.constants.mHe, data.constants.mO, data.constants.mN2,...
    data.constants.mO2, data.constants.mH, data.constants.mN]; % Molecular mass of {He, O, N2, O2, H, N} [g/mol]

alphaN = param_eq.alphaN;
sigmaT = param_eq.sigmaT;
%s = param_eq.s;
Vinf = param_eq.vinf;
Tw = param_eq.Tw;
Tinf = param_eq.Tinf;
rho = param_eq.rho;

rhoTot = rho(6);
rhoAr = rho(5);
rhoOa = rho(9);

rho(5) = 0; % density of Argon (not included in the model)
rho(6) = 0; % total density
rho(9) = 0; % density of anomalous Oxygen (not included in the model)
rho(rho == 0) = [];

% Initialise variables:
ctau_j = zeros(length(rho),length(delta)); % shear stress coefficient
cp_j = zeros(length(rho),length(delta));   % normal pressure coefficient
cd_j = zeros(length(rho),length(delta));   % drag coefficient
cl_j = zeros(length(rho),length(delta));   % lift coefficient
m_j = zeros(length(rho),1);
xi_j = zeros(length(rho),1);

for i = 1: length(rho)
    s(i) = Vinf./sqrt(2*(data.constants.kb/(M_j(i)/data.constants.NA/1000)*Tinf));
    % compute fitted parameters for the species considered [ref: Walker et al. (2014)]
    [beta_fp, gamma_fp, delta_fp, zeta_fp] = Fitted_Parameters(i, alphaN);
    % define aux functions
    x_var = s(i).*cos(delta);
    gamma1 = 1/sqrt(pi).*(x_var.*exp(-x_var.^2) + sqrt(pi)/2*(1+2*x_var.^2).*(1+erf(x_var)));
    gamma2 = 1/sqrt(pi).*(exp(-x_var.^2) + sqrt(pi).*x_var.*(1+erf(x_var)));
    % compute normal and shear stress coefficients [ref: Mostaza-Prieto PhD Thesis 2017]
    if alphaN < 1
        cp_j(i,:) = (1/s(i)^2).*((1+sqrt(1-alphaN)).*gamma1 +...
            0.5*(exp(-beta_fp.*(1-alphaN).^gamma_fp).*(Tw/Tinf).^delta_fp.*(zeta_fp./s(i))).*...
            sqrt(Tw/Tinf)*sqrt(pi).*gamma2);
        ctau_j(i,:) = (sigmaT.*sin(delta))./s(i).*gamma2;
    else % sigmaN = alphaN [ref: Walker et al. (2014)]
        cp_j(i,:) = 1/s(i)^2.*(((2-alphaN).*s(i)/sqrt(pi).*cos(delta)+alphaN/2*(Tw/Tinf)^0.5).*exp(-s(i)^2.*(cos(delta)).^2)+...
            ((2-alphaN).*(0.5+s(i)^2.*(cos(delta)).^2)+alphaN/2.*(Tw/Tinf)^0.5.*sqrt(pi).*s(i).*cos(delta)).*(1+erf(s(i).*cos(delta))));
        
        ctau_j(i,:) = sigmaT.*sin(delta)/(s(i)*sqrt(pi)).*(exp(-s(i)^2.*(cos(delta)).^2)+ s(i).*sqrt(pi).*cos(delta).*(1+erf(s(i).*cos(delta))));
    end
    
    % compute drag coefficient:
    cd_j(i,:) = cp_j(i,:).*cos(delta) + ctau_j(i,:).*sin(delta);
    % compute lift coefficient:
    cl_j(i,:) = cp_j(i,:).*sin(delta) - ctau_j(i,:).*cos(delta);
    
end

% Computing total Cd and Cl
% The total aerodynamic coefficients are computed as the weighted sum of the aerodynamic
% coefficients for each molecular species (O2, N2, O, N, He, H) present in the gas
% mixture [ref: Walker et al]

m_avg = 0;
sum_ctau = 0;
sum_cp = 0;
sum_cd = 0;
sum_cl = 0;

for i = 1: length(rho)
    
    m_j(i) = (M_j(i)/data.constants.NA*rho(i)); % [g]
    % Mole/number fraction of each species in the mixture {He, O, N2, O2, H, N}
    % In the computation the presence of molecular species ignored by the
    % analytical implementation is considered
    xi_j(i) = rho(i) / (sum(rho)+rhoAr+rhoOa);
    m_avg = m_avg + xi_j(i) * m_j(i); % Average mass of the mixture
    
    sum_ctau = sum_ctau + xi_j(i) * m_j(i).* ctau_j(i,:);
    sum_cp = sum_cp + xi_j(i) * m_j(i).* cp_j(i,:);
    sum_cd = sum_cd + xi_j(i) * m_j(i).* cd_j(i,:);
    sum_cl = sum_cl + xi_j(i) * m_j(i).*cl_j(i,:);
    
end

% Total coefficients:
ctau = 1/m_avg*sum_ctau;
cp   = 1/m_avg*sum_cp;
cd   = 1/m_avg*sum_cd;
cl   = 1/m_avg*sum_cl;
end

function [beta, gamma, delta, zeta] = Fitted_Parameters(i, alphaN)
% Fitted parameters for CLL model From A.C. Walker, P.M. Mehta, J. Koller,
% Drag Coefficient Model Using the Cercignani–Lampis–Lord Gas–Surface
% Interaction Model, J. Spacecr. Rockets. 51 (2014) 1544–1563.
% doi:10.2514/1.A32677
beta = zeros(1,length(alphaN));
gamma = zeros(1,length(alphaN));
delta = zeros(1,length(alphaN));
zeta = zeros(1,length(alphaN));

if i == 1 % molecular species: He
    beta(alphaN >= 0.95 & alphaN <= 1) = 6.2;
    gamma(alphaN >= 0.95 & alphaN <= 1) = 0.38;
    delta(alphaN >= 0.95 & alphaN <= 1) = 3.3;
    zeta(alphaN >= 0.95 & alphaN <= 1) = 0.74;
    %
    beta(alphaN >= 0.90 & alphaN <= 0.95) = 3.8;
    gamma(alphaN >= 0.90 & alphaN <= 0.95) = 0.52;
    delta(alphaN >= 0.90 & alphaN <= 0.95) = 3.4;
    zeta(alphaN >= 0.90 & alphaN <= 0.95) = 1.12;
    %
    beta(alphaN >= 0.50  & alphaN <= 0.90) = 3.45;
    gamma(alphaN >= 0.50  & alphaN <= 0.90) = 0.52;
    delta(alphaN >= 0.50  & alphaN <= 0.90) = 2.4;
    zeta(alphaN >= 0.50  & alphaN <= 0.90) = 0.93;
    %
    beta(alphaN >= 0  & alphaN <= 0.50) = 0.08;
    gamma(alphaN >= 0  & alphaN <= 0.50) = 0.52;
    delta(alphaN >= 0  & alphaN <= 0.50) = 4.2;
    zeta(alphaN >= 0  & alphaN <= 0.50) = 1.1;
elseif i == 2 % molecular species: O
    beta = 5.85;
    gamma = 0.2;
    delta = 0.48;
    zeta = 31;
elseif i == 3 % molecular species: N2
    beta = 6.6;
    gamma = 0.22;
    delta = 0.48;
    zeta = 35;
elseif i == 4 % molecular species: O2
    beta = 6.3;
    gamma = 0.26;
    delta = 0.42;
    zeta = 20.5;
elseif i == 5 % molecular species: H
    beta(alphaN >= 0.95 & alphaN <= 1) = 3.9;
    gamma(alphaN >= 0.95 & alphaN <= 1) = 0.195;
    delta(alphaN >= 0.95 & alphaN <= 1) = 1.4;
    zeta(alphaN >= 0.95 & alphaN <= 1) = 0.3;
    %
    beta(alphaN >= 0.90 & alphaN <= 0.95) = 3.5;
    gamma(alphaN >= 0.90 & alphaN <= 0.95) = 0.42;
    delta(alphaN >= 0.90 & alphaN <= 0.95) = 2;
    zeta(alphaN >= 0.90 & alphaN <= 0.95) = 0.72;
    %
    beta(alphaN >= 0.50  & alphaN <= 0.90) = 3.45;
    gamma(alphaN >= 0.50  & alphaN <= 0.90) = 0.52;
    delta(alphaN >= 0.50  & alphaN <= 0.90) = 2.4;
    zeta(alphaN >= 0.50  & alphaN <= 0.90) = 0.93;
    %
    beta(alphaN >= 0  & alphaN <= 0.50) = 0.095;
    gamma(alphaN >= 0  & alphaN <= 0.50) = 0.465;
    delta(alphaN >= 0  & alphaN <= 0.50) = 2.9;
    zeta(alphaN >= 0  & alphaN <= 0.50) = 0.92;
elseif i == 6 % molecular species: N
    beta = 4.9;
    gamma = 0.32;
    delta = 0.42;
    zeta = 8;
end

end

%------------- END OF CODE --------------
