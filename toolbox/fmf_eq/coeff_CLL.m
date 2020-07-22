function [cp, ctau, cd, cl] = coeff_CLL(param_eq, delta)
%COEFF_CLL computes the aerodynamic coefficients using CLL model
%Analytical implementation provided by Walker et al. from the S&C/CLL model
%
% Inputs:
%       param_eq.alphaN : normal thermal energy accommodation coefficient;
%       param_eq.sigmaT : tangential momentum accommodation coefficient;
%       param_eq.Vw     : velocity of the reflected diffuse molecules [Vw = sqrt(pi.*R.*Tw/2)];
%       param_eq.V      : thermal velocity
%       param_eq.rho    : structure containing density value [provided by atmosnrlmsise00 model for the termosphere]
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
%------------- BEGIN CODE --------------

alphaN = param_eq.alphaN;
sigmaT = param_eq.sigmaT;
s = param_eq.s;
Vw = param_eq.Vw;
V = param_eq.V;
rho = param_eq.rho;

rhoTot = rho(6);

rho(5) = 0; % density of Argon (not included in the model)
rho(6) = 0; % total density (saved into a separate variable [rhoTot])
rho(9) = 0; % density of anomalous Oxygen (not included in the model)
rho(rho == 0) = [];

% Initialise variables:
ctau_j = zeros(length(rho),length(delta)); % shear stress coefficient
cp_j = zeros(length(rho),length(delta));   % normal pressure coefficient
cd_j = zeros(length(rho),length(delta));   % drag coefficient
cl_j = zeros(length(rho),length(delta));   % lift coefficient
m_j = zeros(length(rho),1);
n_j = zeros(length(rho),1);
xi_j = zeros(length(rho),1);

for i = 1: length(rho) % compute fitted parameters for the species considered [ref: Walker et al.]
    [beta_fp, gamma_fp, delta_fp, zeta_fp] = Fitted_Parameters(i, alphaN);
    % compute shear stress coefficient
    ctau_j(i,:) = sigmaT./(sqrt(pi)*s).*sin(2*delta).*cos(delta).*exp(-s^2*cos(delta).^2)+ ...
        2*sigmaT./(sqrt(pi)*s).*sin(delta).^3.*exp(-s^2*cos(delta).^2) + ...
        sigmaT.*sin(2*delta).*erf(s*cos(delta));
    % compute normal stress coefficient
    cp_j(i,:) = 4*cos(delta)./(sqrt(pi)*s).*exp(-s^2.*cos(delta).^2) - ...
        2*cos(delta)./(sqrt(pi)*s)*(1-sqrt(1-alphaN)).*exp(-s^2*cos(delta).^2) + ...
        (((1+sqrt(1-alphaN))/s^2) + 4.*cos(delta).^2 - 2*(1-sqrt(1-alphaN)).*cos(delta).^2).*erf(s*cos(delta)) + ...
        exp(-beta_fp*(1-alphaN)^gamma_fp)*2.*cos(delta).*zeta_fp/s*(2*Vw*s/(sqrt(pi)*V))^(2*delta_fp)*Vw/V;
    
    %     ind = find(delta>pi/2);
    %     ctau(ind) = 0;
    %     cp(ind) = 0;
    
    % compute drag coefficient:
    cd_j(i,:) = cp_j(i,:).*cos(delta) + ctau_j(i,:).*sin(delta);
    % compute lift coefficient:
    cl_j(i,:) = cp_j(i,:).*sin(delta) - ctau_j(i,:).*cos(delta);
end

% Computing total Cd and Cl
% The total aerodynamic coefficients are computed as the weighted sum of the aerodynamic
% coefficients for each molecular species (O2, N2, O, N, He, H) present in the gas
% mixture [ref: Walker et al]

M_j = [4.002602, 15.999, 14.0067*2, 15.999*2, 1.00784, 14.0067]; % Molecular mass of {He, O, N2, O2, H, N} [g/mol]
mw_mix = sum(M_j); % Molar weigth of the mixture [g/mol]

% Not precise: rhoTot comprises Argon and anomalous Oxygen at the moment!
mass_mix = 1000*rhoTot*1; % Mass of the gas mixture [g]
n_mix = mass_mix/mw_mix; % Number of moles for the gas mixture [mol]

m_avg = 0;
sum_ctau = 0;
sum_cp = 0;
sum_cd = 0;
sum_cl = 0;

for i = 1: length(rho)
    
    m_j(i) = 1000*rho(i)*1;
    n_j(i) = m_j(i) / M_j(i);  % Number of moles for each species in the mixture {He, O, N2, O2, H, N} [mol]
    xi_j(i) = n_j(i) / n_mix; % Mole Fraction of each species in the mixture {He, O, N2, O2, H, N}
    m_avg = m_avg + xi_j(i) * M_j(i); % Average mass of the mixture
    
    sum_ctau = sum_ctau + xi_j(i) * M_j(i).* ctau_j(i,:);
    sum_cp = sum_cp + xi_j(i) * M_j(i).* cp_j(i,:);
    sum_cd = sum_cd + xi_j(i) * M_j(i).* cd_j(i,:);
    sum_cl = sum_cl + xi_j(i) * M_j(i).*cl_j(i,:);
    
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

if i == 1
    % molecular species: He
    if alphaN >= 0.95 && alphaN <= 1
        beta = 6.2; gamma = 0.38; delta = 3.3; zeta = 0.74;
    elseif alphaN >= 0.90 && alphaN <= 0.95
        beta = 3.8; gamma = 0.52; delta = 3.4; zeta = 1.12;
    elseif  alphaN >= 0.50  && alphaN <= 0.90
        beta = 3.45; gamma = 0.52; delta = 2.4; zeta = 0.93;
    elseif  alphaN >= 0  && alphaN <= 0.50
        beta = 0.08; gamma = 0.52; delta = 4.2; zeta = 1.1;
    end
elseif i == 2
    % molecular species: O
    beta = 4.4; gamma = 0.32; delta = 0.48; zeta = 11;
elseif i == 3
    % molecular species: N2
    beta = 5.5; gamma = 0.18; delta = 0.5; zeta = 51;
elseif i == 4
    % molecular species: O2
    beta = 5.45; gamma = 0.18; delta = 0.5; zeta = 49;
elseif i == 5
    % molecular species: H
    if alphaN >= 0.95 && alphaN <= 1
        beta = 3.9; gamma = 0.195; delta = 1.4; zeta = 0.3;
    elseif alphaN >= 0.90 && alphaN <= 0.95
        beta = 3.5; gamma = 0.42; delta = 2; zeta = 0.72;
    elseif  alphaN >= 0.50  && alphaN <= 0.90
        beta = 3.45; gamma = 0.52; delta = 2.4; zeta = 0.93;
    elseif  alphaN >= 0  && alphaN <= 0.50
        beta = 0.095; gamma = 0.465; delta = 2.9; zeta = 0.92;
    end
elseif i == 6
    % molecular species: N
    beta = 4.75; gamma = 0.24; delta = 0.5; zeta = 20;
end

end

%------------- END OF CODE --------------