function [beta, gamma, delta, zeta] = Fitted_Parameters(i, alphaN)
% Fitted parameters for CLL model From A.C. Walker, P.M. Mehta, J. Koller,
% Drag Coefficient Model Using the CercignaniLampisLord GasSurface
% Interaction Model, J. Spacecr. Rockets. 51 (2014) 15441563.
% doi:10.2514/1.A32677
%
% Original: 
% Modif: 7/09/2022 A.Macario-Rojas, Multimaterial handling
% beta, gamma, delta zeta, from cosntants to vectors to handle
% multimaterial alphaN

% Initalize arrays
temp = ones(1,length(alphaN));
beta = temp; 
gamma = temp;
delta = temp;
zeta = temp;

if i == 1
    % molecular species: He, FLAT PLATE VERIFIED
    iacc = and(alphaN > 0.95,alphaN <= 1);
    beta(iacc) = 6.2; gamma(iacc) = 0.38; delta(iacc) = 3.3; zeta(iacc) = 0.74;

    iacc = and((alphaN > 0.90),(alphaN <= 0.95));
    beta(iacc) = 3.8; gamma(iacc) = 0.52; delta(iacc) = 3.4; zeta(iacc) = 1.12;

    iacc = and((alphaN > 0.50),(alphaN <= 0.90));
    beta(iacc) = 3.45; gamma(iacc) = 0.52; delta(iacc) = 2.4; zeta(iacc) = 0.93;

    iacc = and((alphaN > 0),(alphaN <= 0.50));
    beta(iacc) = 0.08; gamma(iacc) = 0.52; delta(iacc) = 4.2; zeta(iacc) = 1.1;

    % Original code        
    % if and(alphaN >= 0.95,alphaN <= 1)
    %    beta = 6.2; gamma = 0.38; delta = 3.3; zeta = 0.74;
    %elseif and((alphaN >= 0.90),(alphaN <= 0.95))
    %    beta = 3.8; gamma = 0.52; delta = 3.4; zeta = 1.12;
    %elseif  and((alphaN >= 0.50),(alphaN <= 0.90))
    %    beta = 3.45; gamma = 0.52; delta = 2.4; zeta = 0.93;
    %elseif  and((alphaN >= 0),(alphaN <= 0.50))
    %    beta = 0.08; gamma = 0.52; delta = 4.2; zeta = 1.1;
    %end
elseif i == 2
    % molecular species: O, FLAT PLATE VRIFIED 
    beta = 5.85.*beta; gamma = 0.2.*gamma; delta = 0.48.*delta; zeta = 31.*zeta;

    % Original code 
    % beta = 4.4; gamma = 0.32; delta = 0.48; zeta = 11; % SPHERE!!!
elseif i == 3
    % molecular species: N2, FLAT PLATE VERIFIED
    beta = 6.6.*beta; gamma = 0.22.*gamma; delta = 0.48.*delta; zeta = 35.*zeta;

    % Original code 
    %beta = 5.5; gamma = 0.18; delta = 0.5; zeta = 51; % SPHERE!!!
elseif i == 4
    % molecular species: O2, FLAT PLATE VERIFIED
    beta = 6.3.*beta; gamma = 0.26.*gamma; delta = 0.42.*delta; zeta = 20.5.*zeta;

    % Original code 
    %beta = 5.45; gamma = 0.18; delta = 0.5; zeta = 49;
elseif i == 5
    % molecular species: H, FLAT PLATE VERIFIED
    iacc = and(alphaN > 0.95,alphaN <= 1);
    beta(iacc) = 3.9; gamma(iacc) = 0.195; delta(iacc) = 1.4; zeta(iacc) = 0.3;

    iacc = and((alphaN > 0.90),(alphaN <= 0.95));
    beta(iacc) = 3.5; gamma(iacc) = 0.42; delta(iacc) = 2.0; zeta(iacc) = 0.72;

    iacc = and((alphaN > 0.50),(alphaN <= 0.90));
    beta(iacc) = 3.45; gamma(iacc) = 0.52; delta(iacc) = 2.4; zeta(iacc) = 0.93;

    iacc = and((alphaN > 0),(alphaN <= 0.50));
    beta(iacc) = 0.095; gamma(iacc) = 0.465; delta(iacc) = 2.9; zeta(iacc) = 0.92;

    % Original code 
    %if and(alphaN >= 0.95,alphaN <= 1)
    %    beta = 3.9; gamma = 0.195; delta = 1.4; zeta = 0.3;
    %elseif and(alphaN >= 0.90,alphaN <= 0.95)
    %    beta = 3.5; gamma = 0.42; delta = 2; zeta = 0.72;
    %elseif and(alphaN >= 0.50,alphaN <= 0.90)
    %    beta = 3.45; gamma = 0.52; delta = 2.4; zeta = 0.93;
    %elseif and(alphaN >= 0,alphaN <= 0.50)
    %    beta = 0.095; gamma = 0.465; delta = 2.9; zeta = 0.92;
    %end
elseif i == 6
    % molecular species: N, FLAT PLATE VERIFIED
    beta = 4.9.*beta; gamma = 0.32.*gamma; delta = 0.42.*delta; zeta = 8.*zeta;

    % Original code 
    % beta = 4.75; gamma = 0.24; delta = 0.5; zeta = 20; % SPHERE!!!
end
