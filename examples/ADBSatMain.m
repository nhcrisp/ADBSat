%ADBSATMAIN Creates a .mat file(s) in "/inou/results" with the following fields:
%
% Inputs:
%       filename : name of the input model file (.mat) in /inou/models
%       model    : GSI model
%       shadow   : flag for shadown analysis
%       solar    : flag for solar coefficient analysis
%       env      : vector of input environmental parameters 
%       verb     : flag for visual output
%       varargin : GSI parameter input (dependent on GSI model)
%
% Outputs:
%       fileout  : .mat file(s) with the following items in a structure variable
%           aoa      : Angle of attack
%           aos      : Angle of sideslip
%           Aref     : Reference area = Area total/2
%           Lref     : Reference longitude = Xmax-Xmin
%           AreaProj : Projected area
%           aero
%               Cf_w     : Total aerodynamic force coeffients in wind axes [x,y,z] (3x1)
%               Cf_f     : Total aerodynamic force coeffients in flight axes [x,y,z] (3x1)
%               Cm_B     : Total aerodynamic moment coeffients in body axes [x,y,z](3x1)
%           solar
%               C_s      : Total solar force coefficients in 
%               C_sB     : Total solar moment coefficients in 
%           param_eq : Structure containing input GSI parameters
%
% Author: David Mostaza-Prieto
% The University of Manchester
% November 2012
%
%------------- BEGIN CODE --------------

clear

fiName = 'cube'; %name of the .mat file in /inou/models

flag_shadow = 1; % Shadow analysis 1:perform, 0:no
flag_solar = 1; % Include solar coefficient calculation

% Atmopsheric environment
[inparam.V, ~, inparam.s, Rmean, inparam.Tinf] = environment(200e3, 45, 0, 1, 12*3600, 250, 250, ones(1,7)*7);
inparam.Tw = 300; % Surface temperature     
inparam.Vw = sqrt(pi.*Rmean.*inparam.Tw/2); % Velocity of the reflected diffuse molecules

% GSI Model Inputs
inparam.gsi_model = 'sentman'; % ('schaaf', 'storchHyp', 'cook', 'newton', 'CLL')
% Vectors for multiple materials
%sigmaT = 1; % Tangent accomodation coefficient
%sigmaN = 1;   % Normal accomodation coeffcient
inparam.alpha = 1;   % Thermal accomodation coefficient

% Intrinsic rotation angles (angle of attack then sideslip)
aoa = [0]*pi/180; % Angle(s) of attack
aos = [10]*pi/180; % Angle(s) of sideslip

if flag_solar == 1
    inparam.sol_cD = 0.25;
    inparam.sol_cR = 0.15;
end

% Calculate
fileOut = calc_coeff(fiName, aoa, aos, inparam, flag_shadow, flag_solar, 1, 0); 

%------------- END OF CODE --------------