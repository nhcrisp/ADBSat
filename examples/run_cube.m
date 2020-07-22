% Cube Example
clear

modName = 'cube';
% Path to model file
ADBSat_path = ADBSat_dynpath;
modIn = fullfile(ADBSat_path,'inou','obj_files',[modName,'.obj']);

%Input conditions
alt = 200; %km
inc = 51.6; %deg
env = [alt*1e3, inc/2, 0, 106, 0, 65, 65, ones(1,7)*3]; % Environment variables

aoa = [0 45]; % Angle of attack
aos = [0 45]; % Angle of sideslip

% Model parameters
shadow = 1;
inparam.gsi_model = 'sentman';
inparam.alpha = 1; % Accommodation (altitude dependent)
inparam.Tw = 300; % Wall Temperature [K]

solar = 1;
inparam.sol_cR = 0.15; % Specular Reflectivity
inparam.sol_cD = 0.25; % Diffuse Reflectivity

verb = 1;
del = 0;

% Import model
[modName, modOut] = ADBSatImport(modIn, verb);

% Calculate
[ADBout] = ADBSatFcn(modName, inparam, aoa, aos, shadow, solar, env, del, verb);

% Plot
if verb && ~del
    plot_surfq(ADBout, modOut, aoa(1), aos(1), 'cd');
end