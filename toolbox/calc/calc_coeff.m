function [fileOut] = calc_coeff(fiName, aoaS, aosS, param_eq, flag_shad, flag_sol, del, verb)
%CALC_COEFF Calcultes local and global coefficients for the triangular mesh geometry
% Produces one or several output .m files with the following variables
%
%   aoa      : Angle of attack [rad]
%   aos      : Angle of sideslip [rad]
%   tauDir   : Direction of the tangent vectors (3xN)
%   delta    : Angles between the flow and the faces [rad] (1xN)
%   cp       : Face pressure coefficient (1xN)
%   ctau     : Face shear coefficient (1xN)
%   cd       : Face drag coefficient (1xN)
%   cl       : Face lift coefficient (1xN)
%   Cf_w     : Body force coefficient in wind axes (3x1)
%   Cf_f     : Body force coefficient in flight axes (3x1)
%   Cm_B     : Body moment coefficient in body axes (3x1)
%   Aref     : Reference area used in calculations [m^2]
%   AreaProj : Projected area to the flow [m^2]
%   Lref     : Refrence length used for calcualtions [m]
%   param_eq : Structure containing input GSI parameters
%
% where N = number of faces
%
% Inputs:
%       fiName   : Name of the .mat file containing the meshdata structure
%                   XData
%                   YData
%                   ZData
%                   MatID
%                   Areas
%                   SurfN
%                   BariC
%                   Lref
%       aoaS     : Angle(s) of attack [rad]
%       aosS     : Angle(s) of sideslip [rad]
%       eqmodel  : String containing the name of the equation to be used to calculate coefficients
%       param_eq : Parameters associated to "eqmodel"
%       flag_shad: Flag to perform shadow analysis (1=perform, 0=no)
%       flag_sol : Flag to calculate solar wind coefficients
%       pathin   : Path for loop runs (only for loop sims)
%
% Outputs:
%       fileOut  : Output file
%
% Author: David Mostaza-Prieto
% The University of Manchester
% November 2012
%
%------------- BEGIN CODE --------------

% Load mesh parameters
load([fiName,'.mat']);
x = meshdata.XData;
y = meshdata.YData;
z = meshdata.ZData;
areas = meshdata.Areas;
surfN = meshdata.SurfN;
barC = meshdata.BariC;
Lref = meshdata.Lref;
matID = meshdata.MatID;

indexAoA = length(aoaS);
indexAoS = length(aosS);

%Waitbar
if verb
    h = waitbar(0,'Please wait...','Name','Calculating coefficients...');
end

% Create output folder if required
ADBSat_path = ADBSat_dynpath;
if (indexAoA*indexAoS) > 1
    res_path = fullfile(ADBSat_path,'inou','results');
    foldname = strcat(fiName,'_',datestr(now,30),'_',num2str(randi(1000)));
    mkdir(fullfile(res_path,foldname,filesep));
    pathsav = fullfile(res_path,foldname);
else
    pathsav = fullfile(ADBSat_path,'inou','results');
end

% Values to save in output
var_out = {'aoa';'aos';'tauDir';'delta';'cp';'ctau';'cd';'cl';...
    'Cf_w';'Cf_f';'Cm_B';'Cf_s';'Cm_S';'Aref';'AreaProj';'Lref';'param_eq'};
if flag_sol
    var_out = [var_out;{'Cf_s';'Cm_S'}];
end

for ii = 1:indexAoA
    for jj = 1:indexAoS
        
        aoa = aoaS(ii);
        aos = aosS(jj);
        
        L_wb = dcmbody2wind(aoa, aos); % Body to Wind
        L_gb = [1 0 0; 0 -1 0; 0 0 -1]; % Body to Geometric
        L_gw = L_gb*L_wb'; % Wind to Geometric
        L_fb = [-1 0 0; 0 1 0; 0 0 -1]; % Body to Flight
        
        % Flow direction
        vdir = L_gw * [-1;0;0];
        vdir = vdir/norm(vdir);
        
        % Surface normals and angles between faces and flow
        vMatrix = [vdir(1)*ones(1,length(surfN(1,:)));...
            vdir(2)*ones(1,length(surfN(1,:)));...
            vdir(3)*ones(1,length(surfN(1,:)))];
        
        % Angles between flow and normals
        delta = real(acos(dot(-vMatrix,surfN)));
        
        % Local flat plate coefficients
        [cp, ctau, cd, cl] = mainCoeff(param_eq, delta, matID);
        
        if flag_sol
            [cn, cs] = coeff_solar(delta, param_eq);
        end
        
        % Shadow analysis
        indB = find(delta*180/pi>90);
        areaB = areas;
        areaB(indB) = 0;
        
        if flag_shad
            [shadPan] = shadowAnaly(x, y, z, barC, delta, L_gw);
            
            cp(shadPan) = 0;
            ctau(shadPan) = 0;
            cd(shadPan) = 0;
            cl(shadPan) = 0;
            areaB(shadPan) = 0;
            
            cn(shadPan) = 0;
            cs(shadPan) = 0;
        end
        
        % Areas
        AreaProj = areaB*cos(delta)';
        AreaT = sum(areas);
        
        Aref = AreaT/2;
        
        % Shear direction
        tauDir = cross(surfN,cross(vMatrix,surfN)); % direction of the shear coefficient
        ModT = sqrt(dot(tauDir,tauDir));
        tauDir = tauDir./[ModT;ModT;ModT]; % normalise
        
        TF = isnan(tauDir(1,:)); % Check any non defined tauDir (angles 90,180)
        tauDir(:,TF) = 0;
        
        % Force coefficeints
        Cf_g = 1/(Aref).*(tauDir*(ctau'.*areas') - surfN*(cp'.*areas'));
        Cf_w = L_gw' * Cf_g;
        Cf_f = L_fb * L_gb' * Cf_g;
        
        % Moment coefficients
        Cmom_G = 1/(Aref*Lref).*(cross(barC,tauDir)*(ctau'.*areas') +...
            cross(barC,-surfN)*(cp'.*areas'));
        %Cmom_G_b = sum((1/(Aref*Lref).*(cross(barC, (tauDir.*(ctau.*areas) - surfN.*(cp.*areas))))),2);
        Cm_B = L_gb' * Cmom_G;
        
        % Solar coefficients
        if flag_sol
            cstau = cs.*sin(delta); % resolve incident to shear
            cs_n = cs.*cos(delta); % resolve incident to normal
            Cs_G = 1/(Aref).*(tauDir*(cstau'.*areas') - surfN*((cn+cs_n)'.*areas'));
            Cf_s = L_gw'*Cs_G;
            Cms_G = 1/(Aref*Lref).*(cross(barC,tauDir)*(cstau'.*areas') +...
                cross(barC,-surfN)*((cn+cs_n)'.*areas'));
            Cm_S = L_gb'*Cms_G;
        end
        
        % Save file
        fileOut = fullfile(pathsav,[fiName,'_a',mat2str(aoa*180/pi),'s',...
            mat2str(aos*180/pi),'.mat']);
        save(fileOut, var_out{:})
        
        if verb
            percent = (ii-1)*1/indexAoA+jj*1/(indexAoS*indexAoA);
            waitbar(percent,h,sprintf('%6.4g',percent*100));
        end
        
    end
end

if verb
    delete(h); % clean-up waitbar
end

if (indexAoA*indexAoS) > 1
    % Creates a merged aerodynamic database from multiple .mat files
    fileOut  = mergeAEDB(pathsav, fiName, del);
end

%------------- END OF CODE --------------