% Calculates local and global coefficients for the triangular mesh geometry.
% Produces one or several output .m files with the object characteristics.
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
%           aoa      : Angle of attack [rad]
%           aos      : Angle of sideslip [rad]
%           tauDir   : Direction of the tangent vectors (3xN)
%           delta    : Angles between the flow and the faces [rad] (1xN)
%           cp       : Face pressure coefficient (1xN)
%           ctau     : Face shear coefficient (1xN)
%           cd       : Face drag coefficient (1xN)
%           cl       : Face lift coefficient (1xN)
%           Cf_w     : Body force coefficient in wind axes (3x1)
%           Cf_f     : Body force coefficient in flight axes (3x1)
%           Cm_B     : Body moment coefficient in body axes (3x1)
%           AreaRef  : Reference area used in calculations [m^2]
%           AreaProj : Projected area to the flow [m^2]
%           LenRef   : Refrence length used for calcualtions [m]
%           param_eq : Structure containing input GSI parameters
%
% where N = number of faces.
%
% Author: David Mostaza-Prieto
% The University of Manchester
% November 2012
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

function [fileOut] = calc_coeff(fiName, respath, aoaS, aosS, param_eq, flag_shad, flag_sol, del, verb)

[~,matName,~] = fileparts(fiName);

% Load mesh parameters
load(fiName,'meshdata');
x = meshdata.XData;
y = meshdata.YData;
z = meshdata.ZData;
areas = meshdata.Areas;
surfN = meshdata.SurfN;
barC = meshdata.BariC;
LenRef = meshdata.Lref;
matID = meshdata.MatID;

indexAoA = length(aoaS);
indexAoS = length(aosS);

%Waitbar
if verb
    h = waitbar(0,'Please wait...','Name','Calculating coefficients...');
end

% Create output folder if required
if (indexAoA*indexAoS) > 1
    foldname = strcat(matName,'_',datestr(now,30),'_',num2str(randi(1000)));
    mkdir(fullfile(respath,filesep,foldname));
    pathsav = fullfile(respath,foldname);
    aedb = 1;
else
    pathsav = fullfile(respath);
    aedb = 0;
end

% Values to save in output
var_out = {'aoa';'aos';'tauDir';'delta';'cp';'ctau';'cd';'cl';...
    'Cf_w';'Cf_f';'Cm_B';'AreaRef';'AreaProj';'LenRef';'param_eq';'shadow'};
if flag_sol
    var_out = [var_out;{'Cf_s';'Cm_S'}];
end

for ii = 1:indexAoA
    for jj = 1:indexAoS
        
        aoa = aoaS(ii);
        aos = aosS(jj);
        
        %L_wb = dcmbody2wind(aoa, aos); % Body to Wind (Aerospace Toolbox)
        L_wb = [cos(aos)*cos(aoa), sin(aos), sin(aoa)*cos(aos);...
            -sin(aos)*cos(aoa), cos(aos), -sin(aoa)*sin(aos);...
            -sin(aoa), 0, cos(aoa)]; % Body to Wind
        
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
        
        uD = vMatrix; % Unit drag vector
        uL = -cross(cross(uD,surfN),uD)./vecnorm(cross(cross(uD,surfN),uD)); % Unit lift vector for each panel
        % Undefined uL panels (plates normal to flow)
        col = find(all(isnan(uL),1));
        uL(:,col) = -surfN(:,col);
        % Negative dot product of unit drag and lift vectors with surface normal vector (March 2019)
        param_eq.gamma = dot(-uD,surfN);
        param_eq.ell = dot(-uL,surfN);
        
        % Local flat plate coefficients
        [cp, ctau, cd, cl, param_eq] = mainCoeff(param_eq, delta, matID);
        
        if flag_sol
            [cn, cs] = coeff_solar(delta, param_eq);
        end
        
        % Backwards facing panels
        areaB = areas;
        areaB(delta*180/pi>90) = 0;
        
        shadow = zeros(size(areas));
        shadow(areaB == 0) = 1;
        
        % Shadow analysis
        if flag_shad
            [shadPan] = shadowAnaly(x, y, z, barC, delta, L_gw, surfN);
            
            cp(shadPan) = 0;
            ctau(shadPan) = 0;
            cd(shadPan) = 0;
            cl(shadPan) = 0;
            areaB(shadPan) = 0;
            
            cn(shadPan) = 0;
            cs(shadPan) = 0;            
            shadow(shadPan) = 0.5;
        end
  
        % Areas
        AreaProj = areaB*cos(delta)'; % Projected
        AreaT = sum(areas); % Total
        AreaRef = AreaT/2; % ADBSat Reference
        
        % Shear direction
        tauDir = cross(surfN,cross(vMatrix,surfN)); % direction of the shear coefficient
        ModT = sqrt(dot(tauDir,tauDir));
        tauDir = tauDir./[ModT;ModT;ModT]; % normalise
        
        TF = isnan(tauDir(1,:)); % Check any non defined tauDir (angles 90,180)
        tauDir(:,TF) = 0;
        
        % Force coefficeints
        Cf_g = 1/(AreaRef).*(tauDir*(ctau'.*areas') - surfN*(cp'.*areas'));
        Cf_w = L_gw' * Cf_g;
        Cf_f = L_fb * L_gb' * Cf_g;
        
        % Moment coefficients
        Cmom_G = 1/(AreaRef*LenRef).*(cross(barC,tauDir)*(ctau'.*areas') +...
            cross(barC,-surfN)*(cp'.*areas'));
        %Cmom_G_b = sum((1/(AreaRef*LenRef).*(cross(barC, (tauDir.*(ctau.*areas) - surfN.*(cp.*areas))))),2);
        Cm_B = L_gb' * Cmom_G;
        
        % Solar coefficients
        if flag_sol
            cstau = cs.*sin(delta); % resolve incident to shear
            cs_n = cs.*cos(delta); % resolve incident to normal
            Cs_G = 1/(AreaRef).*(tauDir*(cstau'.*areas') - surfN*((cn+cs_n)'.*areas'));
            Cf_s = L_gw'*Cs_G;
            Cms_G = 1/(AreaRef*LenRef).*(cross(barC,tauDir)*(cstau'.*areas') +...
                cross(barC,-surfN)*((cn+cs_n)'.*areas'));
            Cm_S = L_gb'*Cms_G;
        end
        
        % Save file
        if aedb
            fileOut = fullfile(pathsav,[matName,'_a',mat2str(aoa*180/pi),'_s',...
                mat2str(aos*180/pi),'.mat']);
        else
            fileOut = [pathsav,'.mat'];
        end
        
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
    fileOut  = mergeAEDB(pathsav, matName, del);
end

%------------- END OF CODE --------------
