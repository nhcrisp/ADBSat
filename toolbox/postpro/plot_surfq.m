function h = plot_surfq(fileIn, modIn, aoa_deg, aos_deg, param)
% Plots the surface mesh with color proportional to the chosen parameter
%
% Inputs:
%    fileIn     : Name of the file containing the results (fiName_eqmodel)
%    modIn      : Folder containig the file
%    aoa_deg    : Angle of attack [deg]
%    aos_deg    : Angle of sideslip [deg]
%    param      : Surface parameter to plot (cp, ctau, cd, cl)
%
% Outputs:
%   h           : A patch object, containing the shape colour-coded by the chosen parameter
%
% Author: David Mostaza-Prieto
% The University of Manchester
% December 2012
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

% Load model mesh
[~,modName,~] = fileparts(modIn);
load(modIn);
x = meshdata.XData;
y = meshdata.YData;
z = meshdata.ZData;

% Load results for indicated aoa, aos
s = load(fileIn);
if isfield(s, 'aedb')
    disp('Please select a single ADBSat output .mat file')
end

% Convert to Radians
aoa = deg2rad(aoa_deg);
aos = deg2rad(aos_deg);

%L_wb = dcmbody2wind(aoa, aos); % Body to Wind (Aerospace Toolbox)
L_wb = [cos(aos)*cos(aoa), sin(aos), sin(aoa)*cos(aos);...
    -sin(aos)*cos(aoa), cos(aos), -sin(aoa)*sin(aos);...
    -sin(aoa), 0, cos(aoa)]; % Body to Wind

L_gb = [1 0 0; 0 -1 0; 0 0 -1]; % Body to Geometric
L_gw = L_gb*L_wb'; % Wind to Geometric
L_fb = [-1 0 0; 0 1 0; 0 0 -1]; % Body to Flight

ax_F = -L_fb * L_gb';
ax_W = -L_gw';

axlength = max([max(max(x))-min(min(x)), max(max(y))-min(min(y)), max(max(z))-min(min(z))]);

% Reference Point
x0 = [0;0;0]; y0 = [0;0;0]; z0 = [0;0;0];

hFig = figure;
hold on
%Wind
%W = quiver3(x0,y0,z0,L_gw(:,1),L_gw(:,2),L_gw(:,3),axlength, 'b');
W = quiver3(0,0,0,ax_W(1,1),ax_W(1,2),ax_W(1,3),axlength, 'b', 'LineWidth',2);
% Body
B = quiver3(x0,y0,z0,L_gb(:,1),L_gb(:,2),L_gb(:,3),axlength,'r');
%quiver3(0,0,0,L_gb(1,1),L_gb(1,2),L_gb(1,3),axlength, 'r', 'LineWidth',2)
% Geometric
G = quiver3(x0,y0,z0,[1;0;0],[0;1;0],[0;0;1],axlength,'g');
% Flight
%F = quiver3(x0,y0,z0,L_fb(:,1),L_fb(:,2),L_fb(:,3),axlength,'k');
F = quiver3(0,0,0,ax_F(1,1),ax_F(1,2),ax_F(1,3),axlength,'k','LineWidth',2);
axis equal
grid on

h = patch(x, y, z, s.(param));
colorbar
legend([W,F,B,G],'Wind Vector','Flight Vector','Body Axes','Geometric Axes','Location','NorthWest')
% set(h,'EdgeAlpha',0)
string1 = strcat(param,' Surface Distribution');
string2 = strcat('AoA: ',sprintf('%.2f',aoa_deg),' deg,  AoS: ', sprintf('%.2f',aos_deg), ' deg');
xlabel('X'); ylabel('Y'); zlabel('Z')
title(char(string1,string2))

axis equal

% Set custom update function
dcm = datacursormode(hFig);
set(dcm,'UpdateFcn',{@myupdatefcn,s.(param),param});
end

function txt = myupdatefcn(~,evt,data,name)
pos = get(evt,'Position');
ind = ceil(get(evt, 'DataIndex')/3);
txt = { sprintf('(x,y,z): (%g, %g, %g)', pos(1:3)),...
    sprintf('index: %g', ind),...
    sprintf('%s value: %g', name, data(ind))
    };
end

%------------- END OF CODE --------------
