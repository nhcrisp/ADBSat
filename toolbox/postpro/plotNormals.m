function [] = plotNormals( fiName )
% Plots the mesh and normal vectors
%
% Inputs:   
%    fiName : Path to the .mat file containing the meshdata structure
%
% Author: David Mostaza-Prieto
% The University of Manchester
% September 2012
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

load(fiName);
x = meshdata.XData;
y = meshdata.YData;
z = meshdata.ZData;
barC = meshdata.BariC;
surfN = meshdata.SurfN;
matID = meshdata.MatID;

mats = numel(unique(matID));
map = lines(mats);

figure('Name','ADBSat Mesh')
quiver3(barC(1,:),barC(2,:),barC(3,:),surfN(1,:),surfN(2,:),surfN(3,:))
hold on
patch(x, y, z, matID);
xlabel('X')
ylabel('Y')
zlabel('Z')
axis equal
axis vis3d

figure('Name','ADBSat Material ID')
colormap(map);
P = patch(x, y, z, matID);
P.FaceAlpha = 0.7;
P.LineStyle = 'none';
grid on;
xlabel('X')
ylabel('Y')
zlabel('Z')
cb = colorbar;
set(cb,'Ticks',0:1:mats)
axis equal
axis vis3d

%------------- END OF CODE --------------
