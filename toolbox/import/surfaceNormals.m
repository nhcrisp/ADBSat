function [ surfN, areas, bariC ] = surfaceNormals( x, y, z )
%SURFACENORMALS Calculates the surface normals, areas, and baricenters of a triangular mesh
%
% Inputs:   
%       x : X coordinates of the vertices mesh elements (3xN)
%       y : Y coordinates of the vertices mesh elements (3xN)
%       z : Z coordinates of the vertices mesh elements (3xN)
%
% Outputs:
%       surfN   : Normalised normal vectors (3xN) 
%       areas   : Aera of each triangle (1xN)
%       bariC   : Baricenter of each triangle (3xN)
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
  
xV1 = x(2,:)-x(1,:);
xV2 = x(3,:)-x(1,:);

yV1 = y(2,:)-y(1,:);
yV2 = y(3,:)-y(1,:);

zV1 = z(2,:)-z(1,:);
zV2 = z(3,:)-z(1,:);

V1 = [xV1;yV1;zV1];
V2 = [xV2;yV2;zV2];

% Baricenters
bariC = 1/3.*[x(1,:)'+x(2,:)'+x(3,:)',y(1,:)'+y(2,:)'+y(3,:)',z(1,:)'+ ...
    z(2,:)'+z(3,:)'];
bariC = bariC';

% Normals
surfN = cross(V1,V2);
ModN  = sqrt(dot(surfN,surfN));
surfN = surfN./[ModN;ModN;ModN]; %Normalise

% Areas
areas = 0.5.*ModN;

%------------- END OF CODE --------------
