function [V,F,X,Y,Z,M] = obj_fileTri2patch(fileIn)
% Reads a .obj trianglar mesh file and outputs the vertex coordinates
%
% Inputs:
%   fileIn : filename string of .obj filepath
%
% Ouputs:
%    V  : Vertices
%    F  : Faces
%    X  : X coordinates of the vertices mesh elements (3xN)
%    Y  : Y coordinates of the vertices mesh elements (3xN)
%    Z  : Z coordinates of the vertices mesh elements (3xN)
%    M  : Material identifier of each mesh element (1xN)
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

V = zeros(0,3);
F = zeros(0,3);
M = zeros(0,1);
vertex_index = 1;
face_index = 1;
mat_id = 0;
fid = fopen(fileIn,'r');
line = fgets(fid);

while ischar(line)
    vertex = sscanf(line,'v %f %f %f');
    face = sscanf(line,'f %d %d %d');
    face_long = sscanf(line,'f %d/%d/%d %d/%d/%d %d/%d/%d');
    face_long2 = sscanf(line,'f %d/%d %d/%d %d/%d');
    face_long3 = sscanf(line,'f %d//%d %d//%d %d//%d');
    material = sscanf(line,'usemtl %d;');
    
    % see if line is vertex command if so add to vertices
    if(size(vertex)>0)
        V(vertex_index,:) = vertex;
        vertex_index = vertex_index+1;
        
        % see if line is a material identifier
    elseif(size(material)>0)
        mat_id = material;
        
        % see if line is simple face command if so add to faces
    elseif(size(face)>0)
        
        facelength   = [length(face),length(face_long),length(face_long2),...
            length(face_long3)];
        [~,indFfor] = max(facelength);
        
        
        switch indFfor
            case 1
                F(face_index,:) = face;
                M(face_index,:) = mat_id;
                % see if line is a long face command if so add to faces
            case 2
                % remove normal and texture indices
                face_long = face_long(1:3:end);
                F(face_index,:) = face_long;
            case 3
                face_long = face_long(1:3:end);
                F(face_index,:) = face_long;
            case 4
                face_long = face_long3(1:2:end);
                F(face_index,:) = face_long;
        end
        M(face_index,:) = mat_id;
        face_index = face_index+1;
    else
        %       fprintf('Ignored: %s',line);
    end
    
    line = fgets(fid);
end
fclose(fid);

% Convert to patch structure

FI = F';
X = [V(FI(1,:),1)';V(FI(2,:),1)';V(FI(3,:),1)'];
Y = [V(FI(1,:),2)';V(FI(2,:),2)';V(FI(3,:),2)'];
Z = [V(FI(1,:),3)';V(FI(2,:),3)';V(FI(3,:),3)'];
M = M';

%------------- END OF CODE --------------
