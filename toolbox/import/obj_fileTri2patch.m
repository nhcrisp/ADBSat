function [V,F,M] = obj_fileTri2patch(fileIn)
% Reads a .obj trianglar mesh file and outputs the vertex coordinates
%
% Inputs:
%   fileIn : filename string of .obj filepath
%
% Ouputs:
%    V  : Vertices
%    F  : Faces
%    M  : Material identifier of each mesh element (1xN)
%
% Author: David Mostaza-Prieto
% The University of Manchester
% September 2012
%
%--- Copyright notice ---%
% Copyright (C) 2021 The University of Manchester
% Written by David Mostaza Prieto,  Nicholas H. Crisp, Luciana Sinpetru, 
% Sabrina Livadiotti and Joseph Tucker
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

% This opens the file. If the file location is broken or the file is
% corrupted in some way, an error will now be thrown
fid = fopen(fileIn,'r');
if fid == -1
    error('Could not open file: %s', fileIn);
end

% The new preallocation creates a massive array rather than a small one,
% removing the penalty that is associated with a dynamically changing array
% size due to constant allocation of memory, copying and deleting in the
% RAM
chunk_size = 50000;
V = zeros(chunk_size, 3);
F = zeros(chunk_size, 3);
M = zeros(chunk_size, 1);

vertex_index = 1;
face_index = 1;
mat_id = 0;

% While loop that runs through every line one by one in the file. fgetl is
% marginally faster than fgets as it ignores newline characters
while ~feof(fid)
    line = fgetl(fid);

    % Skip empty lines or comments instantly. This removing the time
    % wasting done by the parsing algorithm on irrelevant lines in the file
    if isempty(line) || line(1) == '#'
        continue;
    end

    % Similar to an if statement. This system examines the first character
    % only to establish what each row represents and THEN scans it. This
    % prevent each row being scanned multiple times
    switch line(1)
        case 'v'
            if line(2) == ' '

                % This grabs the vertex coords and adds them to our matrix
                vertex = sscanf(line(3:end), '%f %f %f');
                if length(vertex) >= 3

                    % if there are more than 50000 vertices, it
                    % automatically adds 50000 rows rather than 1 row at a
                    % time
                    if vertex_index > size(V,1)
                        V(end+chunk_size, 3) = 0;
                    end
                    V(vertex_index,:) = vertex(1:3);
                    vertex_index = vertex_index + 1;
                end
            end
        
        % this now repeats with the face lines. There are multiple
        % different for formats for face lines. Rather than scanning for
        % each format and then waiting for one of them to be true, it
        % counts the number of slashes and then automatically finds the
        % right format. It then retrieves the useful data.
        case 'f'
            if line(2) == ' '
                
                num_slashes = sum(line == '/');

                if num_slashes == 0       % f v1 v2 v3

                    face = sscanf(line(3:end), '%d %d %d');

                elseif num_slashes == 3   % f v1/vt1 v2/vt2 v3/vt3

                    face_data = sscanf(line(3:end), '%d/%d %d/%d %d/%d');
                    face = face_data(1:2:end);

                elseif num_slashes == 6

                    if contains(line, '//') % f v1//vn1 v2//vn2 v3//vn3

                        face_data = sscanf(line(3:end), '%d//%d %d//%d %d//%d');
                        face = face_data(1:2:end);
                    
                    else                    % f v1/vt1/vn1 v2/vt2/vn2 v3/vt3/vn3
                        
                        face_data = sscanf(line(3:end), '%d/%d/%d %d/%d/%d %d/%d/%d');
                        face = face_data(1:3:end);
                    
                    end
                
                else
                    
                    error(".obj file is unsupported. All elements must have 3 vertices with" + ...
                        " an equal amount of data for each vertex");
                
                end

                % This auto extends the face and material arrays in the
                % event that there are more lines than the current chunk
                % size
                if length(face) >= 3

                    if face_index > size(F,1)

                        F(end+chunk_size, 3) = 0;
                        M(end+chunk_size, 1) = 0;

                    end
                       
                    % this takes the current data and stores it to the
                    % master matrices and then furthers the loop
                    F(face_index, :) = face(1:3);
                    M(face_index, 1) = mat_id;
                    face_index = face_index + 1;

                end
            end

        % repeat the same process for mtl lines
        case 'u'

            if strncmp(line, 'usemtl ', 7)
                mat_id = mat_id + 1;

            end
    end
end

fclose(fid);

% This shrinks the arrays back down to their necessary sizes so they can be
% passed back to the importobjtri function
V = V(1:vertex_index-1, :);
F = F(1:face_index-1, :);
M = M(1:face_index-1, :);

% Vectorized patch conversion. This previously also calculated the x,y,z
% coordinates of each vertex however as part of the mesh processing
% improvements, this has been moved to the importobjtri

M = M';

%------------- END OF CODE --------------
