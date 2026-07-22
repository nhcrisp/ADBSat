function [pathSav] = importobjtri_app(fileIn, pathOut, struName, verb, meshParam, axs)
% IMPORTOBJTRI_APP Imports a triangular mesh from a .obj file
%
% Inputs:
%       fileIn     : input filepath
%       pathOut    : output filepath (Overridden dynamically inside script execution context)
%       struName   : output name for model file 
%       verb       : flag for command window output
%       axs        : (Optional) Struct of target UIAxes handles for inline app
%                    rendering. Fields: mesh1 (pre-subdivision mesh), mesh2
%                    (post-subdivision/final mesh), quality (mesh quality map)
%
% Outputs:
%       pathSav    : output file path
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

if nargin < 6
    axs = [];
end

axMesh1   = [];
axMesh2   = [];
axQuality = [];
if isstruct(axs)
    if isfield(axs,'mesh1')   && ~isempty(axs.mesh1),   axMesh1   = axs.mesh1;   end
    if isfield(axs,'mesh2')   && ~isempty(axs.mesh2),   axMesh2   = axs.mesh2;   end
    if isfield(axs,'quality') && ~isempty(axs.quality), axQuality = axs.quality; end
end

% Overwrite pathOut variable parameter references cleanly - resolved
% dynamically via ADBSat_dynpath rather than hardcoded relative to
% MATLAB's current working folder (pwd). A literal
% 'ADBSat-Master/inou/models' string is only reachable when pwd happens
% to be the exact parent of "ADBSat-Master" - it silently fails (or
% points at the wrong place) from any other working folder, even though
% the target directory genuinely exists on disk.
pathOut = fullfile(ADBSat_dynpath(), 'inou', 'models');

scaling = meshParam.subdivisions;

[V,F,M] = obj_fileTri2patch(fileIn);

% Auto Centering
if isfield(meshParam, 'centering') && meshParam.centering
    minCorner = min(V, [], 1);
    maxCorner = max(V, [], 1);
    
    % Calculate the exact center point of this bounding box
    bbCenter = 0.5 * (minCorner + maxCorner);
    
    % Translate the mesh so that this center point sits exactly at [0, 0, 0]
    V = V - bbCenter;
    
    if verb
        fprintf('Bounding Box Bounds:\n');
        fprintf('  X: [%.3f, %.3f] -> Center: %.3f\n', minCorner(1), maxCorner(1), bbCenter(1));
        fprintf('  Y: [%.3f, %.3f] -> Center: %.3f\n', minCorner(2), maxCorner(2), bbCenter(2));
        fprintf('  Z: [%.3f, %.3f] -> Center: %.3f\n', minCorner(3), maxCorner(3), bbCenter(3));
        fprintf('Mesh automatically centered at [0, 0, 0].\n');
    end
end

% Transformations
if isfield(meshParam, 'rotationAngles') && ~isempty(meshParam.rotationAngles)
    % Convert degrees to radians
    rad = deg2rad(meshParam.rotationAngles);
    cx = cos(rad(2)); sx = sin(rad(2)); % Roll (X-axis)
    cy = cos(rad(1)); sy = sin(rad(1)); % Pitch (Y-axis)
    cz = cos(rad(3)); sz = sin(rad(3)); % Yaw (Z-axis)

    % Create individual 3D rotation matrices
    Rx = [1,  0,   0; 0, cx, -sx; 0, sx,  cx];
    Ry = [cy, 0,  sy; 0,  1,   0; -sy, 0, cy];
    Rz = [cz, -sz, 0; sz, cz,  0; 0,   0,  1];

    % Combined rotation matrix (Z * Y * X order)
    R = Rz * Ry * Rx;

    % Apply rotation to all vertices
    V = V * R';
end

if isfield(meshParam, 'translation') && ~isempty(meshParam.translation)
    % Apply translation: shift each row by [dx, dy, dz]
    V = V + meshParam.translation; 
end

numFaces = size(F,1);

% Mesh Upscaling
if scaling >= 1 && scaling <= 5
    if verb
        plotMesh_app(V, F, 0, axMesh1)
    end
    
    % Loop through each subdivision level requested
    for s = 1:scaling
        % Identify all unique edges to place new midpoint vertices
        edgesSorted = sort([F(:,1), F(:,2); F(:,2), F(:,3); F(:,3), F(:,1)], 2);
        [uniqueEdges, ~, edgeMap] = unique(edgesSorted, 'rows');
        
        % Compute coordinates of midpoints for each unique edge
        midpoints = 0.5 * (V(uniqueEdges(:,1), :) + V(uniqueEdges(:,2), :));
        
        % Append new vertices onto our vertex array
        newVertStartIdx = size(V, 1);
        V = [V; midpoints]; %#ok<AGROW>
        
        % Map back to retrieve the global index of each midpoint
        m1 = newVertStartIdx + edgeMap(1:numFaces);          % Midpoint of edge 1-2
        m2 = newVertStartIdx + edgeMap(numFaces+1:2*numFaces);% Midpoint of edge 2-3
        m3 = newVertStartIdx + edgeMap(2*numFaces+1:end);    % Midpoint of edge 3-1
        
        % Each triangle is split into 4 new triangles
        f1 = [F(:,1), m1, m3];
        f2 = [m1, F(:,2), m2];
        f3 = [m3, m2, F(:,3)];
        f4 = [m1, m2, m3];
        
        F = [f1; f2; f3; f4];
        numFaces = size(F, 1);
    end
    if verb
        plotMesh_app(V, F, 1, axMesh2)
    end

elseif scaling == 0
    if verb
        fprintf('No mesh upscaling was performed\n');
        plotMesh_app(V, F, 0, axMesh1)
        plotMesh_app(V, F, 0, axMesh2)
    end

elseif scaling > 5
    answer = questdlg(['Scaling of greater than 5 requires significant computational resources.' ...
        ' Are you sure you want to continue?'], ...
    	'Warning', ...
    	'Yes, fry my processor please','nonononono take me back','nonononono take me back');
    switch answer
        case 'Yes, fry my processor please'
            if verb
                plotMesh_app(V, F, 0, axMesh1)
            end
            % Loop through each subdivision level requested
            for s = 1:scaling
                % Identify all unique edges to place new midpoint vertices
                edgesSorted = sort([F(:,1), F(:,2); F(:,2), F(:,3); F(:,3), F(:,1)], 2);
                [uniqueEdges, ~, edgeMap] = unique(edgesSorted, 'rows');

                % Compute coordinates of midpoints for each unique edge
                midpoints = 0.5 * (V(uniqueEdges(:,1), :) + V(uniqueEdges(:,2), :));

                % Append new vertices onto our vertex array
                newVertStartIdx = size(V, 1);
                V = [V; midpoints]; %#ok<AGROW>

                % Map back to retrieve the global index of each midpoint
                m1 = newVertStartIdx + edgeMap(1:numFaces);          % Midpoint of edge 1-2
                m2 = newVertStartIdx + edgeMap(numFaces+1:2*numFaces);% Midpoint of edge 2-3
                m3 = newVertStartIdx + edgeMap(2*numFaces+1:end);    % Midpoint of edge 3-1

                % Each triangle is split into 4 new triangles
                f1 = [F(:,1), m1, m3];
                f2 = [m1, F(:,2), m2];
                f3 = [m3, m2, F(:,3)];
                f4 = [m1, m2, m3];

                F = [f1; f2; f3; f4];
                numFaces = size(F, 1);
            end
            if verb
                plotMesh_app(V, F, 1, axMesh2)
            end
        case 'nonononono take me back'
            error('ADBSat:importobjtri_app:SubdivisionCancelled', ...
                'Mesh import cancelled: subdivision level %d exceeds the recommended maximum of 5 and was not confirmed.', scaling);
        otherwise
            % Covers the dialog being closed directly (answer returned as
            % empty) - treat the same as a decline rather than silently
            % continuing with the un-subdivided mesh.
            error('ADBSat:importobjtri_app:SubdivisionCancelled', ...
                'Mesh import cancelled: no response given to the subdivision confirmation dialog.');
    end
else
    error("Scaling value cannot be negative. Recommended values are between 0 and 5")
end

% Scale materials configuration array matrices maps cleanly
M = repelem(M, 1, 4^scaling);

FI = F';
X = [V(FI(1,:),1)'; V(FI(2,:),1)'; V(FI(3,:),1)'];
Y = [V(FI(1,:),2)'; V(FI(2,:),2)'; V(FI(3,:),2)'];
Z = [V(FI(1,:),3)'; V(FI(2,:),3)'; V(FI(3,:),3)'];

A = V(F(:,1), :);
B = V(F(:,2), :);
C = V(F(:,3), :);

E1 = B - A; l1_sq = sum(E1.^2, 2);
E2 = C - B; l2_sq = sum(E2.^2, 2);
E3 = A - C; l3_sq = sum(E3.^2, 2);

cross_prods = cross(E1, -E3, 2);
areas = 0.5 * sqrt(sum(cross_prods.^2, 2));

quality = (4 * sqrt(3) * areas) ./ (l1_sq + l2_sq + l3_sq);

if verb
    plotMeshQuality_app(F, V, quality, 0.5, axQuality);
end

[surfN, areas, bariC] = surfaceNormals(X, Y, Z);
Lref = max(max(X))-min(min(X));

meshdata.XData = X;
meshdata.YData = Y;
meshdata.ZData = Z;
meshdata.MatID = M;
meshdata.Areas = areas;
meshdata.SurfN = surfN;
meshdata.BariC = bariC;
meshdata.Lref  = Lref;

if ~exist(pathOut, 'dir')
    mkdir(pathOut);
end

pathSav = fullfile(pathOut,[struName,'.mat']);
save(pathSav, 'meshdata')

nfaces  = length(X(1,:));
Tarea   = sum(areas);
MaxArea = max(areas);
MinArea = min(areas);
nMat = max(M);

if ~all(areas)
    warndlg('Zero area element detected. Check input mesh quality')
end

if verb
    fprintf('Import finished!\n');
    fprintf('******************************************\n');
    fprintf('Number of elements: %u\n', nfaces);
    fprintf('Total area  : %f\n', Tarea);
    fprintf('Maximum element area: %f\n', MaxArea);
    fprintf('Minumum element area: %f\n', MinArea);
    fprintf('Reference length (maxX-minX): %f\n', Lref);
    fprintf('Number of material references: %u\n', nMat);
    fprintf('Created file ''%s'' \n', pathSav);
end

end