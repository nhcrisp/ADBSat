function [pathSav] = importobjtri(fileIn, pathOut, struName, verb, meshParam)
% Imports a triangular mesh from a .obj file
%
% Inputs:
%       fileIn     : input filepath
%       struName   : output name for model file 
%       verb       : flag for command window output
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
    % Convert degrees to radians (Input vector stays: [Pitch, Roll, Yaw])
    rad = deg2rad(meshParam.rotationAngles);
    
    cx = cos(rad(2)); sx = sin(rad(2)); % Roll  (X-axis) -> rad(2)
    cy = cos(rad(1)); sy = sin(rad(1)); % Pitch (Y-axis) -> rad(1)
    cz = cos(rad(3)); sz = sin(rad(3)); % Yaw   (Z-axis) -> rad(3)

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
numVertices = size(V,1);

% Mesh Validation
results = validateMesh(V,F);

orientable = results.orientable;
watertight = results.watertight;
edgeManifold = results.edgeManifold;
vertexManifold = results.vertexManifold;
hasUnreferencedVertices = results.hasUnreferencedVertices;
hasDegenerateFaces = results.hasDegenerateFaces;

if orientable && watertight
    if verb
        fprintf('Mesh is orientable and watertight. Proceeding with upscaling...\n');
    end
elseif ~orientable
    error("Not all face normals point outwards. Check input mesh quality");
elseif ~watertight
    if ~edgeManifold
        error("Mesh is not edge manifold. Check input mesh quality. If using SolidWorks, consider switching to" + ...
            " blender or other software that can handle .obj files natively for better results")
    elseif ~vertexManifold
        error("Mesh is not vertex manifold. Check input mesh qualityIf using SolidWorks, consider switching to" + ...
            " blender or other software that can handle .obj files natively for better results");
    end
end

% Mesh Upscaling. This upscales the mesh by subdividing each element into 4
% smaller elements. This is repeated as many time as is specified by the 
% user. Some checks are employed to ensure that the scaling factor is
% feasible
if scaling >= 1 && scaling <= 5
    if verb
        divided = 0;
        plotMesh(V, F, divided)
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
        divided = 1;
        plotMesh(V, F, divided)
    end

elseif scaling == 0
    if verb
        fprintf('No mesh upscaling was performed\n'); 
    end

elseif scaling > 5

    answer = questdlg(['Scaling of greater than 5 requires significant computational resources.' ...
        ' Are you sure you want to continue?'], ...
    	'Warning', ...
    	'Yes, fry my processor please','nonononono take me back','nonononono take me back');
    switch answer

        case 'Yes, fry my processor please'

            mesh = surfaceMesh(V,F);
            if verb
                divided = 0;
                plotMesh(V, F, divided)
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
                V = [V; midpoints]; 

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
                divided = 1;
                plotMesh(V, F, divided)
            end

        case 'nonononono take me back'
            return
    end
else
    error("Scaling value cannot be negative. Recommended values are between 0 and 5")
end

% After this, it is returned to the arrays that it started with. These are
% now processed  by splitting the vertices into their x,y,z coordinates.
% The material matrix does not scale with the mesh, creating problems down
% the line when array sizes do not match. To solve this, the repelem
% function is used, which scales the matrix without changing the ordering
% of the numbers. note that any number^0 = 1
M = repelem(M, 1, 4^scaling);

FI = F';
X = [V(FI(1,:),1)'; V(FI(2,:),1)'; V(FI(3,:),1)'];
Y = [V(FI(1,:),2)'; V(FI(2,:),2)'; V(FI(3,:),2)'];
Z = [V(FI(1,:),3)'; V(FI(2,:),3)'; V(FI(3,:),3)'];

% Element quality checks. Now that the mesh has reached its final form, so
% to speak, the function will check that the mesh elements are not
% horrendously skewed or have poor aspect ratios. It plots this data on a
% graph so it is clear to the user whether or not their mesh is well
% constructed

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
    plotMeshQuality(F, V, quality, 0.5);
end

% the surfaceNormals function is given the x,y,z coords in order to find
% the surface normals, areas and barycentres of each element. This is
% important for the shadowing algorithm
[surfN, areas, bariC] = surfaceNormals(X, Y, Z);
Lref = max(max(X))-min(min(X));

% puts all relevant data for the model into a struct so that it can be
% easily saved to a .mat file and moved around different script and
% functions
meshdata.XData = X;
meshdata.YData = Y;
meshdata.ZData = Z;
meshdata.MatID = M;
meshdata.Areas = areas;
meshdata.SurfN = surfN;
meshdata.BariC = bariC;
meshdata.Lref  = Lref;

pathSav = fullfile(pathOut,[struName,'.mat']);

save(pathSav, 'meshdata')

% this is mostly QoL stuff. Outputs to the user so they can see mesh size
% and quality. Also useful for debugging if number of material references
% shown is incorrect for example
nfaces  = length(X(1,:));
Tarea   = sum(areas);
MaxArea = max(areas);
MinArea = min(areas);
nMat = max(M);

% Final QC for the mesh after all steps
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

%------------- END OF CODE --------------
