function [pathSav] = stpToMesh(modIn, pathOut, verb, meshParam)
% Imports a triangular mesh from a .stp file
%
% Inputs:
%       modIn      : input filepath
%       pathOut    : output filepath 
%       verb       : flag for command window output
%       elementSize: size of each element in the generated mesh
%
% Outputs:
%       pathSav    : output filepath
%
% Author: Joseph Tucker
% The University of Manchester
% July 2026
%
%--- Copyright notice ---%
% Copyright (C) 2026 The University of Manchester
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

elementSize = meshParam.elementSize;

model = createpde();

[modPath,modName,ext] = fileparts(modIn);
fileIn = fullfile(modPath,[modName,ext]); 

gm = importGeometry(model,modIn);

% This is a really ugly plot but it is useful to verify that the geometry
% is being imported correctly
if verb
    figure('Color', 'w', 'Name', 'Imported Geometry');
    pdegplot(model, 'FaceLabels', 'on', 'FaceAlpha', 0.5, 'Lighting','on');

    axis equal;
    view(3);
    grid on;             % Standard readability grid
    box off;             % Drop ugly outer bounding frame lines

    % Apply your explicit axis line thickness and minor ticks
    ax = gca;
    ax.XMinorTick = 'on';
    ax.YMinorTick = 'on';
    ax.ZMinorTick = 'on';
    ax.LineWidth = 1.5;   % Exact line thickness matching plotMeshQuality

    % Exact typography match using LaTeX and Times New Roman
    title('Geometry Imported from .stp File', 'Interpreter', 'latex', 'FontName', 'Times New Roman', 'FontSize', 16);
    xlabel('X Position', 'Interpreter', 'latex', 'FontName', 'Times New Roman', 'FontSize', 14);
    ylabel('Y Position', 'Interpreter', 'latex', 'FontName', 'Times New Roman', 'FontSize', 14);
    zlabel('Z Position', 'Interpreter', 'latex', 'FontName', 'Times New Roman', 'FontSize', 14);

end

% Element size is a fixed scalar and is completely independent of the size
% of the input geometry
mesh = generateMesh(model, 'Hmax', elementSize, 'GeometricOrder', 'linear'); 

% The generateMesh function will produce a 3D tetrahedral mesh. This is not
% suitable for our purposes and creates problems down the line as we try to
% produce a 2D triangular mesh using 4 coordinates per face. This solves
% that problem by only extracting face coordinates for faces that belong to
% one element only
F = boundaryFacets(mesh); 

% Get the full list of 3D node coordinates (X, Y, Z)
V = mesh.Nodes; 

% Auto Centering
if isfield(meshParam, 'centering') && meshParam.centering
    minCorner = min(V, [], 2);
    maxCorner = max(V, [], 2);
    
    bbCenter = 0.5 * (minCorner + maxCorner);
    V = V - bbCenter;
    
    if verb
        fprintf('Bounding Box Bounds:\n');
        fprintf('  X: [%.3f, %.3f] -> Center: %.3f\n', minCorner(1), maxCorner(1), bbCenter(1));
        fprintf('  Y: [%.3f, %.3f] -> Center: %.3f\n', minCorner(2), maxCorner(2), bbCenter(2));
        fprintf('  Z: [%.3f, %.3f] -> Center: %.3f\n', minCorner(3), maxCorner(3), bbCenter(3));
        fprintf('Mesh automatically centered at [0, 0, 0].\n');
    end
end

% Mesh Transformations
if isfield(meshParam, 'rotationAngles') && ~isempty(meshParam.rotationAngles)
    % Convert degrees to radians (Input vector stays: [Pitch, Roll, Yaw])
    rad = deg2rad(meshParam.rotationAngles);
    
    cx = cos(rad(2)); sx = sin(rad(2)); % Roll  (X-axis) -> rad(2)
    cy = cos(rad(1)); sy = sin(rad(1)); % Pitch (Y-axis) -> rad(1)
    cxz = cos(rad(3)); sz = sin(rad(3)); % Yaw  (Z-axis) -> rad(3)

    % Create individual 3D rotation matrices
    Rx = [1,  0,   0; 0, cx, -sx; 0, sx,  cx];
    Ry = [cy, 0,  sy; 0,  1,   0; -sy, 0, cy];
    Rz = [cxz, -sz, 0; sz, cxz,  0; 0,   0,  1];

    % Combined rotation matrix (Z * Y * X order)
    R = Rz * Ry * Rx;

    % Apply rotation. Since V is 3 x N, we multiply R * V directly
    V = R * V;
end

% Check if translation vector is provided in meshParam [dx, dy, dz]
if isfield(meshParam, 'translation') && ~isempty(meshParam.translation)
    % Ensure translation is a 3x1 column vector for implicit expansion with V (3 x N)
    transVec = meshParam.translation(:); 
    V = V + transVec;
end

% ========================================================================
% STEP APPEARANCE (COLOR) PARSER - ENTITY GRAPH BASED
% ========================================================================
%
% DISCLAIMER
% This section was produced using AI, based on code I had written or
% repurposed from other functions. I cannot claim all of it to be mine or
% even really claim to understand each and every line but so far it has
% produced correct results so I am not questioning it yet
%
% Colors in a STEP file are attached to specific ADVANCED_FACE entities
% through a chain of style entities (STYLED_ITEM /
% OVER_RIDING_STYLED_ITEM -> PRESENTATION_STYLE_ASSIGNMENT ->
% SURFACE_STYLE_USAGE -> SURFACE_STYLE_FILL_AREA -> FILL_AREA_STYLE ->
% FILL_AREA_STYLE_COLOUR -> COLOUR_RGB). This section walks that chain
% explicitly and maps each resolved color to the PDE Toolbox face ID it
% actually belongs to, instead of re-deriving "faces" from mesh triangle
% geometry and cycling through colors arbitrarily.
if verb
    fprintf('Parsing STEP file for explicit face-color relationships...\n');
end

% Read the raw STEP text database
fid = fopen(fileIn, 'r');
if fid == -1
    error('Could not open the STEP file. Verify fileIn path: %s', fileIn);
end
stepText = fread(fid, '*char')';
fclose(fid);

% Extract Triangle Vertices Matrix (3 x numFacets for X, Y, Z coordinates).
% This process has been carried across from previous versions of this that
% were designed for obj or stl
X = [V(1, F(1, :)); V(1, F(2, :)); V(1, F(3, :))];
Y = [V(2, F(1, :)); V(2, F(2, :)); V(2, F(3, :))];
Z = [V(3, F(1, :)); V(3, F(2, :)); V(3, F(3, :))];
numFacets = size(F, 2);

% Mesh Quality check and visualization. Identical to .obj process
A = V(:, F(1, :))';
B = V(:, F(2, :))';
C = V(:, F(3, :))';

E1 = B - A; l1_sq = sum(E1.^2, 2);
E2 = C - B; l2_sq = sum(E2.^2, 2);
E3 = A - C; l3_sq = sum(E3.^2, 2);

cross_prods = cross(E1, -E3, 2);
areas = 0.5 * sqrt(sum(cross_prods.^2, 2));

quality = (4 * sqrt(3) * areas) ./ (l1_sq + l2_sq + l3_sq);

if verb
    plotMeshQuality(F, V, quality, 0.5);
end

% This replaces the surfaceNormals function by doing it inside the script.
% This is deemed acceptable as this only needs to be performed once per run
% of the code and therefore does not need to be repeated anywhere else
v1 = [X(2,:)-X(1,:); Y(2,:)-Y(1,:); Z(2,:)-Z(1,:)];
v2 = [X(3,:)-X(1,:); Y(3,:)-Y(1,:); Z(3,:)-Z(1,:)];
crossProd = [v1(2,:).*v2(3,:) - v1(3,:).*v2(2,:); ...
             v1(3,:).*v2(1,:) - v1(1,:).*v2(3,:); ...
             v1(1,:).*v2(2,:) - v1(2,:).*v2(1,:)];
crossMag = sqrt(sum(crossProd.^2, 1));
areas = 0.5 * crossMag;
surfN = crossProd ./ [crossMag; crossMag; crossMag];
bariC = [mean(X, 1); mean(Y, 1); mean(Z, 1)];
Lref = max(max(X)) - min(min(X));

% Tokenize every STEP entity: #ID = TYPE(args);
% If i'm honest i did not write this. .STP files are not easily readable
% directly by a human in the same way that .obj files are. This menat that
% I did not really know what to look for so U handed that off to a LLM to
% determine the best approach.
entTokens = regexp(stepText, '#(\d+)\s*=\s*([A-Z0-9_]+)\(([^;]*)\)\s*;', 'tokens');
entityType = containers.Map('KeyType','double','ValueType','char');
entityArgs = containers.Map('KeyType','double','ValueType','char');
for k = 1:numel(entTokens)
    id = str2double(entTokens{k}{1});
    entityType(id) = entTokens{k}{2};
    entityArgs(id) = entTokens{k}{3};
end

% Get the STEP face order the way MATLAB's importer traverses it: the
% order faces are listed inside each CLOSED_SHELL. Each shell's own
% entity ID is tracked alongside its faces so that, for assemblies
% with multiple bodies, each body's un-styled faces fall back to
% THAT body's default color rather than a single global fallback.
shellTokens = regexp(stepText, '#(\d+)\s*=\s*CLOSED_SHELL\(''[^'']*'',\(([^)]*)\)\)', 'tokens');
stepFaceOrder = [];
stepFaceShell = [];
for k = 1:numel(shellTokens)
    shellID = str2double(shellTokens{k}{1});
    ids = regexp(shellTokens{k}{2}, '#(\d+)', 'tokens');
    faceIDs = cellfun(@(c) str2double(c{1}), ids);
    stepFaceOrder = [stepFaceOrder, faceIDs]; %#ok<AGROW>
    stepFaceShell = [stepFaceShell, repmat(shellID, 1, numel(faceIDs))]; %#ok<AGROW>
end

% if the parser reads a different number of faces to what the native matlab
% step interpreter produces, it will proceed anyway but throw a warning
% that it may not produce correct results
nPDEFaces = model.Geometry.NumFaces;
if numel(stepFaceOrder) ~= nPDEFaces
    warning(['Number of STEP ADVANCED_FACE entities (%d) does not match ' ...
             'the number of PDE geometry faces (%d). Color mapping may ' ...
             'be misaligned - check against pdegplot(...,''FaceLabels'',''on'').'], ...
             numel(stepFaceOrder), nPDEFaces);
end

% MANIFOLD_SOLID_BREP entities wrap a CLOSED_SHELL, and STYLED_ITEMs
% that set a whole body's default color usually target the BREP, not
% the shell directly. Map BREP id -> the shell id it wraps so those
% defaults can be attributed to the right body.
brepToShell = containers.Map('KeyType','double','ValueType','double');
brepIDs = keys(entityType);
for k = 1:numel(brepIDs)
    id = brepIDs{k};
    if strcmp(entityType(id), 'MANIFOLD_SOLID_BREP')
        refs = regexp(entityArgs(id), '#(\d+)', 'tokens');
        if ~isempty(refs)
            brepToShell(id) = str2double(refs{1}{1});
        end
    end
end

% Walk every STYLED_ITEM / OVER_RIDING_STYLED_ITEM in file order,
% resolving each one's color. Colors targeting an ADVANCED_FACE are
% per-face; colors targeting a BREP/shell/anything else are treated
% as that body's default (or the global default, if the target can't
% be tied to a specific body). Later entries overwrite earlier ones,
% matching how OVER_RIDING_STYLED_ITEM works in the STEP standard.
styleTokens = regexp(stepText, ...
    '#\d+\s*=\s*(?:OVER_RIDING_)?STYLED_ITEM\(''[^'']*'',\(([^)]*)\),#(\d+)', 'tokens');

faceColour        = containers.Map('KeyType','double','ValueType','any');
shellDefaultColour = containers.Map('KeyType','double','ValueType','any');
globalDefaultColour = [0.7, 0.7, 0.7]; % last-resort fallback

for k = 1:numel(styleTokens)
    styleRefs = regexp(styleTokens{k}{1}, '#(\d+)', 'tokens');
    targetID  = str2double(styleTokens{k}{2});
    if isempty(styleRefs)
        continue
    end
    rgb = resolveColour(str2double(styleRefs{1}{1}), entityType, entityArgs, 0);
    if isempty(rgb)
        continue
    end
    if ~isKey(entityType, targetID)
        continue
    end
    targetType = entityType(targetID);
    if strcmp(targetType, 'ADVANCED_FACE')
        faceColour(targetID) = rgb;
    elseif strcmp(targetType, 'MANIFOLD_SOLID_BREP') && isKey(brepToShell, targetID)
        shellDefaultColour(brepToShell(targetID)) = rgb;
    elseif strcmp(targetType, 'CLOSED_SHELL')
        shellDefaultColour(targetID) = rgb;
    else
        globalDefaultColour = rgb; % e.g. product/representation-level style
    end
end

% Build the RGB color for each PDE face ID (face color, else that
% body's default, else the global default), then collapse to a
% compact palette (uniqueColors) for plotting.
faceRGB = zeros(nPDEFaces, 3);
for i = 1:min(nPDEFaces, numel(stepFaceOrder))
    stepFaceID = stepFaceOrder(i);
    shellID    = stepFaceShell(i);
    if isKey(faceColour, stepFaceID)
        faceRGB(i,:) = faceColour(stepFaceID);
    elseif isKey(shellDefaultColour, shellID)
        faceRGB(i,:) = shellDefaultColour(shellID);
    else
        faceRGB(i,:) = globalDefaultColour;
    end
end
for i = (numel(stepFaceOrder)+1):nPDEFaces
    faceRGB(i,:) = globalDefaultColour; % only hit if the count mismatch warning fired above
end

[uniqueColors, ~, faceColorIdx] = unique(faceRGB, 'rows', 'stable');
numColors = size(uniqueColors, 1);

% Assign each mesh triangle to a PDE geometric face using the mesh's
% own region membership (not re-derived plane geometry), then color
% it according to that face's resolved color index.
% To keep consistency with ADBSat, the appearance and colour data is placed
% inside an M array corresponding to the material data of .obj files. This
% is not perfect as ideally it would be able to read material data and map
% that instead of appearances, however for now this is as much as i can
% achieve
M = zeros(1, numFacets);
for i = 1:nPDEFaces
    faceNodes = findNodes(mesh, 'region', 'Face', i);
    inFace = false(1, size(V,2));
    inFace(faceNodes) = true;
    triInFace = all(inFace(F), 1);
    M(triInFace) = faceColorIdx(i);
end

numPhysicalFaces = nPDEFaces;

if any(M == 0)
    warning('%d facet(s) could not be matched to any geometric face.', sum(M==0));
end

if verb
    fprintf('Resolved colors for %d STEP face(s) into %d distinct color(s).\n', ...
        nPDEFaces, numColors);
end
% ========================================================================

% Back to human-written code to place all required data for ADBSat into the
% meshdata struct and save it to the models folder in inou

meshdata.XData = X;
meshdata.YData = Y;
meshdata.ZData = Z;
meshdata.MatID = M; 
meshdata.Areas = areas;
meshdata.SurfN = surfN;
meshdata.BariC = bariC;
meshdata.Lref  = Lref;

if ~exist(pathOut, 'dir')
    if verb
        fprintf('Output folder does not exist yet, creating: %s\n', pathOut);
    end
    mkdir(pathOut);
end

pathSav = fullfile(pathOut,[modName,'.mat']);
try
    save(pathSav, 'meshdata')
catch saveErr
    error('Failed to save meshdata to %s\nReason: %s', pathSav, saveErr.message);
end
if exist(pathSav, 'file')
    if verb
        fprintf('Saved mesh data to: %s\n', pathSav);
    end
else
    warning('save() did not raise an error, but %s was not found afterward.', pathSav);
end
matOut = pathSav;

% Final quality control and debugging calculations
nfaces  = length(X(1,:));
Tarea   = sum(areas);
MaxArea = max(areas);
MinArea = min(areas);
nMat    = max(M);

if ~all(areas)
    warndlg('Zero area element detected. Check input mesh quality')
end

% final readouts for QoL
if verb
    fprintf('Import finished!\n');
    fprintf('******************************************\n');
    fprintf('Number of elements   : %u\n', nfaces);
    fprintf('Total area           : %f\n', Tarea);
    fprintf('Material/Color Refs  : %u\n', nMat);
    fprintf('******************************************\n');
    
    % Render final multi-colored surface mesh plot. Shows colours as they
    % are registered in the .stp file, with the mesh on top. Probably my
    % favourite plot in this script
    if verb
        figure('Color', 'w', 'Name', 'Processed Model');
        trisurf(F', V(1,:), V(2,:), V(3,:), M, ...
            'EdgeColor', [0.1, 0.1, 0.1], ...
            'FaceAlpha', 0.9);
        colormap(uniqueColors);

    % Geometry layout normalization
    axis equal;
    view(3);
    grid on;             % Standard readability grid
    box off;             % Drop ugly outer bounding frame lines

    % Apply your explicit axis line thickness and minor ticks
    ax = gca;
    ax.XMinorTick = 'on';
    ax.YMinorTick = 'on';
    ax.ZMinorTick = 'on';
    ax.LineWidth = 1.5;   % Exact line thickness matching plotMeshQuality

    % Exact typography match using LaTeX and Times New Roman
    title('Processed Mesh', 'Interpreter', 'latex', 'FontName', 'Times New Roman', 'FontSize', 16);
    xlabel('X Position', 'Interpreter', 'latex', 'FontName', 'Times New Roman', 'FontSize', 14);
    ylabel('Y Position', 'Interpreter', 'latex', 'FontName', 'Times New Roman', 'FontSize', 14);
    zlabel('Z Position', 'Interpreter', 'latex', 'FontName', 'Times New Roman', 'FontSize', 14);
    end
    % Produces the traditional ADBSat plots to verify that data has been
    % stored correctly and in a way that the rest of the software can
    % interpret
    % plotNormals(matOut); % Plots the surface mesh with the normals
end
end
%% Functions
function rgb = resolveColour(id, entityType, entityArgs, depth)
    % Follows a STEP style-entity reference chain (STYLED_ITEM ->
    % PRESENTATION_STYLE_ASSIGNMENT -> SURFACE_STYLE_USAGE ->
    % SURFACE_STYLE_FILL_AREA -> FILL_AREA_STYLE -> FILL_AREA_STYLE_COLOUR
    % -> COLOUR_RGB) down to the base color. Each of these entity types has
    % exactly one relevant outgoing reference, so simply following the
    % first "#<id>" found in an entity's argument list walks the whole
    % chain generically without hardcoding every entity type's argument
    % layout. Returns [] if the chain dead-ends (e.g. NULL_STYLE, which has
    % no numeric reference to follow) or exceeds the recursion depth guard.
    rgb = [];
    if depth > 12 || ~isKey(entityType, id)
        return
    end
    typ  = entityType(id);
    args = entityArgs(id);
    if strcmp(typ, 'COLOUR_RGB')
        nums = regexp(args, '[-\d.eE+]+', 'match');
        rgb = str2double(nums(end-2:end));
        return
    end
    if strcmp(typ, 'DRAUGHTING_PRE_DEFINED_COLOUR') || strcmp(typ, 'PRE_DEFINED_COLOUR')
        % Some exporters reference one of the 8 standard named colors
        % (ISO 10303-46) instead of an explicit RGB triple.
        nameTok = regexp(args, '''([^'']*)''', 'tokens');
        if ~isempty(nameTok)
            rgb = predefinedColourLookup(nameTok{1}{1});
        end
        return
    end
    refs = regexp(args, '#(\d+)', 'tokens');
    if isempty(refs)
        return % e.g. PRESENTATION_STYLE_ASSIGNMENT((NULL_STYLE(.NULL.)))
    end
    rgb = resolveColour(str2double(refs{1}{1}), entityType, entityArgs, depth+1);
end

function rgb = predefinedColourLookup(name)
    % The 8 standard colors defined by ISO 10303-46 for
    % DRAUGHTING_PRE_DEFINED_COLOUR / PRE_DEFINED_COLOUR entities.
    switch lower(name)
        case 'black',   rgb = [0 0 0];
        case 'red',     rgb = [1 0 0];
        case 'green',   rgb = [0 1 0];
        case 'blue',    rgb = [0 0 1];
        case 'yellow',  rgb = [1 1 0];
        case 'magenta', rgb = [1 0 1];
        case 'cyan',    rgb = [0 1 1];
        case 'white',   rgb = [1 1 1];
        otherwise,      rgb = [];
    end
end