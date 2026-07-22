function [pathSav] = stpToMesh_app(modIn, pathOut, verb, meshParam, axs)
% STPTOMESH_APP Imports a triangular mesh from a .stp file
%
% Inputs:
%       modIn      : input filepath
%       pathOut    : output filepath (Overridden dynamically inside script execution context)
%       verb       : flag for command window output
%       axs        : (Optional) Struct of target UIAxes handles for inline app
%                    rendering. Fields: mesh1 (imported geometry preview),
%                    mesh2 (final processed mesh), quality (mesh quality map)
%
% Outputs:
%       pathSav    : output filepath

if nargin < 5
    axs = [];
end

axGeom    = [];
axMesh2   = [];
axQuality = [];
if isstruct(axs)
    if isfield(axs,'mesh1')   && ~isempty(axs.mesh1),   axGeom    = axs.mesh1;   end
    if isfield(axs,'mesh2')   && ~isempty(axs.mesh2),   axMesh2   = axs.mesh2;   end
    if isfield(axs,'quality') && ~isempty(axs.quality), axQuality = axs.quality; end
end

pathOut = fullfile(ADBSat_dynpath(), 'inou', 'models');

elementSize = meshParam.elementSize;

model = createpde();

[modPath,modName,ext] = fileparts(modIn);
fileIn = fullfile(modPath,[modName,ext]); 

gm = importGeometry(model,modIn);

if verb
    if isempty(axGeom)
        figure('Color', 'w', 'Name', 'Imported Geometry');
        targetAx = gca;
    else
        targetAx = axGeom;
        cla(targetAx);
    end
    
    % Force the geometry preview plotting straight into target container handle context
    pdegplot(model, 'FaceLabels', 'on', 'FaceAlpha', 0.5, 'Lighting','on', 'Parent', targetAx);

    axis(targetAx, 'equal');
    view(targetAx, 3);
    grid(targetAx, 'on');             
    box(targetAx, 'off');             

    targetAx.XMinorTick = 'on';
    targetAx.YMinorTick = 'on';
    targetAx.ZMinorTick = 'on';
    targetAx.LineWidth = 1.5;   

    title(targetAx, 'Geometry Imported from .stp File', 'Interpreter', 'latex', 'FontName', 'Times New Roman', 'FontSize', 16);
    xlabel(targetAx, 'X Position', 'Interpreter', 'latex', 'FontName', 'Times New Roman', 'FontSize', 14);
    ylabel(targetAx, 'Y Position', 'Interpreter', 'latex', 'FontName', 'Times New Roman', 'FontSize', 14);
    zlabel(targetAx, 'Z Position', 'Interpreter', 'latex', 'FontName', 'Times New Roman', 'FontSize', 14);
end

mesh = generateMesh(model, 'Hmax', elementSize, 'GeometricOrder', 'linear'); 

F = boundaryFacets(mesh); 
V = mesh.Nodes; 

% Auto Centering
if isfield(meshParam, 'centering') && meshParam.centering
    minCorner = min(V, [], 2);
    maxCorner = max(V, [], 2);
    bbCenter = 0.5 * (minCorner + maxCorner);
    V = V - bbCenter;
end

% Mesh Transformations
if isfield(meshParam, 'rotationAngles') && ~isempty(meshParam.rotationAngles)
    rad = deg2rad(meshParam.rotationAngles);
    cx = cos(rad(2)); sx = sin(rad(2)); % Roll (X-axis)
    cy = cos(rad(1)); sy = sin(rad(1)); % Pitch (Y-axis)
    cxz = cos(rad(3)); sz = sin(rad(3)); % Yaw (Z-axis)

    Rx = [1,  0,   0; 0, cx, -sx; 0, sx,  cx];
    Ry = [cy, 0,  sy; 0,  1,   0; -sy, 0, cy];
    Rz = [cxz, -sz, 0; sz, cxz,  0; 0,   0,  1];

    R = Rz * Ry * Rx;
    V = R * V;
end

if isfield(meshParam, 'translation') && ~isempty(meshParam.translation)
    transVec = meshParam.translation(:); 
    V = V + transVec;
end

if verb
    fprintf('Parsing STEP file for explicit face-color relationships...\n');
end

fid = fopen(fileIn, 'r');
if fid == -1
    error('Could not open the STEP file. Verify fileIn path: %s', fileIn);
end
stepText = fread(fid, '*char')';
fclose(fid);

X = [V(1, F(1, :)); V(1, F(2, :)); V(1, F(3, :))];
Y = [V(2, F(1, :)); V(2, F(2, :)); V(2, F(3, :))];
Z = [V(3, F(1, :)); V(3, F(2, :)); V(3, F(3, :))];
numFacets = size(F, 2);

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
    plotMeshQuality_app(F, V, quality, 0.5, axQuality);
end

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

entTokens = regexp(stepText, '#(\d+)\s*=\s*([A-Z0-9_]+)\(([^;]*)\)\s*;', 'tokens');
entityType = containers.Map('KeyType','double','ValueType','char');
entityArgs = containers.Map('KeyType','double','ValueType','char');
for k = 1:numel(entTokens)
    id = str2double(entTokens{k}{1});
    entityType(id) = entTokens{k}{2};
    entityArgs(id) = entTokens{k}{3};
end

shellTokens = regexp(stepText, '#(\d+)\s*=\s*CLOSED_SHELL\(''[^'']*'',\(([^)]*)\)\)', 'tokens');
stepFaceOrder = [];
stepFaceShell = [];
for k = 1:numel(shellTokens)
    shellID = str2double(shellTokens{k}{1});
    ids = regexp(shellTokens{k}{2}, '#(\d+)', 'tokens');
    faceIDs = cellfun(@(c) str2double(c{1}), ids);
    stepFaceOrder = [stepFaceOrder, faceIDs]; 
    stepFaceShell = [stepFaceShell, repmat(shellID, 1, numel(faceIDs))]; 
end

nPDEFaces = model.Geometry.NumFaces;

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

styleTokens = regexp(stepText, ...
    '#\d+\s*=\s*(?:OVER_RIDING_)?STYLED_ITEM\(''[^'']*'',\(([^)]*)\),#(\d+)', 'tokens');

faceColour        = containers.Map('KeyType','double','ValueType','any');
shellDefaultColour = containers.Map('KeyType','double','ValueType','any');
globalDefaultColour = [0.7, 0.7, 0.7]; 

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
        globalDefaultColour = rgb; 
    end
end

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
    faceRGB(i,:) = globalDefaultColour; 
end

[uniqueColors, ~, faceColorIdx] = unique(faceRGB, 'rows', 'stable');

M = zeros(1, numFacets);
for i = 1:nPDEFaces
    faceNodes = findNodes(mesh, 'region', 'Face', i);
    inFace = false(1, size(V,2));
    inFace(faceNodes) = true;
    triInFace = all(inFace(F), 1);
    M(triInFace) = faceColorIdx(i);
end

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

pathSav = fullfile(pathOut,[modName,'.mat']);
save(pathSav, 'meshdata');
matOut = pathSav;

nfaces  = length(X(1,:));
Tarea   = sum(areas);
nMat    = max(M);

if ~all(areas)
    warndlg('Zero area element detected. Check input mesh quality')
end

if verb
    fprintf('Import finished!\n');
    fprintf('******************************************\n');
    fprintf('Number of elements   : %u\n', nfaces);
    fprintf('Total area           : %f\n', Tarea);
    fprintf('Material/Color Refs  : %u\n', nMat);
    fprintf('******************************************\n');
    
    if isempty(axMesh2)
        figure('Color', 'w', 'Name', 'Processed Model');
        finalAx = gca;
    else
        finalAx = axMesh2;
        cla(finalAx);
    end
    
    % Enforce drawing strictly into app layout view pane boundaries via the Parent property
    trisurf(F', V(1,:), V(2,:), V(3,:), M, ...
        'EdgeColor', [0.1, 0.1, 0.1], ...
        'FaceAlpha', 0.9, ...
        'Parent', finalAx);
    colormap(finalAx, uniqueColors);

    axis(finalAx, 'equal');
    view(finalAx, 3);
    grid(finalAx, 'on');             
    box(finalAx, 'off');             

    finalAx.XMinorTick = 'on';
    finalAx.YMinorTick = 'on';
    finalAx.ZMinorTick = 'on';
    finalAx.LineWidth = 1.5;   

    title(finalAx, 'Processed Mesh', 'Interpreter', 'latex', 'FontName', 'Times New Roman', 'FontSize', 16);
    xlabel(finalAx, 'X Position', 'Interpreter', 'latex', 'FontName', 'Times New Roman', 'FontSize', 14);
    ylabel(finalAx, 'Y Position', 'Interpreter', 'latex', 'FontName', 'Times New Roman', 'FontSize', 14);
    zlabel(finalAx, 'Z Position', 'Interpreter', 'latex', 'FontName', 'Times New Roman', 'FontSize', 14);
end
end

function rgb = resolveColour(id, entityType, entityArgs, depth)
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
        nameTok = regexp(args, '''([^'']*)''', 'tokens');
        if ~isempty(nameTok)
            rgb = predefinedColourLookup(nameTok{1}{1});
        end
        return
    end
    refs = regexp(args, '#(\d+)', 'tokens');
    if isempty(refs)
        return 
    end
    rgb = resolveColour(str2double(refs{1}{1}), entityType, entityArgs, depth+1);
end

function rgb = predefinedColourLookup(name)
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