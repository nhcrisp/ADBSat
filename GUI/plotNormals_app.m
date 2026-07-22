function [] = plotNormals_app( fiName, axNormals, axMat )
% PLOTNORMALS_APP Plots mesh normals and Material IDs with enabled interaction

load(fiName, 'meshdata');
x = meshdata.XData;
y = meshdata.YData;
z = meshdata.ZData;
barC = meshdata.BariC;
surfN = meshdata.SurfN;
matID = meshdata.MatID;

mats = numel(unique(matID));
map = lines(mats);

%% --- PLOT 1: SURFACE NORMALS ---
if nargin < 2 || isempty(axNormals)
    figure('Name','ADBSat Mesh', 'Color', 'w');
    axNormals = gca;
else
    cla(axNormals);
end

hold(axNormals, 'on');
quiver3(axNormals, barC(1,:), barC(2,:), barC(3,:), surfN(1,:), surfN(2,:), surfN(3,:), 0.5);
patch(x, y, z, matID, 'Parent', axNormals, 'EdgeColor', 'none'); 

xlabel(axNormals, 'X'); ylabel(axNormals, 'Y'); zlabel(axNormals, 'Z');
axis(axNormals, 'equal');
axis(axNormals, 'tight');
grid(axNormals, 'on');
view(axNormals, 3);

if isprop(axNormals, 'Interactions')
    enableDefaultInteractivity(axNormals);
end

hold(axNormals, 'off');

%% --- PLOT 2: MATERIAL ID ---
if nargin < 3 || isempty(axMat)
    figure('Name','ADBSat Material ID', 'Color', 'w');
    axMat = gca;
else
    cla(axMat);
end

colormap(axMat, map);
P = patch(x, y, z, matID, 'Parent', axMat);
P.FaceAlpha = 0.7;
P.LineStyle = 'none';

grid(axMat, 'on');
xlabel(axMat, 'X'); ylabel(axMat, 'Y'); zlabel(axMat, 'Z');
axis(axMat, 'equal');
axis(axMat, 'tight');
view(axMat, 3);

if isprop(axMat, 'Interactions')
    enableDefaultInteractivity(axMat);
end

cb = colorbar(axMat);
set(cb, 'Ticks', 0:1:mats);
end