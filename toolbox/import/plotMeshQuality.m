function plotMeshQuality(F, V, quality, threshold)
%% Fig 1

if size(F, 2) ~= 3, F = F'; end
if size(V, 2) ~= 3, V = V'; end
if size(quality, 2) ~= 1, quality = quality'; end

figure('Color', 'w', 'Name', '3D Mesh Quality Map');

h = trisurf(F, V(:,1), V(:,2), V(:,3), quality, ...
    'EdgeColor', 'none', ...
    'FaceLighting', 'gouraud');

% Geometric presentation adjustments
axis equal;
grid on;
view(3);
camlight headlit; 
material dull;

% Setup actionable color grading
colormap(turbo);
clim([0 1]);

% Add a descriptive colorbar
cb = colorbar;
cb.Label.Interpreter = 'latex';
cb.Label.String = 'Element Quality ($0 = $ Degenerate, $1 = $ Equilateral)';
cb.Label.FontSize = 12;

title('Spatial Element Quality Distribution', 'Interpreter', 'latex', 'FontSize', 14);

%% Fig 2
if size(F, 2) ~= 3, F = F'; end
if size(V, 2) ~= 3, V = V'; end

% 1. Find indices of faces failing the threshold
badFaceIdx = quality < threshold;
F_bad = F(badFaceIdx, :);
quality_bad = quality(badFaceIdx);

if isempty(F_bad)
    fprintf('Congratulations! No elements found below threshold of %0.2f\n', threshold);
    return;
end

figure('Color', 'w', 'Name', 'Problem Mesh Elements');

% Plot the ghost silhouette of the full mesh in faint transparent gray 
% so you maintain spatial context of where you are on the satellite
trisurf(F, V(:,1), V(:,2), V(:,3), ...
    'FaceColor', [0.9 0.9 0.9], ...
    'EdgeColor', 'none', ...
    'FaceAlpha', 0.15); 
hold on;

% Overlay the problematic elements in bright solid red with visible edges
trisurf(F_bad, V(:,1), V(:,2), V(:,3), ...
    'FaceColor', 'r', ...
    'EdgeColor', [0.4 0 0], ...
    'LineWidth', 1);

axis equal; view(3); grid on; camlight;
title(sprintf('Isolating Bad Elements (Quality $< %0.2f$, Count: %d)', ...
    threshold, size(F_bad,1)), 'Interpreter', 'latex', 'FontSize', 13);
