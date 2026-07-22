function plotMeshQuality_app(F, V, quality, threshold, ax)


if size(F, 2) ~= 3, F = F'; end
if size(V, 2) ~= 3, V = V'; end
if size(quality, 2) ~= 1, quality = quality'; end

if nargin < 5 || isempty(ax)
    figure('Color', 'w', 'Name', '3D Mesh Quality Map');
    ax = gca;
else
    cla(ax);
end

trisurf(F, V(:,1), V(:,2), V(:,3), quality, ...
    'EdgeColor', 'none', ...
    'FaceLighting', 'none', ... 
    'Parent', ax);

axis(ax, 'equal');
grid(ax, 'on');
view(ax, 3);

colormap(ax, turbo);
clim(ax, [0 1]);

% Generate colorbar object bound to the active axes context
cb = colorbar(ax);

set(cb.Label, 'String', 'Element Quality ($0 = $ Degenerate, $1 = $ Equilateral)', ...
              'Interpreter', 'latex', ...
              'FontSize', 12);
% ----------------------------------------

title(ax, 'Spatial Element Quality Distribution', 'Interpreter', 'latex', 'FontSize', 14);


if nargin < 5 || isempty(ax)
    badFaceIdx = quality < threshold;
    F_bad = F(badFaceIdx, :);

    if isempty(F_bad)
        fprintf('Congratulations! No elements found below threshold of %0.2f\n', threshold);
        return;
    end

    figure('Color', 'w', 'Name', 'Problem Mesh Elements');

    trisurf(F, V(:,1), V(:,2), V(:,3), ...
        'FaceColor', [0.9 0.9 0.9], ...
        'EdgeColor', 'none', ...
        'FaceAlpha', 0.15); 
    hold on;

    trisurf(F_bad, V(:,1), V(:,2), V(:,3), ...
        'FaceColor', 'r', ...
        'EdgeColor', [0.4 0 0], ...
        'LineWidth', 1);

    axis equal; view(3); grid on;
    title(sprintf('Isolating Bad Elements (Quality $< %0.2f$, Count: %d)', ...
        threshold, size(F_bad,1)), 'Interpreter', 'latex', 'FontSize', 13);
end

end