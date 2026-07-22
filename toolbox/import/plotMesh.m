function plotMesh(V, F, divided)

if ~divided
    figure('Color', 'w' ,'Name', 'Imported Mesh'); % Enforce custom white background layout
    colororder("gem12");  % Apply your specified toolkit color palette

    % Render the base triangular mesh
    trisurf(F, V(:,1), V(:,2), V(:,3), ...
        'FaceColor', 'w', ...
        'EdgeColor', 'b', ...
        'LineWidth', 0.5);

    % Geometry layout normalization
    axis equal;
    view(3);
    grid on;             % Standard readability grid
    box off;             % Drop ugly outer bounding frame lines

    % Apply explicit axis line thickness and minor ticks
    ax = gca;
    ax.XMinorTick = 'on';
    ax.YMinorTick = 'on';
    ax.ZMinorTick = 'on';
    ax.LineWidth = 1.5;   % Exact line thickness matching plotMeshQuality

    % Exact typography match using LaTeX and Times New Roman
    title('Original Mesh Geometry', 'Interpreter', 'latex', 'FontName', 'Times New Roman', 'FontSize', 14);
    xlabel('X Position', 'Interpreter', 'latex', 'FontName', 'Times New Roman', 'FontSize', 14);
    ylabel('Y Position', 'Interpreter', 'latex', 'FontName', 'Times New Roman', 'FontSize', 14);
    zlabel('Z Position', 'Interpreter', 'latex', 'FontName', 'Times New Roman', 'FontSize', 14);

else
    figure('Color', 'w', 'Name', 'Subdivided Mesh'); % Enforce custom white background layout
    colororder("gem12");  % Apply specified toolkit color palette

    % Render the base triangular mesh
    trisurf(F, V(:,1), V(:,2), V(:,3), ...
        'FaceColor', 'w', ...
        'EdgeColor', 'b', ...
        'LineWidth', 0.5);

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
    title('Subdivided Mesh Geometry', 'Interpreter', 'latex', 'FontName', 'Times New Roman', 'FontSize', 14);
    xlabel('X Position', 'Interpreter', 'latex', 'FontName', 'Times New Roman', 'FontSize', 14);
    ylabel('Y Position', 'Interpreter', 'latex', 'FontName', 'Times New Roman', 'FontSize', 14);
    zlabel('Z Position', 'Interpreter', 'latex', 'FontName', 'Times New Roman', 'FontSize', 14);
end