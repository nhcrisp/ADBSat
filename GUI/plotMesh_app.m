function plotMesh_app(V, F, divided, ax)

if nargin < 4 || isempty(ax)
    if ~divided
        figure('Color', 'w' ,'Name', 'Imported Mesh'); 
    else
        figure('Color', 'w', 'Name', 'Subdivided Mesh');
    end
    ax = gca;
else
    cla(ax);
end

colororder(ax, "gem12");

trisurf(F, V(:,1), V(:,2), V(:,3), ...
    'FaceColor', 'w', ...
    'EdgeColor', 'b', ...
    'LineWidth', 0.5, ...
    'Parent', ax);

axis(ax, 'equal');
view(ax, 3);
grid(ax, 'on');             
box(ax, 'off');             

ax.XMinorTick = 'on';
ax.YMinorTick = 'on';
ax.ZMinorTick = 'on';
ax.LineWidth = 1.5;   

if ~divided
    title(ax, 'Original Mesh Geometry', 'Interpreter', 'latex', 'FontName', 'Times New Roman', 'FontSize', 14);
else
    title(ax, 'Subdivided Mesh Geometry', 'Interpreter', 'latex', 'FontName', 'Times New Roman', 'FontSize', 14);
end

xlabel(ax, 'X Position', 'Interpreter', 'latex', 'FontName', 'Times New Roman', 'FontSize', 14);
ylabel(ax, 'Y Position', 'Interpreter', 'latex', 'FontName', 'Times New Roman', 'FontSize', 14);
zlabel(ax, 'Z Position', 'Interpreter', 'latex', 'FontName', 'Times New Roman', 'FontSize', 14);

end