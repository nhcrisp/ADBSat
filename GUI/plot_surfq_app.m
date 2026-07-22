function h = plot_surfq_app(fileIn, modIn, aoa_deg, aos_deg, param, ax)
% PLOT_SURFQ_APP Plots the surface mesh with color proportional to the chosen parameter
%
% Inputs:
%    fileIn     : Name of the file containing the results (fiName_eqmodel)
%    modIn      : Folder containig the file
%    aoa_deg    : Angle of attack [deg]
%    aos_deg    : Angle of sideslip [deg]
%    param      : Surface parameter to plot (cp, ctau, cd, cl)
%    ax         : (Optional) Target UIAxes handle for inline app rendering

[~,modName,~] = fileparts(modIn);
load(modIn);
x = meshdata.XData;
y = meshdata.YData;
z = meshdata.ZData;

s = load(fileIn);
if isfield(s, 'aedb')
    disp('Please select a single ADBSat output .mat file')
end

aoa = deg2rad(aoa_deg);
aos = deg2rad(aos_deg);

L_wb = [cos(aos)*cos(aoa), sin(aos), sin(aoa)*cos(aos);...
    -sin(aos)*cos(aoa), cos(aos), -sin(aoa)*sin(aos);...
    -sin(aoa), 0, cos(aoa)]; 

L_gb = [1 0 0; 0 -1 0; 0 0 -1]; 
L_gw = L_gb*L_wb'; 
L_fb = [-1 0 0; 0 1 0; 0 0 -1]; 

ax_F = -L_fb * L_gb';
ax_W = -L_gw';

axlength = max([max(max(x))-min(min(x)), max(max(y))-min(min(y)), max(max(z))-min(min(z))]);

x0 = [0;0;0]; y0 = [0;0;0]; z0 = [0;0;0];

if nargin < 6 || isempty(ax)
    hFig = figure;
    ax = gca;
else
    cla(ax);
    hFig = [];
end

hold(ax, 'on');
% Enforce parent handling bounds strictly inside quiver annotations calls
W = quiver3(ax, 0,0,0,ax_W(1,1),ax_W(1,2),ax_W(1,3),axlength, 'b', 'LineWidth',2);
B = quiver3(ax, x0,y0,z0,L_gb(:,1),L_gb(:,2),L_gb(:,3),axlength,'r');
G = quiver3(ax, x0,y0,z0,[1;0;0],[0;1;0],[0;0;1],axlength,'g');
F = quiver3(ax, 0,0,0,ax_F(1,1),ax_F(1,2),ax_F(1,3),axlength,'k','LineWidth',2);

h = patch(x, y, z, s.(param), 'Parent', ax, 'EdgeColor', 'none');
colorbar(ax);
colormap(ax, 'cool');
legend(ax, [W,F,B,G],'Wind Vector','Flight Vector','Body Axes','Geometric Axes','Location','NorthWest')

string1 = strcat(param,' Surface Distribution');
string2 = strcat('AoA: ',sprintf('%.2f',aoa_deg),' deg,  AoS: ', sprintf('%.2f',aos_deg), ' deg');
xlabel(ax, 'X'); ylabel(ax, 'Y'); zlabel(ax, 'Z')
title(ax, char(string1,string2))

axis(ax, 'equal');
grid(ax, 'on');
hold(ax, 'off');

if ~isempty(hFig)
    dcm = datacursormode(hFig);
    set(dcm,'UpdateFcn',{@myupdatefcn,s.(param),param});
end
end

function txt = myupdatefcn(~,evt,data,name)
pos = get(evt,'Position');
ind = ceil(get(evt, 'DataIndex')/3);
txt = { sprintf('(x,y,z): (%g, %g, %g)', pos(1:3)),...
    sprintf('index: %g', ind),...
    sprintf('%s value: %g', name, data(ind))
    };
end