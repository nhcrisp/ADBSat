function [] = plotNormals( fiName )
%PLOTNORMALS Plots the mesh and normal vectors
%
% Inputs:   
%    fiName : Path to the .mat file containing the meshdata structure
%
% Author: David Mostaza-Prieto
% The University of Manchester
% September 2012
%
%------------- BEGIN CODE --------------

load(fiName);
x = meshdata.XData;
y = meshdata.YData;
z = meshdata.ZData;
barC = meshdata.BariC;
surfN = meshdata.SurfN;
matID = meshdata.MatID;

figure 
quiver3(barC(1,:),barC(2,:),barC(3,:),surfN(1,:),surfN(2,:),surfN(3,:))
hold on

patch(x, y, z, matID);

axis equal

%------------- END OF CODE --------------