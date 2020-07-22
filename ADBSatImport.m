function [ modName, pathSav ] = ADBSatImport( fileIn, verb )
%ADBSATIMPORT Creates a file .mat in "/inou/models" containing a structure with the following fields
%     XData     : 3xN matrix X coordinates of the vertices (triangular)
%     YData     : 3xN matrix Y coordinates of the vertices (triangular)
%     ZData     : 3xN matrix Z coordinates of the vertices (triangular)
%     Areas     : Areas of the triangular faces  
%     SurfN     : Surface normals 
%     BariC     : Surfaces baricenters 
%      Lref     : Refrence longitude
%
% Inputs:
%       fileIn  : Input filename string should include either .stl or .obj extension
%       verb    : Verbose flag
%
% Outputs:
%       modName : Name of the output model
%       pathSav : Full path to the output file   
%
% Author: David Mostaza-Prieto
% The University of Manchester
% November 2012
%
%------------- BEGIN CODE --------------

[pathstr,modName,ext] = fileparts(fileIn); % Name of the .obj or .stl file in /inou/obj_files

% Create .obj file from .stl if required (using meshlabserver)
if strcmpi(ext,'.stl')
    [status] = stl2obj(modName);
    if status
        objname = [modName,'.obj'];
    end
else
    objname = [modName,'.obj'];
end

objpath = fullfile(pathstr,objname); % Input:

pathSav = importobjtri(objpath, modName, verb);

if verb
    plotNormals(pathSav); % Plots the surface mesh with the normals
end