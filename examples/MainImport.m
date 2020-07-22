% Creates a file .mat in "/inou/models" containing a structure with the following fields:

%     XData: 3xN matrix X coordinates of the vertices (triangular)
%     YData: 3xN matrix Y coordinates of the vertices (triangular)
%     ZData: 3xN matrix Z coordinates of the vertices (triangular)
%     Areas: Areas of the triangular faces  
%     SurfN: Surface normals 
%     BariC: Surfaces baricenters 
%      Lref: Refrence longitude 

filename   = 'cube.obj'; % Input: Name of the file in /inou/obj_files

[~,name,ext] = fileparts(filename);

if strcmpi(ext,'.stl')
    [status] = stl2obj(name);
    filename = [name,'.obj'];
end

objname = filename; % Input: Name of the "obj" file in /inou/obj_files
[~,name,~] = fileparts(filename);
struName = name;     % Output: Name of the matlab file mesh in /inou/models

importobjtri(filename,struName)

h = plotNormals(struName); % Plots the surface mesh with the normals