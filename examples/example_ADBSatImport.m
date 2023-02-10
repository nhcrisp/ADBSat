% Creates a file .mat in "/inou/models" containing a structure with the following fields:
%     XData: 3xN matrix X coordinates of the vertices (triangular)
%     YData: 3xN matrix Y coordinates of the vertices (triangular)
%     ZData: 3xN matrix Z coordinates of the vertices (triangular)
%     Areas: Areas of the triangular faces  
%     SurfN: Surface normals 
%     BariC: Surfaces baricenters 
%      Lref: Refrence longitude 
%
%--- Copyright notice ---%
% Copyright (C) 2021 The University of Manchester
% Written by David Mostaza Prieto,  Nicholas H. Crisp, Luciana Sinpetru and Sabrina Livadiotti
%
% This file is part of the ADBSat toolkit.
%
% This program is free software: you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or (at
% your option) any later version.
%
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
% Public License for more details.
%
% You should have received a copy of the GNU General Public License along
% with this program. If not, see <http://www.gnu.org/licenses/>.
%----------- BEGIN CODE ------------%

ADBSat_path = ADBSat_dynpath;

filename = 'cube.obj'; % Input: Name of the file in /inou/obj_files

modIn = fullfile(ADBSat_path,'inou','obj_files',filename);
[modPath,modName,ext] = fileparts(filename);

% Create .obj file from .stl if required (using meshlabserver)
if strcmpi(ext,'.stl')
    [err] = stl2obj(modIn);
    if ~err
        objname = [modName,'.obj'];
    end
else
    objname = [modName,'.obj'];
end

objpath = fullfile(modPath,objname);

pathOut = fullfile(ADBSat_path,'inou','models');

matOut = importobjtri(modIn, pathOut, modName, verb);

plotNormals(matOut); % Plots the surface mesh with the normals

%------------END CODE------------%
