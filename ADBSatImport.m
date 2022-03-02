function [ pathSav ] = ADBSatImport( modIn, pathOut, verb )
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
%------------- BEGIN CODE --------------

[modPath,modName,ext] = fileparts(modIn); % Path to the .obj or .stl file

% Create .obj file from .stl if required (using meshlabserver)
if strcmpi(ext,'.stl')
    [err] = stl2obj(modName);
    if ~err
        objname = [modName,'.obj'];
    end
else
    objname = [modName,'.obj'];
end

objpath = fullfile(modPath,objname); % Input:

pathSav = importobjtri(objpath, pathOut, modName, verb);

if verb
    plotNormals(pathSav); % Plots the surface mesh with the normals
end
