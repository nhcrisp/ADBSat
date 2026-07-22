function [ matOut ] = ADBSatImport_app( modIn, pathOut, verb, meshParam, axs)
%ADBSATIMPORT_APP Creates a file .mat containing a structure with the following fields
%     XData     : 3xN matrix X coordinates of the vertices (triangular)
%     YData     : 3xN matrix Y coordinates of the vertices (triangular)
%     ZData     : 3xN matrix Z coordinates of the vertices (triangular)
%     Areas     : Areas of the triangular faces  
%     SurfN     : Surface normals 
%     BariC     : Surfaces baricenters 
%      Lref     : Reference longitude
%
% Inputs:
%       fileIn  : Input filename string should include either .stl or .obj extension
%       verb    : Verbose flag
%       axs     : (Optional) Struct of target UIAxes handles for inline app rendering.
%                 Recognised fields: mesh1 (original mesh), mesh2 (final/
%                 subdivided mesh or processed STEP mesh), quality (mesh
%                 quality map), normals (surface normals), material
%                 (material ID map). Any field left empty/absent falls
%                 back to opening a standalone figure.
%
% Outputs:
%       matOut : Full path to the output MAT file   
%
% Author: David Mostaza-Prieto
% The University of Manchester
% November 2012
%
%--- Copyright notice ---%
% Copyright (C) 2021 The University of Manchester
% Written by David Mostaza Prieto,  Nicholas H. Crisp, Luciana Sinpetru, 
% Sabrina Livadiotti and Joseph Tucker
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

if nargin < 5
    axs = [];
end

[modPath,modName,ext] = fileparts(modIn); % Path to the .obj or .stl file

% Overwrite incoming pathOut to enforce the fixed "inou/models" storage
% location. This is resolved dynamically via ADBSat_dynpath (the ADBSat
% toolbox root) rather than hardcoded as a path relative to MATLAB's
% current working folder (pwd) - a hardcoded relative string like
% 'ADBSat-Master/inou/models' only resolves correctly if pwd happens to
% be the exact parent folder of "ADBSat-Master", which will not be true
% in general (e.g. when launching the app from ADBSat-Master/GUI, or
% from any other working folder). Using ADBSat_dynpath keeps this
% consistent with how the rest of the toolkit locates its own root.
ADBSat_path = ADBSat_dynpath();
pathOut = fullfile(ADBSat_path, 'inou', 'models');

if ~exist(pathOut, 'dir')
    mkdir(pathOut);
end

% Create .obj file from .stl if required (using meshlabserver)
if strcmpi(ext,'.stl')
    [err] = stl2obj(modIn);
    if ~err
        objname = [modName,'.obj'];
        objpath = fullfile(modPath,objname); % Input:
        matOut = importobjtri_app(objpath, pathOut, modName, verb, meshParam, axs);
    end
elseif strcmpi(ext,'.obj')
    objname = [modName,'.obj'];
    objpath = fullfile(modPath,objname); % Input:
    matOut = importobjtri_app(objpath, pathOut, modName, verb, meshParam, axs);
elseif strcmpi(ext,'.stp') || strcmpi(ext,'.step')
    matOut = stpToMesh_app(modIn, pathOut, verb, meshParam, axs);
end

if verb
    % Route the normals/material-ID plots to their own dedicated axes when
    % running inside the app, instead of always popping standalone figures
    axNormals  = [];
    axMaterial = [];
    if isstruct(axs)
        if isfield(axs,'normals')  && ~isempty(axs.normals),  axNormals  = axs.normals;  end
        if isfield(axs,'material') && ~isempty(axs.material), axMaterial = axs.material; end
    end
    plotNormals_app(matOut, axNormals, axMaterial); % Plots the surface mesh with the normals
end

end