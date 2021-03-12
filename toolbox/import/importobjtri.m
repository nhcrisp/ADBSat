function [ pathSav ] = importobjtri( fileIn, struName, verb )
% Imports a triangular mesh from a .obj file
%
% Inputs:
%       fileIn     : input filepath
%       struName   : output name for model file 
%       verb       : flag for command window output
%
% Outputs:
%       pathSav    : output file path
%
% Author: David Mostaza-Prieto
% The University of Manchester
% September 2012
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

[V,F,X,Y,Z,M] = obj_fileTri2patch(fileIn);
[surfN, areas, bariC] = surfaceNormals(X, Y, Z);

Lref = max(max(X))-min(min(X));

meshdata.XData = X;
meshdata.YData = Y;
meshdata.ZData = Z;
meshdata.MatID = M;
meshdata.Areas = areas;
meshdata.SurfN = surfN;
meshdata.BariC = bariC;
meshdata.Lref  = Lref;

ADBSat_path = ADBSat_dynpath;

pathSav = fullfile(ADBSat_path,'inou','models',[struName,'.mat']);

save(pathSav, 'meshdata')

nfaces  = length(X(1,:));
Tarea   = sum(areas);
MaxArea = max(areas);
MinArea = min(areas);
nMat = max(M);

if verb
    fprintf('Import finished!\n');
    fprintf('******************************************\n');
    fprintf('Number of elements: %u\n', nfaces);
    fprintf('Total area  : %f\n', Tarea);
    fprintf('Maximum element area: %f\n', MaxArea);
    fprintf('Minumum element area: %f\n', MinArea);
    fprintf('Reference length (maxX-minX): %f\n', Lref);
    fprintf('Number of material references: %u\n', nMat);
    fprintf('Created file ''%s.mat'' \n', pathSav);
end

%------------- END OF CODE --------------
