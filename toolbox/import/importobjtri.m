function [ pathSav ] = importobjtri( fileIn, struName, verb )
%IMPORTOBJTRI Imports a triangular mesh from a .obj file
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