function [cp, ctau, cd, cl] = mainCoeff(param_eq, delta, matID)
% Passes the input parameters to the correct GSI model
%
% Inputs:
%       delta       : angle between surface normal and the flow [rad]
%       matID       : vector of material identifiers for each mesh element
%       eqmodel     : GSI model
%       param_eq    : structure of GSI model parameter inputs
%
% Outputs:
%       cp     : normal pressure coefficient;
%       ctau   : shear stress coefficient;
%       cd     : drag coefficient;
%       cl     : lift coefficient
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
%
%------------- BEGIN CODE --------------

% Check material references
nMat = max(matID); % Number of materials defined in .obj
if nMat == 0
    % if zero material references
    matID(matID == 0) = 1;
    nMat = 1;
end

%    warning('Material ID not correctly defined in .obj file for multiple surface types')
    
% Apply correct material properties per mesh face
if isfield(param_eq,'alpha')
    nParam = numel(param_eq.alpha);
    if ~isequal(nMat,nParam)
        error('Number of materials defined in .obj does not match number of materials defined in GSI parameter input')
    end
    param_eq.alpha = param_eq.alpha(matID);
end

if isfield(param_eq,'alphaN')
    nParam = length(param_eq.alphaN);
    if ~isequal(nMat,nParam)
        error('Number of materials defined in .obj does not match number of materials defined in GSI parameter input')
    end
    param_eq.alphaN = param_eq.alphaN(matID);
end
    
if isfield(param_eq,'sigmaN')
    nParam = length(param_eq.sigmaN);
    if ~isequal(nMat,nParam)
        error('Number of materials defined in .obj does not match number of materials defined in GSI parameter input')
    end
    param_eq.sigmaN = param_eq.sigmaN(matID);
end

if isfield(param_eq,'sigmaT')
    nParam = length(param_eq.sigmaT);
    if ~isequal(nMat,nParam)
        error('Number of materials defined in .obj does not match number of materials defined in GSI parameter input')
    end
    param_eq.sigmaT = param_eq.sigmaT(matID);
end

% Calculate GSIM
gsimodel = str2func(['coeff_',param_eq.gsi_model]);
[cp, ctau, cd, cl] = gsimodel(param_eq, delta);

%------------- END OF CODE --------------
