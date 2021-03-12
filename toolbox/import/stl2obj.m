function [status] = stl2obj(filename)
%STL2OBJ converts a .STL file to a .OBJ file using meshlabserver
%
% Inputs:
%   filename: mesh filename
%
% Outputs:
%   status  : error flag [1/0]
%
% Author: Nicholas H. Crisp
% The University of Manchester
% February 2019
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

[ADBSat_path] = ADBSat_dynpath;

filein = fullfile(ADBSat_path,'inou','stl_files',[filename,'.stl']);
fileout = fullfile(ADBSat_path,'inou','obj_files',[filename,'.obj']);

script = fullfile(ADBSat_path,'toolbox','import','meshlab_reset_origin.mlx');

% System Call
[status] = system(['meshlabserver -i ', filein,' -o ', fileout, ' -s ', script]);

% Shell Escape
%!meshlabserver -i fileout -o fileout -s meshlab_reset_origin.mlx

if status ~= 0
    error('.OBJ file not correctly created. Check input file and that meshlabserver is available on the system PATH.')
end

%------------- END OF CODE --------------
