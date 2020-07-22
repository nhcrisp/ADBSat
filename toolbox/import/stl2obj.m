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