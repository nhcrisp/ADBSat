% ADBSat instalation file.
% Adds the relative paths of the toolbox to the MATLAB path.
%
% Inputs:
%      '' (empty) : installs the toolbox
%      'uninstall': uninstalls the toolbox
%
% Outputs:
%       None.
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
%--------- BEGIN CODE -----------%

function ADBSat_install(input)

if nargin == 0
    install = 1;
elseif nargin > 0
    if strcmpi(input,'install')
        install = 1;
    elseif strcmpi(input,'uninstall')
        install = 0;
    else
        warning('Install function input not recognised')
        return
    end
end

%Get full path to this file as a string
currentpath = string(mfilename('fullpath'));
%Find folder separators
pathparts = split(currentpath, filesep);
%ADBSat base folder
ADBSat_path = char(join(pathparts(1:end-2),filesep));

if install % Install toolbox
    
    p = genpath(ADBSat_path);
    addpath(p)
    savepath
    disp(['ADBSat_path = ', ADBSat_path])
    
    % Create folders for examples
    if not(isfolder([ADBSat_path,filesep,'inou',filesep,'stl_files']))
        mkdir([ADBSat_path,filesep,'inou',filesep,'stl_files'])
    end
    if not(isfolder([ADBSat_path,filesep,'inou',filesep,'models']))
        mkdir([ADBSat_path,filesep,'inou',filesep,'models'])
    end
    if not(isfolder([ADBSat_path,filesep,'inou',filesep,'results']))
        mkdir([ADBSat_path,filesep,'inou',filesep,'results'])
    end
    
    disp('ADBSat installed successfully')
  
elseif ~install % Uninstall toolbox
   
    % Dynamic path file
    ADBSat_path = ADBSat_dynpath;
    
    p = genpath(ADBSat_path);
    rmpath(p)
    savepath
    
    disp('ADBSat uninstalled successfully')
    
end
%------- END CODE --------%
