function ADBSat_install(input)

%  ADBSat instalation
%
%    Adds the relative paths of the toolbox to the MATLAB path.
%
%    input:
%      '' (empty) : installs the toolbox
%      'uninstall': uninstalls the toolbox
%
%*********************************************************************


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

% Absolute path

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
    
    disp(['ADBS_path = ', ADBSat_path])
    
    disp('ADBSat installed successfully')
  
elseif ~install % Uninstall toolbox
   
    % Dynamic path file
    ADBSat_path = ADBSat_dynpath;
    
    p = genpath(ADBSat_path);
    rmpath(p)
    savepath
    
    disp('ADBSat uninstalled successfully')
    
end