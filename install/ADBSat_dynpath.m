function [ADBSat_path] = ADBSat_dynpath( )
%Dynamic Path to ADBSat Base Folder
% Must not be moved from relative folder

%Get full path to this file as a string
currentpath = (mfilename('fullpath'));
%Find folder separators
pathparts = split(currentpath, filesep);
%ADBSat base folder
ADBSat_path = char(join(pathparts(1:end-2),filesep));