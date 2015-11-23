%  Initialize Matlab's path to include MTSG_Toolbox
%
%
%
% Copyright 2015 The Regents of the University of California
%
% Implemented by Jeff Irion (Adviser: Dr. Naoki Saito)



% version info
MTSGVERSION = 001;


% declare MTSGPATH
MTSGPATH = [fileparts(mfilename('fullpath')),filesep];


% generate the path string ==> includes the folder and multiple levels of subfolders below it
V = version('-date');
V = str2double(V(end-3:end));

% make sure the Matlab version is 6.x or above
if V < 2006
    fprintf('Warning: This version is only supported on Matlab 6.x or above');
    mtsgp=genpath(MTSGPATH,1);
else
    mtsgp=genpath(MTSGPATH);
end

addpath(mtsgp);


% generate datasets
Acquire_Data;

% print a "success" message
fprintf('\n\n\nMTSG Toolbox v %g was successfully installed.\n\n\n', MTSGVERSION);


% clear variables
clear MTSGVERSION MTSGPATH V mtsgp