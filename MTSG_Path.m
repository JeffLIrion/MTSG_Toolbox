%  MTSG_Path -- initialize Matlab's path to include MTSG_Toolbox
%
%
%
% Copyright 2015 The Regents of the University of California
%
% Implemented by Jeff Irion (Adviser: Dr. Naoki Saito)



%% version info

MTSGVERSION = 001;
fprintf('\nWelcome to the MTSG Toolbox v %g\n\n', MTSGVERSION);


%% identify the type of computer and declare PATHNAMESEPARATOR, MTSGPATH, and MATLABPATHSEPARATOR accordingly

Friend = computer;
if strcmp(Friend,'MAC2')
    PATHNAMESEPARATOR = ':';
    MTSGPATH = ['Macintosh HD:Build:MTSG_Toolbox', num2str(MTSGVERSION), PATHNAMESEPARATOR];
    MATLABPATHSEPARATOR = ';';
elseif isunix
    PATHNAMESEPARATOR = '/';
    MTSGPATH = [pwd, PATHNAMESEPARATOR];
    MATLABPATHSEPARATOR = ':';
elseif strcmp(Friend(1:2),'PC')
    PATHNAMESEPARATOR = '\';
    MTSGPATH = [pwd, PATHNAMESEPARATOR];
    MATLABPATHSEPARATOR = ';';
else
    fprintf('\n\nI don''t recognize this computer.\nPathnames not set.\n\nSolution: edit MTSG_Path.m\n\n');
end


%% generate the path string ==> includes the folder and multiple levels of subfolders below it

V = version('-date');
V = str2double(V(end-3:end));

if V < 2006
    fprintf('Warning: This version is only supported on Matlab 6.x or above');
    mtsgp=genpath(MTSGPATH,1);
else
    mtsgp=genpath(MTSGPATH);
end

addpath(mtsgp);