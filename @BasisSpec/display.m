function display(BS)
% Display info about the BasisSpec ojbect
%
% Input
%   BS      a BasisSpec object
%
%
%
% Copyright 2015 The Regents of the University of California
%
% Implemented by Jeff Irion (Adviser: Dr. Naoki Saito)



%% Input name
fprintf('\n\nBasisSpec object ''%s''\n\n\n',inputname(1));


%% levlist
if ~isempty(BS.levlist)
    levlist = sprintf('|    "levlist"?  Yes (length = %d)',length(BS.levlist));
else
    levlist = sprintf('|    "levlist"?  No');
end


%% levlengths
if ~isempty(BS.levlengths)
    levlengths = sprintf('| "levlengths"?  Yes');
else
    levlengths = sprintf('| "levlengths"?  No');
end


%% c2f
if BS.c2f
    dictionary = sprintf('|   Dictionary:  coarse-to-fine');
else
    dictionary = sprintf('|   Dictionary:  fine-to-coarse');
end


%% description
if ~isempty(BS.description)
    description = sprintf('|  Description:  %s',BS.description);
else
    description = [];
end


%% print the info
tablewidth = max([length(levlist), length(levlengths), length(dictionary), length(description)])+2;

fprintf(' %s \n',repmat('-',1,tablewidth-2));
fprintf('|%s|\n',blanks(tablewidth-2));
fprintf('%s%s|\n',levlist,blanks(tablewidth-length(levlist)-1));
fprintf('%s%s|\n',levlengths,blanks(tablewidth-length(levlengths)-1));
fprintf('|%s|\n',blanks(tablewidth-2));
fprintf('%s%s|\n',dictionary,blanks(tablewidth-length(dictionary)-1));
fprintf('|%s|\n',blanks(tablewidth-2));
if ~isempty(description)
    fprintf('%s%s|\n',description,blanks(tablewidth-length(description)-1));
    fprintf('|%s|\n',blanks(tablewidth-2));
end
fprintf(' %s ',repmat('-',1,tablewidth-2));
    
fprintf('\n\n\n');


end