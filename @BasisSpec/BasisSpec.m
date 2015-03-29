function BS = BasisSpec(levlist,levlengths,c2f,description)
% Constructor for a BasisSpec object
%
% Input
%   levlist         the integer sequence that specifies a particular basis
%   levlengths      the integer sequence that specifies the length of each
%                   basis block in "levlist" (optional)
%   c2f             if true (default), this indicates that the basis comes
%                   from the coarse-to-fine dictionary; if false, this
%                   indicates that the basis comes from the fine-to-coarse
%                   dictionary
%   description     a description of the specified basis
%
% Output
%   BS              a struct with four fields: levlist, levlengths, c2f,
%                   and description
%
%
%
% Copyright 2015 The Regents of the University of California
%
% Implemented by Jeff Irion (Adviser: Dr. Naoki Saito)



%% if there is no meaningful rs data
if ~exist('levlist','var')
    datastruct.levlist = [];
else
    datastruct.levlist = levlist;
end


%% if there is no meaningful ind data
if ~exist('levlengths','var') || length(datastruct.levlist) ~= length(levlengths)
    datastruct.levlengths = [];
else
    datastruct.levlengths = levlengths;
end


%% if c2f is not specified
if ~exist('c2f','var') || ~islogical(c2f)
    datastruct.c2f = true;
else
    datastruct.c2f = c2f;
end


%% if c2f is not specified
if ~exist('description','var') || ~ischar(description)
    datastruct.description = [];
else
    datastruct.description = description;
end



BS = class(datastruct,'BasisSpec');


end