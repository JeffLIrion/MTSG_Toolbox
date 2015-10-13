function [levlist,levlengths,c2f,description] = ExtractData(BS,GP)
% Returns data in a GraphSig object
%
% Input
%   BS              the input BasisSpec object whose data is to be 
%                   extracted
%   GP              a GraphPart object for filling in levlengths, if needed
%
% Output
%   levlist         the integer sequence that specifies a particular basis
%   levlengths      the integer sequence that specifies the length of each
%                   basis block in "levlist" (optional)
%   c2f             if true (default), this indicates that the partition
%                   refers to the coarse-to-fine dictionary; if false, this
%                   indicates that the partition refers to the fine-to-
%                   coarse dictionary
%   description     a description of the specified basis
%
%
%
% Copyright 2015 The Regents of the University of California
%
% Implemented by Jeff Irion (Adviser: Dr. Naoki Saito)



% fill in levlengths if it is missing and if GP is provided
if isempty(BS.levlengths) && exist('GP','var')
    BS = levlist2levlengths(GP,BS);
end

levlist = BS.levlist;
levlengths = BS.levlengths;
c2f = BS.c2f;
description = BS.description;


end