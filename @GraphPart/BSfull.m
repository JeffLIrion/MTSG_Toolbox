function [levlistfull,levlengthsfull,transfull] = BSfull(GP,BS,trans)
% Given a BasisSpec object, return the full-length, redundant levlist, 
% levlengths, and trans descriptions.  
%
% Input
%   GP              a GraphPart object
%   BS              a BasisSpec object
%   trans           a specification of the transforms used for the
%                   HGLET-GHWT hybrid transform
%
% Output
%   levlistfull     the full-length, redundant levels list description
%   levlengthsfull  the full-length, redundant levels lengths description
%   transfull       the full-length, redundant trans description
%
%
%
% Copyright 2015 The Regents of the University of California
%
% Implemented by Jeff Irion (Adviser: Dr. Naoki Saito)



%% 0. Preliminaries

% extract data
[levlist,levlengths] = ExtractData(BS);
if isempty(levlengths)
    BS = levlist2levlengths(GP,BS);
    [levlist,levlengths] = ExtractData(BS);
end

% allocate space
N = length(GP.ind);
levlistfull = zeros(N,1,'uint8');

if nargout > 1
    levlengthsfull = zeros(N,1,class(GP.ind));
end

if nargout == 3 && exist('trans','var')
    [~,cols] = size(trans);
    transfull = false(N,cols);
end


%% 1. Fill out the redundant descriptions
ind = 0;
for row = 1:length(levlist)
    levlistfull(ind+1:ind+levlengths(row)) = levlist(row);
    
    if nargout > 1
        levlengthsfull(ind+1:ind+levlengths(row)) = levlengths(row);
    end
    
    if exist('transfull','var')
        transfull(ind+1:ind+levlengths(row),:) = repmat(trans(row,:),levlengths(row),1);
    end
    
    ind = ind+levlengths(row);
end


end