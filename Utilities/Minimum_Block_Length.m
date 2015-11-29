function dmatrix = Minimum_Block_Length(dmatrix,GP,Nmin,value)
% Set all coefficients in blocks that are less than 'Nmin' in length 
% (coarse-to-fine dictionary) to 'value' (default = Inf).
%
% Input
%   dmatrix     the matrix of expansion coefficients
%   GP          a GraphPart object
%   Nmin        entries in blocks of length < Nmin will be set to value
%   value       the value to which the entries will be set (default = Inf)
%
% Output
%   dmatrix     the coefficients with entries in blocks of length < Nmin
%               set to Inf
%
%
%
% Copyright 2015 The Regents of the University of California
%
% Implemented by Jeff Irion (Adviser: Dr. Naoki Saito)



[~,jmax,~] = size(dmatrix);

if ~exist('value','var')
    value = Inf;
end

for j = 1:jmax
    BS = LevelBasisSpec(GP,j-1);
    [~,levlengthsfull] = BSfull(GP,BS);
    dmatrix(levlengthsfull < Nmin,j,:) = value;
end


end