function BS = DWTBasisSpec(GP,j0)
% Specify the basis that corresponds to the Haar discrete wavelet transform
% (DWT) performed up to level j0
%
% Input
%   GP      a GraphPart object
%   j0      the level to which the DWT is carried out; j0 >= 0
% 
% Output
%   BS      a BasisSpec object corresponding to the basis corresponding to
%           the Haar DWT performed up to level j0
%
%
%
% Copyright 2015 The Regents of the University of California
%
% Implemented by Jeff Irion (Adviser: Dr. Naoki Saito)



% extract data
[~,jmax] = size(GP.rs);

% make sure the specified j0 value is appropriate
if ~exist('j0','var') || j0 < 0 || j0 >= jmax
    j0 = jmax-1;
end

% allocate space for levlist
levlist = zeros(j0+1,1,'uint8');

% fill in levlist for the DWT basis
levlist(1) = j0+1;
for row = 2:j0+1
    levlist(row) = j0+3-row;
end

% make a BasisSpec object
BS = BasisSpec(levlist,[],false,sprintf('DWT level j=%d basis',j0));

% fill in the levlengths field of the BasisSpec object
BS = levlist2levlengths(GP,BS);


end