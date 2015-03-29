function BS = HaarBasisSpec(GP)
% Specify the Haar basis for a given graph partitioning
%
% Input
%   GP      a GraphPart object
% 
% Output
%   BS      a BasisSpec object corresponding to the Haar basis
%
%
%
% Copyright 2015 The Regents of the University of California
%
% Implemented by Jeff Irion (Adviser: Dr. Naoki Saito)



% determine jmax
[~,jmax] = size(GP.rs);

% allocate space for levlist
levlist = zeros(jmax,1,'uint8');

% fill in levlist for the Haar basis
levlist(1) = jmax;
for row = 2:jmax
    levlist(row) = jmax+2-row;
end

% make a BasisSpec object
BS = BasisSpec(levlist,[],false,'Haar basis');

% fill in the levlengths field of the BasisSpec object
BS = levlist2levlengths(GP,BS);


end