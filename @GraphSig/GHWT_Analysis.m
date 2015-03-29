function [dmatrix,GP] = GHWT_Analysis(G,GP)
% For a GraphSig object 'G', generate the matrix of GHWT expansion 
% coefficients
%
% Input
%   G               a GraphSig object
%   GP              a GraphPart object
%
% Output
%   dmatrix         the matrix of expansion coefficients
%   GP              a GraphPart object
%
%
%
% Copyright 2015 The Regents of the University of California
%
% Implemented by Jeff Irion (Adviser: Dr. Naoki Saito)



%% 0. Preliminaries

if ~exist('GP','var')
    GP = PartitionTreeFiedler(G);
end

[ind,rs] = ExtractData(GP);

[N,jmax] = size(rs);
N = N-1;

% allocate space for the expansion coefficients
[~,fcols] = size(G.f);
dmatrix = zeros(N,jmax,fcols);
dmatrix(:,jmax,:) = G.f(ind,:);


%% 1. Perform the transform

% generate expansion coefficients, tag, and compinfo
[GP,dmatrix] = GHWT_Core(GP,dmatrix);

% fill in the remaining (i.e. fine-to-coarse) GraphPart fields
if nargout > 1
    GP = FineToCoarse(GP);
end



end