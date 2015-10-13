function dmatrix = dvec2dmatrix(dvec,GP,BS)
% Given a vector of expansion coefficients, convert it to a matrix.
%
% Inputs
%   dvec        a vector of expansion coefficients
%   GP          a GraphPart object
%   BS          a BasisSpec object
%
% Outputs
%   dmatrix     a matrix of expansion coefficients
%
%
%
% Copyright 2015 The Regents of the University of California
%
% Implemented by Jeff Irion (Adviser: Dr. Naoki Saito)



%% 0. Preliminaries

% extract data
[levlist,levlengths] = ExtractData(BS,GP);

% constants
[N,jmax] = size(GP.rs);
N = N-1;
[~,fcols] = size(dvec);

% allocate space
dmatrix = zeros(N,jmax,fcols);


%% 1. Put the entries in the vector into the correct places in the matrix
n = 1;
for row = 1:length(levlist)
    dmatrix(n:n+levlengths(row)-1,levlist(row),:) = dvec(n:n+levlengths(row)-1,:);
    n = n+levlengths(row);
end


end