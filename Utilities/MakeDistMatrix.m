function D = MakeDistMatrix(xy,distfun)
% Generate a pairwise distance matrix
%
% Input
%   xy           the N-by-d matrix representing N points in dimension d
%   distfun     the distance metric
%
% Output
%   D           the N-by-N pairwise distance matrix
%
%
%
% Copyright 2015 The Regents of the University of California
%
% Implemented by Jeff Irion (Adviser: Dr. Naoki Saito)



if nargin < 2
    D = squareform(pdist(xy));
else
    D = squareform(pdist(xy,distfun));
end


end