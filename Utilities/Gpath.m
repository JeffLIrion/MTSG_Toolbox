function G = Gpath(N,f,name)
% Generate a GraphSig object for a 1-D path of length N
%
% Input
%   N       the length of the path
%   f       the signal to be used
%   name    the name for the GraphSig object
%
% Output
%   G       the GraphSig ojbect
%
%
%
% Copyright 2015 The Regents of the University of California
%
% Implemented by Jeff Irion (Adviser: Dr. Naoki Saito)



% the x coordinates
xy = (1:N)';

% the weight matrix
e = ones(N,1);
W = spdiags([ e 0*e e], [-1 0 1],N,N);

% the signal on the graph
if ~exist('f','var')
    f = ones(N,1);
elseif isscalar(f)
    f = f*ones(N,1);
end

% the name of the graph
if ~exist('name','var')
    name = sprintf('Path of length %d',N);
end

% create the GraphSig object
G = GraphSig(W,xy,f,name);


end