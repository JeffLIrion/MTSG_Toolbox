function y = norm(G,p)
% Take the norm of G
%
% Input
%   G       a GraphSig object
%   p       the type of norm (default is the 2-norm)
%
% Output
%   y       the norm of G
%
%
%
% Copyright 2015 The Regents of the University of California
%
% Implemented by Jeff Irion (Adviser: Dr. Naoki Saito)



if ~exist('p','var')
    y = norm(G.f);
else
    y = norm(G.f,p);
end


end