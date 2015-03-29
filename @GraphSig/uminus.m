function G = uminus(G)
% Negate the signal
%
% Input
%   G       a GraphSig object
%
% Output
%   G       a GraphSig object with data vector -G.f
%
%
%
% Copyright 2015 The Regents of the University of California
%
% Implemented by Jeff Irion (Adviser: Dr. Naoki Saito)



G.f = -G.f;

end