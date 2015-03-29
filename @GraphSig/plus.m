function G = plus(G1,G2)
% Take the sum of G1 and G2
%
% Input
%   G1      GraphSig object #1
%   G2      GraphSig object #2
%
% Output
%   G       a GraphSig object with data vector G1.f + G2.f
%
%
%
% Copyright 2015 The Regents of the University of California
%
% Implemented by Jeff Irion (Adviser: Dr. Naoki Saito)



if ~isa(G2,'GraphSig')
    G = G1;
    G.f = G.f + G2;
elseif ~isa(G1,'GraphSig')
    G = G2;
    G.f = G.f + G1;
else
    G = G1;
    G.f = G1.f+G2.f;
end


end