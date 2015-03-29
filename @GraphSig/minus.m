function G = minus(G1,G2)
% Take the difference of the data from G1 and G2
%
% Input
%   G1      GraphSig object #1
%   G2      GraphSig object #2
%
% Output
%   G       a GraphSig object with data vector G1.f-G2.f
%
%
%
% Copyright 2015 The Regents of the University of California
%
% Implemented by Jeff Irion (Adviser: Dr. Naoki Saito)



G = G1;
G.f = G1.f-G2.f;

end