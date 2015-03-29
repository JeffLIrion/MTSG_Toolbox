function y = dot(G1,G2)
% Take the dot product of G1 and G2
%
% Input
%   G1      GraphSig object #1
%   G2      GraphSig object #2
%
% Output
%   y       the dot product of G1 and G2
%
%
% 
% Copyright 2015 The Regents of the University of California
%
% Implemented by Jeff Irion (Adviser: Dr. Naoki Saito)



y = full(dot(G1.f,G2.f));


end