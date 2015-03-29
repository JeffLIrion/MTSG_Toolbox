function G = power(G1,G2)
% Perform an element-wise power operation
%
% Input
%   G1      GraphSig object #1
%   G2      GraphSig object #2 OR a constant
%
% Output
%   G       a GraphSig object with data vector G1.f .^ G2 or G1.f .^ G2.f
%
%
%
% Copyright 2015 The Regents of the University of California
%
% Implemented by Jeff Irion (Adviser: Dr. Naoki Saito)



G = G1;

if ~strcmp(class(G2),'GraphSig')
    G.f = G1.f .^ G2;
else
    G.f = G1.f .^ G2.f;
end


end