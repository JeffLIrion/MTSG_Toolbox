function G = mrdivide(G1,G2)
% Divide a (graph signal/constant) by a (graph signal/constant)
%
% Input
%   G1      GraphSig object #1 OR a constant
%   G2      GraphSig object #2 OR a constant
%
% Output
%   G       a GraphSig object with data vector G1.f / G2.f
%
%
%
% Copyright 2015 The Regents of the University of California
%
% Implemented by Jeff Irion (Adviser: Dr. Naoki Saito)



if strcmp(class(G1),'GraphSig') && ~strcmp(class(G2),'GraphSig')
    G = G1;
    G.f = G1.f / G2;
elseif strcmp(class(G2),'GraphSig') && ~strcmp(class(G1),'GraphSig')
    G = G2;
    G.f = G1 * G2.f .^(-1);
else
    G = G1;
    G.f = G1.f ./ G2.f;
end


end