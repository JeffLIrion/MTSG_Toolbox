function GInvEuc = Adj2InvEuc(G)
% Given a GraphSig object 'G' with a binary weight (adjacency) matrix, 
% generate a GraphSig object 'GInvEuc' with an inverse Euclidean distance 
% matrix
%
% Input
%   G           a GraphSig object
%
% Output
%   GInvEuc     a GraphSig object with inverse Euclidean weights
%
%
% 
% Copyright 2015 The Regents of the University of California
%
% Implemented by Jeff Irion (Adviser: Dr. Naoki Saito)



% find the edges
[rows,cols] = find(G.W);

% remove duplicate points
for j = length(rows):-1:1
    if j <= length(rows) && norm(G.xy(rows(j),:)-G.xy(cols(j),:),2) < 10^3*eps
        G.W(rows(j),:) = [];
        G.W(:,rows(j)) = [];
        G.xy(rows(j),:) = [];
        G.f(rows(j))   = [];
        [rows,cols] = find(G.W);
    end
end

% compute the new W
for j = 1:length(rows)
    G.W(rows(j),cols(j)) = 1/norm(G.xy(rows(j),:)-G.xy(cols(j),:),2);
end

GInvEuc = GraphSig(G.W, G.xy, G.f, strcat(G.name,' (inverse Euclidean weight matrix)'), G.plotspecs);


end