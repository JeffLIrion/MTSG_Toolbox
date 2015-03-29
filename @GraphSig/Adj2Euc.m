function GEuc = Adj2Euc(G)
% Given a GraphSig object 'G' with a binary weight (adjacency) matrix, 
% generate a GraphSig object 'GEuc' with a Euclidean distance matrix
%
% Input
%   G           a GraphSig object
%
% Output
%   GEuc        a GraphSig object with Euclidean weights
%
% 
% 
% Copyright 2015 The Regents of the University of California
%
% Implemented by Jeff Irion (Adviser: Dr. Naoki Saito)



% find the edges in the graph
[rows,cols] = find(G.W);

% compute the new W
for j = 1:length(rows)
    G.W(rows(j),cols(j)) = norm(G.xy(rows(j),:)-G.xy(cols(j),:),2);
end

GEuc = GraphSig(G.W, G.xy, G.f, strcat(G.name,' (Euclidean weight matrix)'), G.plotspecs);


end