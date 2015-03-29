function [GGauss,e] = Adj2Gauss(G,e)
% Given a GraphSig object 'G' with a binary weight (adjacency) matrix,
% generate a GraphSig object 'GGauss' with a Gaussian weight matrix
%
% Input
%   G           a GraphSig object
%   e           (optional) the Gaussian parameter: 
%                   exp(-(dist(v_i,v_j)/e)^2)
%
% Output
%   GGauss      a GraphSig object with Gaussian weights
%   e           the Gaussian parameter: exp(-(dist(v_i,v_j)/e)^2)
%
%
% 
% Copyright 2015 The Regents of the University of California
%
% Implemented by Jeff Irion (Adviser: Dr. Naoki Saito)



% find the edges in the graph
[rows,cols] = find(G.W);

% compute e, if necessary
if nargin == 1
    e = median(pdist(G.xy));
end

% compute the new W
for j = 1:length(rows)
    G.W(rows(j),cols(j)) = exp(-(norm(G.xy(rows(j),:)-G.xy(cols(j),:),2)/e)^2);
end

GGauss = GraphSig(G.W, G.xy, G.f, sprintf('%s (Gaussian weight matrix with e = %f)',G.name,e), G.plotspecs);


end