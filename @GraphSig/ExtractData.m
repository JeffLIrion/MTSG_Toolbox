function [W,xy,f,name,plotspecs] = ExtractData(G)
% Return the data in a GraphSig object.
%
% Input
%   G           the input GraphSig object whose data is to be extracted
%
% Output
%   W           the edge weight matrix in 'G'
%   xy          the coordinates of the vertices (not necessarily 2D) in 'G'
%   f           the values of the function/data at the vertices in 'G'
%   name        the name of the GraphSig object
%   plotspecs   the plotspecs of the GraphSig object
%
%
%
% Copyright 2015 The Regents of the University of California
%
% Implemented by Jeff Irion (Adviser: Dr. Naoki Saito)



W = G.W;
xy = G.xy;
f = G.f;
name = G.name;
plotspecs = G.plotspecs;


end